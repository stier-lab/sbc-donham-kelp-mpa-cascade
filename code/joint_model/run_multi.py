"""Fit and simulate the multi-taxon model (stage 2, experimental, review only).

Examples (from the repository root):
  python3 code/joint_model/run_multi.py fit --curve rw --cpus 4 --warmup 3000 --samples 6000 --thin 4
  python3 code/joint_model/run_multi.py fit --no-paths                      # taxa side by side, no links
  python3 code/joint_model/run_multi.py simulate --scenario paths --seed 1  # known path strengths, real design
  python3 code/joint_model/run_multi.py fit --taxa lob,pur                  # two taxa only
"""
import argparse, json, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).parent))
import run as single                      # sets the host device count and 64-bit mode before JAX starts
from run import write_retry, attach_counts, TRUTHS, CPUS, USES
import jax
import numpyro
from numpyro.infer import MCMC, NUTS, init_to_uniform
from data_multi import build_multi, TAXA, EDGES, CONFIG
from model_multi import model, basis_row

SITES = ["E10_direct", "E10_total", "g_knots", "rho", "sigma_process", "sigma_region_year", "sigma_drift", "sigma_drift_region",
         "sigma_visit", "sigma_b", "inv_sqrt_phi", "tau_D", "tau_L", "sigma_shock", "mu_lambda", "sigma_lambda", "sigma_curve"]

# Simulation truths. Direct response increments per knot interval (0-4, 4-8, 8-12, 12-16, 16-20 years).
SIM_TRUTH = {"lob": "lobster", "shp": "lobster", "pur": "urchin", "red": "urchin", "kelp": "urchin"}
SIM_DIRECT = {"lob": [1.0, 0.7, 0.45, 0.4, 0.4], "shp": [0.1, 0.1, 0.05, 0.0, 0.0], "pur": [-0.15, -0.35, -0.30, -0.10, 0.0],
              "red": [-0.10, -0.10, -0.10, -0.10, -0.10], "kelp": [0.3, 0.3, 0.2, 0.0, 0.0]}
SIM_BETA = {"paths": {("pur", "lob"): -0.30, ("pur", "shp"): -0.20, ("red", "lob"): -0.20, ("red", "shp"): -0.10,
                      ("kelp", "pur"): -0.30, ("kelp", "red"): -0.15},
            "nopaths": {},
            # no paths, but a local driver shared by taxa at a station (tests whether a common cause fakes a path)
            "confounded": {}}
SIM_DRIFT = dict(sigma_drift=0.10, sigma_drift_region=0.10)


def simulate(data, edges, scenario, seed, lag=1, direct=True, common_sd=0.0):
    """Independent NumPy simulator of the linked model on the real survey designs.
    common_sd > 0 adds a station-level AR(1) driver shared by every taxon (enters each taxon's log density with a
    loading 1); the fitted model has no such term."""
    rng = np.random.default_rng(seed)
    betas = SIM_BETA[scenario]; dev, truth = {}, {}
    all_st = sorted({s for D in data.values() for s in D["stations"]})
    y0 = min(int(D["years"][0]) for D in data.values()); y1 = max(int(D["years"][-1]) for D in data.values())
    common = np.zeros((y1 - y0 + 1, len(all_st)))
    if common_sd > 0:
        r = 0.9; common[0] = rng.normal(0, common_sd, len(all_st))
        for t in range(1, len(common)):
            common[t] = r * common[t - 1] + common_sd * np.sqrt(1 - r ** 2) * rng.normal(size=len(all_st))
    c_ix = {s: i for i, s in enumerate(all_st)}
    ys = {}
    for g, D in data.items():
        tr = TRUTHS[SIM_TRUTH[g]]
        S, M, R, T = len(D["stations"]), len(D["mpas"]), len(D["regions"]), len(D["years"])
        inc = np.array(SIM_DIRECT[g], dtype=float) * (1.0 if direct else 0.0)
        lam = np.where(D["mpa_ba"] == 1, tr["mu_lambda"] + tr["sigma_lambda"] * rng.normal(size=M), rng.normal(0, 1.0, size=M))
        b = D["log_mean_density"] + tr["level"] + tr["sigma_b"] * rng.normal(size=S) + D["st_prot"] * lam[D["st_mpa"]]
        gam = tr["sigma_region_year"] * rng.normal(size=(R, T))
        a = np.cumsum(SIM_DRIFT["sigma_drift"] * rng.normal(size=T)); gam = gam + (a - a.mean())[None, :]
        br = np.cumsum(SIM_DRIFT["sigma_drift_region"] * rng.normal(size=(R, T)), axis=1); gam = gam + br - br.mean(1, keepdims=True)
        delta = np.zeros((T, S)); delta[0] = rng.normal(0, tr["sigma_process"] / np.sqrt(1 - tr["rho"] ** 2), size=S)
        for t in range(1, T):
            delta[t] = tr["rho"] * delta[t - 1] + tr["sigma_process"] * rng.normal(size=S)
        Dm = tr["tau_D"] * rng.normal(size=M); Lm = tr["tau_L"] * rng.normal(size=M); omega = tr["sigma_shock"] * rng.normal(size=T)
        m = D["st_mpa"]
        F = D["grid_after"] * (D["grid_B"] @ inc + Lm[m][None, :] + Dm[m][None, :] * D["grid_age"] / 10.0 + omega[:, None])
        path = np.zeros((T, S))
        for k, j, lk in edges:
            if k == g:
                path += betas.get((k, j), 0.0) * lk["valid"] * dev[j][lk["src_t"]][:, lk["src_s"]]
        cm = common[int(D["years"][0]) - y0: int(D["years"][0]) - y0 + T][:, [c_ix[s] for s in D["stations"]]]
        d = b[None, :] + delta + F + path + cm
        n = D["span"].sum(0)
        dev[g] = D["span"] * (d - (D["span"] * d).sum(0) / np.where(n > 0, n, 1.0))
        nu = tr["sigma_visit"] * rng.normal(size=len(D["visit_prog"]))
        eta = d[D["t"], D["s"]] + gam[D["st_region"][D["s"]], D["t"]] + nu[D["v"]] + np.array(tr["zeta"])[D["zone"]]
        mu = np.exp(eta + D["log_area"])[D["unit_cell"]]; phi = 1.0 / tr["inv_sqrt_phi"] ** 2
        ys[g] = rng.negative_binomial(phi, phi / (phi + mu))
        truth[g] = dict(inc=inc.tolist(), zero_fraction=float((ys[g] == 0).mean()))
    # true direct and total effects at ten years
    def tot(g, age):
        return float(basis_row(age) @ np.array(truth[g]["inc"])
                     + sum(betas.get((k, j), 0.0) * tot(j, max(age - lag, 0.0)) for k, j, _ in edges if k == g))
    for g in data:
        truth[g].update(E10_direct=float(basis_row(10.0) @ np.array(truth[g]["inc"])), E10_total=tot(g, 10.0))
    truth["beta"] = [betas.get((k, j), 0.0) for k, j, _ in edges]
    truth["common_sd"] = common_sd
    return ys, truth


def fit(data, edges, ys, config, curve="rw", lag=1, beta_scale=0.5, paths=True, warmup=600, samples=600, chains=4, seed=0,
        depth=10, accept=0.99, step0=0.01, thin=1, **common):
    data = {g: attach_counts(dict(D), ys[g]) for g, D in data.items()}
    kw = dict(curve=curve, lag=lag, beta_scale=beta_scale, paths=paths, **common)
    kern = NUTS(model, init_strategy=init_to_uniform(radius=0.5), target_accept_prob=accept, max_tree_depth=depth, step_size=step0)
    mc = MCMC(kern, num_warmup=warmup, num_samples=samples, num_chains=chains, thinning=thin,
              chain_method="parallel" if CPUS > 1 else "sequential", progress_bar=False)
    t0 = time.time()
    mc.run(jax.random.PRNGKey(seed), data, edges, config, **kw,
           extra_fields=("diverging", "num_steps", "potential_energy", "adapt_state.step_size"))
    return mc, time.time() - t0


def keep_sites(mc):
    """Sites worth keeping from one MCMC run, as arrays shaped (chains, draws, ...), plus sampler diagnostics."""
    grouped = mc.get_samples(group_by_chain=True)
    ex = mc.get_extra_fields(group_by_chain=True)
    keep = {k: np.asarray(v) for k, v in grouped.items()
            if k == "beta" or k.split(":")[-1] in SITES or k.split(":")[-1].startswith("E10_via") or k.split(":")[-1] == "increments"}
    extra = dict(diverging=np.asarray(ex["diverging"]), num_steps=np.asarray(ex["num_steps"]),
                 potential_energy=np.asarray(ex["potential_energy"]), step_size=np.asarray(ex["adapt_state.step_size"])[:, -1])
    return keep, extra


def summarize(keep, ex, edges):
    from numpyro.diagnostics import summary
    s = summary({k: v for k, v in keep.items() if not k.endswith("increments")}, prob=0.9, group_by_chain=True)
    out = {}; worst, least = 1.0, np.inf; worst_site = None
    for k, v in s.items():
        out[k] = {f: np.asarray(v[f]).round(4).tolist() for f in ("mean", "5.0%", "95.0%", "n_eff", "r_hat")}
        r = float(np.nanmax(v["r_hat"]))
        if r > worst:
            worst, worst_site = r, k
        least = min(least, float(np.nanmin(v["n_eff"])))
    out["worst_rhat"] = round(worst, 4); out["worst_rhat_site"] = worst_site; out["min_neff"] = round(least, 1)
    out["divergences"] = int(ex["diverging"].sum())
    q = {}
    for k, v in keep.items():
        if k == "beta" or k.split(":")[-1].startswith("E10"):
            v2 = v.reshape(v.shape[0] * v.shape[1], -1)
            q[k] = dict(median=np.median(v2, 0).round(4).tolist(), q025=np.quantile(v2, .025, 0).round(4).tolist(),
                        q975=np.quantile(v2, .975, 0).round(4).tolist(), prob_negative=(v2 < 0).mean(0).round(4).tolist())
    out["intervals95"] = q
    out["edges"] = [f"{k}<-{j}" for k, j, _ in edges]
    pc = {k: v.reshape(v.shape[0], v.shape[1], -1).mean(1).round(3).tolist()
          for k, v in keep.items() if k == "beta" or k.split(":")[-1] in ("E10_total", "E10_direct", "sigma_process", "rho", "inv_sqrt_phi")}
    pc["potential_energy"] = ex["potential_energy"].mean(1).round(1).tolist()
    pc["step_size"] = ex["step_size"].round(5).tolist()
    pc["mean_leapfrog_steps"] = ex["num_steps"].mean(1).round(0).tolist()
    pc["divergences"] = ex["diverging"].sum(1).tolist()
    out["per_chain"] = pc
    return out


def run_chains(a, data, edges, ys, config, opts):
    """One process per chain. Running the five-taxon model with NumPyro's parallel chains in one process leaks memory
    (about 3.4 MB per second per chain on 2026-11-11: 20 GB in 25 minutes with four chains; one chain alone does not
    leak). So each chain is a separate single-chain process with its own seed; this process waits and merges."""
    import subprocess, tempfile
    tmp = Path(tempfile.mkdtemp(prefix="mpa_multi_"))
    argv = [x for x in sys.argv[1:]]
    for flag in ("--cpus", "--chains"):
        if flag in argv:
            i = argv.index(flag); del argv[i:i + 2]
    procs = [subprocess.Popen([sys.executable, __file__] + argv + ["--chains", "1", "--worker", str(k), "--worker-dir", str(tmp)])
             for k in range(a.chains)]
    codes = [p.wait() for p in procs]
    if any(codes):
        raise SystemExit(f"chain processes failed: {codes}")
    parts = [np.load(tmp / f"chain{k}.npz") for k in range(a.chains)]
    keep = {k: np.concatenate([p[k] for p in parts], 0) for k in parts[0].files if not k.startswith("__")}
    ex = {k[2:]: np.concatenate([p[k] for p in parts], 0) for k in parts[0].files if k.startswith("__")}
    secs = float(max(p["__seconds"][0] for p in parts)); ex.pop("seconds")
    for f in tmp.iterdir():
        f.unlink()
    tmp.rmdir()
    return keep, ex, secs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["simulate", "fit"])
    ap.add_argument("--taxa", default="lob,shp,pur,red,kelp", help="tags, sources before targets: " + ", ".join(f"{k}={v}" for k, v in TAXA.items()))
    ap.add_argument("--edges", default=None, help='links as "target<source,...", default: ' + ",".join(f"{k}<{j}" for k, j in EDGES))
    ap.add_argument("--no-paths", action="store_true", help="fit the taxa side by side with no links")
    ap.add_argument("--lag", type=int, default=1, choices=[0, 1, 2])
    ap.add_argument("--beta-scale", type=float, default=0.5, help="prior SD of the path coefficients")
    ap.add_argument("--curve", default="rw", choices=["piecewise", "rw"])
    ap.add_argument("--drift", default="regional", choices=["none", "global", "regional"])
    ap.add_argument("--inc-scale", type=float, default=2.0)
    ap.add_argument("--regions", default="three", choices=["seven", "three", "one"])
    ap.add_argument("--scenario", default="paths", choices=list(SIM_BETA))
    ap.add_argument("--no-direct", action="store_true", help="simulate: no direct protection effects")
    ap.add_argument("--common-sd", type=float, default=0.0, help="simulate: SD of a local driver shared by all taxa")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--warmup", type=int, default=600)
    ap.add_argument("--samples", type=int, default=600)
    ap.add_argument("--thin", type=int, default=1)
    ap.add_argument("--chains", type=int, default=4)
    ap.add_argument("--cpus", type=int, default=1)
    ap.add_argument("--depth", type=int, default=10)
    ap.add_argument("--step0", type=float, default=0.01)
    ap.add_argument("--accept", type=float, default=0.99)
    ap.add_argument("--omit-mpa", action="append", default=[])
    ap.add_argument("--save-draws", action="store_true")
    ap.add_argument("--one-process", action="store_true", help="NumPyro parallel chains in one process (leaks memory with five taxa)")
    ap.add_argument("--worker", type=int, default=None, help=argparse.SUPPRESS)
    ap.add_argument("--worker-dir", default=None, help=argparse.SUPPRESS)
    ap.add_argument("--name", default=None)
    ap.add_argument("--uses", default=USES)
    ap.add_argument("--out", default="outputs/joint_model_20261109")
    a = ap.parse_args()
    tags = a.taxa.split(",")
    ed = None if a.edges is None else [tuple(e.split("<")) for e in a.edges.split(",") if e]
    data, edges = build_multi(a.uses, tags, ed, lag=a.lag, regions=a.regions, omit_mpas=a.omit_mpa)
    for D in data.values():
        D.pop("cells")
    config = {g: CONFIG[g] for g in tags}
    opts = dict(curve=a.curve, lag=a.lag, beta_scale=a.beta_scale, paths=not a.no_paths, warmup=a.warmup, samples=a.samples,
                chains=a.chains, seed=a.seed, depth=a.depth, accept=a.accept, step0=a.step0, thin=a.thin, drift=a.drift, inc_scale=a.inc_scale)
    meta = dict(mode=a.mode, chain_processes=(a.chains > 1 and not a.one_process), taxa={g: TAXA[g] for g in tags}, regions=a.regions, omit_mpa=a.omit_mpa,
                links={f"{k}<-{j}": {x: lk[x] for x in ("n_target_station_years", "n_linked", "n_stations_linked")} for k, j, lk in edges},
                n_units={g: int(len(D["y"])) for g, D in data.items()}, stations={g: len(D["stations"]) for g, D in data.items()}, **opts)
    if a.mode == "simulate":
        name = a.name or f"msim_{a.scenario}_{a.curve}_seed{a.seed}"
        ys, truth = simulate(data, edges, a.scenario, a.seed, lag=a.lag, direct=not a.no_direct, common_sd=a.common_sd)
    else:
        name = a.name or f"mfit_{a.curve}{'_nopaths' if a.no_paths else ''}"
        ys = {g: D["y"] for g, D in data.items()}
    if a.worker is not None:          # one chain, own seed; hand the kept draws back to the parent process
        opts.update(chains=1, seed=a.seed * 1000 + a.worker)
        mc, secs = fit(data, edges, ys, config, **opts)
        keep, ex = keep_sites(mc)
        np.savez(Path(a.worker_dir) / f"chain{a.worker}.npz", **keep, **{f"__{k}": v for k, v in ex.items()}, __seconds=np.array([secs]))
        return
    if a.chains > 1 and not a.one_process:
        keep, ex, secs = run_chains(a, data, edges, ys, config, opts)
    else:
        mc, secs = fit(data, edges, ys, config, **opts)
        keep, ex = keep_sites(mc)
    res = summarize(keep, ex, edges)
    if a.mode == "simulate":
        res["truth"] = truth; res["scenario"] = a.scenario
    else:
        res["mpas"] = {g: D["mpas"] for g, D in data.items()}
        if a.save_draws:
            try:
                np.savez_compressed(Path(a.out) / f"draws_{name}.npz", **{k.replace(":", "__"): v for k, v in keep.items()
                                    if k == "beta" or k.split(":")[-1].startswith("E10") or k.split(":")[-1] in ("g_knots", "increments")})
            except OSError as e:
                print(f"draws not saved ({e})", flush=True)
    res.update(meta, seconds=round(secs))
    out = Path(a.out)
    try:
        out.mkdir(parents=True, exist_ok=True)
    except OSError:
        pass
    write_retry(out / f"{name}.json", json.dumps(res, indent=1))
    print(name, "seconds", round(secs), "div", res["divergences"], "worst_rhat", res["worst_rhat"], res["worst_rhat_site"],
          "min_neff", res["min_neff"], "beta", res.get("beta", {}).get("mean"))


if __name__ == "__main__":
    main()
