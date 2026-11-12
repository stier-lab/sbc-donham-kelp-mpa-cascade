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
from data_mass import build_mass, find_beside
from model_multi import model, basis_row

SITES = ["sigma_cluster", "sigma_trend", "mu_pretrend", "sigma_pretrend", "pretrend", "E10_mass", "E10_biomass", "m_g_knots", "m_rho", "m_sigma_process", "m_sigma_region_year", "m_sigma_drift", "m_sigma_drift_region",
         "m_sigma_station", "m_shape", "m_tau_D", "m_mu_lambda", "m_sigma_lambda", "m_sigma_curve", "m_zone", "m_theta",
         "E10_direct", "E10_total", "g_knots", "rho", "sigma_process", "sigma_region_year", "sigma_drift", "sigma_drift_region",
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
FLUCT = {}      # simulated station-level density fluctuation per taxon (for the coupled mass scenario)


def simulate(data, edges, scenario, seed, lag=1, direct=True, common_sd=0.0, station_trend_sd=0.0, cluster_trend_sd=0.0, pretrend=(0.0, 0.0)):
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
        # heterogeneity the recommended structure has no term for (tests of the v3 structure options)
        yrs = np.asarray(D["years"], dtype=float)
        st_tr = station_trend_sd * rng.normal(size=S)[None, :] * ((yrs - yrs.mean()) / 10.0)[:, None]
        chi = np.cumsum(cluster_trend_sd * rng.normal(size=(M, T)), axis=1); chi = (chi - chi.mean(1, keepdims=True))[D["st_mpa"]].T
        # differential linear trend at protected stations of before-and-after MPAs (mean, SD across MPAs; per decade)
        pre = (pretrend[0] + pretrend[1] * rng.normal(size=M))[D["st_mpa"]][None, :] * D["grid_rel"] if any(pretrend) else 0.0
        d = b[None, :] + delta + F + path + cm + st_tr + chi + pre
        n = D["span"].sum(0)
        dev[g] = D["span"] * (d - (D["span"] * d).sum(0) / np.where(n > 0, n, 1.0))
        lv = b[None, :] + delta
        FLUCT[g] = D["span"] * (lv - (D["span"] * lv).sum(0) / np.where(n > 0, n, 1.0))
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
    truth["station_trend_sd"] = station_trend_sd; truth["cluster_trend_sd"] = cluster_trend_sd; truth["pretrend"] = list(pretrend)
    return ys, truth


SIM_MASS = {"lob": [0.10, 0.08, 0.04, 0.0, 0.0], "shp": [0.08, 0.06, 0.02, 0.0, 0.0], "pur": [-0.05, -0.05, 0.0, 0.0, 0.0],
            "red": [0.05, 0.05, 0.04, 0.0, 0.0], "kelp": [0.0, 0.0, 0.0, 0.0, 0.0]}
MASS_TRUTH = dict(shape=1.5, sigma_station=0.4, mu_lambda=0.1, sigma_lambda=0.2, sigma_region_year=0.10, sigma_drift=0.04,
                  sigma_drift_region=0.04, rho=0.7, sigma_process=0.2, tau_D=0.2, zone=[0.0, 0.05, 0.10, 0.15])


MASS_SCENARIOS = ["std", "null", "trend", "bimodal", "visitnoise", "coupled"]


def simulate_mass(mass, truth, seed, scenario="std"):
    """Simulate group mean masses on the real measurement design at known values.
    std: from the fitted mass model. null: same with no protection effect on mass.
    trend: no protection effect, but a secular trend in mean mass shared by protected and reference stations
        (+2% a year, regions differing by SD 1% a year).
    Scenarios that break an assumption of the fitted model (protection effect as in std):
    bimodal: individual masses come from two cohorts (40% small, 60% about 4.5 times heavier) instead of a Gamma;
    visitnoise: every group's mean carries extra multiplicative noise (SD 0.25 on the log scale), as if the animals
        measured on a visit were not a random sample of the station;
    coupled: mean mass falls when that station's density is above its usual level (coefficient -0.3), as with a
        recruitment pulse; the fitted model treats density and mass as independent."""
    rng = np.random.default_rng(seed + 7919); tr = MASS_TRUTH; out = {}
    for g, M in mass.items():
        if M is None:
            out[g] = None; continue
        S, Mn = len(M["mass_stations"]), len(M["mpa_ba"]); T = M["span"].shape[0]; R = int(M["st_region"].max()) + 1
        inc = np.array(SIM_MASS[g], dtype=float) * (0.0 if scenario in ("null", "trend") else 1.0)
        lam = np.where(M["mpa_ba"] == 1, tr["mu_lambda"] + tr["sigma_lambda"] * rng.normal(size=Mn), rng.normal(0, 0.3, size=Mn))
        c = M["log_mean_mass"] + tr["sigma_station"] * rng.normal(size=S) + M["st_prot"] * lam[M["st_mpa"]]
        gam = tr["sigma_region_year"] * rng.normal(size=(R, T))
        al = np.cumsum(tr["sigma_drift"] * rng.normal(size=T)); br = np.cumsum(tr["sigma_drift_region"] * rng.normal(size=(R, T)), axis=1)
        gam = gam + (al - al.mean())[None, :] + br - br.mean(1, keepdims=True)
        if scenario == "trend":
            gam = gam + (0.02 + 0.01 * rng.normal(size=R))[:, None] * (np.arange(T) - (T - 1) / 2)[None, :]
        e = np.zeros((T, S)); e[0] = rng.normal(0, tr["sigma_process"] / np.sqrt(1 - tr["rho"] ** 2), size=S)
        for t in range(1, T):
            e[t] = tr["rho"] * e[t - 1] + tr["sigma_process"] * rng.normal(size=S)
        Dm = tr["tau_D"] * rng.normal(size=Mn)
        logM = c[None, :] + gam[M["st_region"]].T + e + M["grid_after"] * (M["grid_B"] @ inc + Dm[M["st_mpa"]][None, :] * M["grid_age"] / 10.0)
        if scenario == "coupled":
            logM = logM - 0.3 * FLUCT[g][:, M["dens_index"]]

        def draw(t_ix, s_ix, zone, n):
            mean = np.exp(logM[t_ix, s_ix] + np.array(tr["zone"])[zone])
            if scenario == "bimodal":
                sd = 0.3; mu = np.array([0.0, 1.5]); pr = np.array([0.4, 0.6]); norm = float(pr @ np.exp(mu + sd ** 2 / 2))
                y = np.empty(len(n))
                for i, k in enumerate(n.astype(int)):
                    comp = rng.random(k) < pr[1]
                    y[i] = np.exp(mu[comp.astype(int)] + sd * rng.normal(size=k)).mean() * mean[i] / norm
                return y
            conc = n * tr["shape"]
            y = rng.gamma(conc, mean / conc)
            if scenario == "visitnoise":
                y = y * np.exp(rng.normal(-0.25 ** 2 / 2, 0.25, size=len(y)))
            return np.maximum(y, 1e-9)

        M2 = dict(M); M2["g_mean"] = draw(M["g_t"], M["g_s"], M["g_zone"], M["g_n"])
        if "h_s" in M:
            M2["h_mean"] = draw(M["h_t"], M["h_s"], M["h_zone"], M["h_n"])
        out[g] = M2
        e10m = float(basis_row(10.0) @ inc)
        truth[g].update(E10_mass=e10m, E10_biomass=truth[g]["E10_total"] + e10m, mass_scenario=scenario)
    return out


def mass_checks(mc, data, edges, config, kw, n_draws, seed):
    """Predictive checks of the mass layer from one chain's draws (arrays; the parent process merges chains).
    Replicates every fitted group's mean mass from the posterior and compares summary statistics with the data;
    for held-out groups, counts how many replicates fall below the observed mean (their predictive rank)."""
    from numpyro.infer import Predictive
    mass = kw.get("mass") or {}
    tags = [g for g, M in mass.items() if M is not None]
    if not tags:
        return {}
    post = mc.get_samples(); n = len(next(iter(post.values()))); idx = np.linspace(0, n - 1, min(n_draws, n)).astype(int)
    thin = {k: v[idx] for k, v in post.items()}
    sites = [f"{g}:m_mu_groups" for g in tags] + [f"{g}:m_mu_holdout" for g in tags if "h_s" in mass[g]]
    pred = Predictive(model, posterior_samples=thin, return_sites=sites)(jax.random.PRNGKey(seed + 11), data, edges, config, **{**kw, "expose": True})
    rng = np.random.default_rng(seed + 13); out = {}
    for g in tags:
        M = mass[g]; k = np.asarray(thin[f"{g}:m_shape"]); mu = np.asarray(pred[f"{g}:m_mu_groups"])
        progs = sorted(set(M["g_prog"].tolist())); w = M["g_n"]; prot = M["g_prot"] > 0; post_ = M["g_post"] > 0; age = M["g_age"]

        def stats(y):
            ly = np.log(y); o = []
            for p_ in progs:
                m = M["g_prog"] == p_
                o += [ly[m].std(), np.quantile(ly[m], .05), np.quantile(ly[m], .95)]
            def lr(mask):
                a_, b_ = mask & prot, mask & ~prot
                if a_.sum() == 0 or b_.sum() == 0:
                    return np.nan
                return np.average(ly[a_], weights=w[a_]) - np.average(ly[b_], weights=w[b_])
            o += [lr(post_) - lr(~post_), lr(age >= 10) - lr((age >= 0) & (age <= 5))]
            return o
        names = [f"{s_}|prog{p_}" for p_ in progs for s_ in ("sd_log", "q05_log", "q95_log")] + ["contrast_after_minus_before", "contrast_age10plus_minus_age0to5"]
        conc = w[None, :] * k[:, M["g_prog"]]
        reps = np.array([stats(np.maximum(rng.gamma(conc[d], mu[d] / conc[d]), 1e-12)) for d in range(len(idx))])
        out[f"{g}|names"] = np.array(names); out[f"{g}|obs"] = np.array(stats(M["g_mean"])); out[f"{g}|reps"] = reps
        if "h_s" in M:
            mh = np.asarray(pred[f"{g}:m_mu_holdout"]); ch = M["h_n"][None, :] * k[:, M["h_prog"]]
            yh = rng.gamma(ch, mh / ch)
            out[f"{g}|h_below"] = (yh < M["h_mean"][None, :]).sum(0).astype(float); out[f"{g}|h_draws"] = np.array([float(len(idx))])
            out[f"{g}|h_prog"] = M["h_prog"].astype(float)
    return out


def merge_checks(parts, programs):
    """Combine per-chain predictive checks into the tables stored in the result file."""
    ppc, hold = {}, {}
    tags = sorted({k.split("|")[0] for k in parts[0]})
    for g in tags:
        names = [str(x) for x in parts[0][f"{g}|names"]]; obs = parts[0][f"{g}|obs"]; reps = np.concatenate([p[f"{g}|reps"] for p in parts], 0)
        ppc[g] = {}
        for i, nm in enumerate(names):
            if "|prog" in nm:
                st, pi = nm.split("|prog"); nm = f"{st}|{programs[g][int(pi)]}"
            r = reps[:, i]
            ppc[g][nm] = dict(observed=round(float(obs[i]), 4), rep_q05=round(float(np.nanquantile(r, .05)), 4), rep_median=round(float(np.nanmedian(r)), 4),
                              rep_q95=round(float(np.nanquantile(r, .95)), 4), p_rep_ge_obs=round(float(np.nanmean(r >= obs[i])), 3))
        if f"{g}|h_below" in parts[0]:
            pit = sum(p[f"{g}|h_below"] for p in parts) / sum(float(p[f"{g}|h_draws"][0]) for p in parts); hp = parts[0][f"{g}|h_prog"].astype(int)
            hold[g] = dict(n_groups=int(len(pit)), coverage90=round(float(((pit >= .05) & (pit <= .95)).mean()), 3), coverage50=round(float(((pit >= .25) & (pit <= .75)).mean()), 3),
                           mean_rank=round(float(pit.mean()), 3),
                           by_program={programs[g][int(q)]: dict(n=int((hp == q).sum()), coverage90=round(float(((pit[hp == q] >= .05) & (pit[hp == q] <= .95)).mean()), 3)) for q in sorted(set(hp.tolist()))})
    return ppc, hold


def fit(data, edges, ys, config, curve="rw", lag=1, beta_scale=0.5, paths=True, split=False, mass=None, warmup=600, samples=600, chains=4, seed=0,
        depth=10, accept=0.99, step0=0.01, thin=1, **common):
    data = {g: attach_counts(dict(D), ys[g]) for g, D in data.items()}
    kw = dict(curve=curve, lag=lag, beta_scale=beta_scale, paths=paths, split=split, mass=mass, **common)
    kern = NUTS(model, init_strategy=init_to_uniform(radius=0.5), target_accept_prob=accept, max_tree_depth=depth, step_size=step0)
    mc = MCMC(kern, num_warmup=warmup, num_samples=samples, num_chains=chains, thinning=thin,
              chain_method="parallel" if CPUS > 1 else "sequential", progress_bar=False)
    t0 = time.time()
    mc.run(jax.random.PRNGKey(seed), data, edges, config, **kw,
           extra_fields=("diverging", "num_steps", "potential_energy", "adapt_state.step_size"))
    return mc, time.time() - t0, data, kw


def keep_sites(mc):
    """Sites worth keeping from one MCMC run, as arrays shaped (chains, draws, ...), plus sampler diagnostics."""
    grouped = mc.get_samples(group_by_chain=True)
    ex = mc.get_extra_fields(group_by_chain=True)
    keep = {k: np.asarray(v) for k, v in grouped.items()
            if k in ("beta", "beta_prot") or k.split(":")[-1] in SITES or k.split(":")[-1].startswith("E10_via") or k.split(":")[-1] == "increments"}
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
        if k in ("beta", "beta_prot") or k.split(":")[-1].startswith("E10"):
            v2 = v.reshape(v.shape[0] * v.shape[1], -1)
            q[k] = dict(median=np.median(v2, 0).round(4).tolist(), q025=np.quantile(v2, .025, 0).round(4).tolist(),
                        q975=np.quantile(v2, .975, 0).round(4).tolist(), prob_negative=(v2 < 0).mean(0).round(4).tolist())
    out["intervals95"] = q
    out["edges"] = [f"{k}<-{j}" for k, j, _ in edges]
    pc = {k: v.reshape(v.shape[0], v.shape[1], -1).mean(1).round(3).tolist()
          for k, v in keep.items() if k in ("beta", "beta_prot") or k.split(":")[-1] in ("E10_total", "E10_direct", "sigma_process", "rho", "inv_sqrt_phi")}
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
    checks = [dict(np.load(tmp / f"checks{k}.npz", allow_pickle=False)) for k in range(a.chains) if (tmp / f"checks{k}.npz").exists()]
    for f in tmp.iterdir():
        f.unlink()
    tmp.rmdir()
    return keep, ex, secs, checks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["simulate", "fit"])
    ap.add_argument("--taxa", default="lob,shp,pur,red,kelp", help="tags, sources before targets: " + ", ".join(f"{k}={v}" for k, v in TAXA.items()))
    ap.add_argument("--edges", default=None, help='links as "target<source,...", default: ' + ",".join(f"{k}<{j}" for k, j in EDGES))
    ap.add_argument("--no-paths", action="store_true", help="fit the taxa side by side with no links")
    ap.add_argument("--split", action="store_true", help="separate link coefficients for station-level fluctuation and for the protection function")
    ap.add_argument("--mass", action="store_true", help="add the mean individual mass layer (biomass = density x mean mass)")
    ap.add_argument("--mass-variant", default="main", choices=["main", "b_lo", "b_hi", "ge25"],
                    help="b_lo/b_hi: allometric exponent -0.1/+0.1; ge25: urchin mass from animals >= 25 mm only (KFM, LTER)")
    ap.add_argument("--mass-drop", action="append", default=[], help="leave out a program's mass groups (KFM, LTER, PISCO, KFM-roving)")
    ap.add_argument("--size25", action="store_true", help="urchins: density on the >= 25 mm definition for KFM and LTER (use with --mass-variant ge25)")
    ap.add_argument("--omit-station", action="append", default=[], help='leave out a station, e.g. "KFM|a-k-01"')
    ap.add_argument("--mass-scenario", default="std", choices=MASS_SCENARIOS, help="simulate: how group mean masses are generated")
    ap.add_argument("--mass-ppc", type=int, default=0, help="posterior draws per chain for predictive checks of the mass layer")
    ap.add_argument("--mass-holdout", type=int, default=None, help="seed: hold out 20%% of interior station-years of mass data and predict them")
    ap.add_argument("--placebo-year", type=int, default=None, help="time placebo: KFM years before the real designation, fake designation year")
    ap.add_argument("--pretrend", action="store_true", help="sensitivity: differential linear trend at protected stations of before-and-after MPAs, extended through the after years")
    ap.add_argument("--evidence", default="all", choices=["all", "ba", "post"],
                    help="one taxon only: keep MPAs with before-and-after surveys (ba) or without (post); the split is fixed by the survey record")
    ap.add_argument("--permute-seed", type=int, default=None, help="randomization placebo: deal the protected labels out at random within each MPA and its references")
    ap.add_argument("--space-placebo-seed", type=int, default=None, help="space placebo: reference stations only, one relabelled as protected per MPA")
    ap.add_argument("--mass-scale", type=float, default=1.0, help="multiplier on the prior scale of the mass response increments")
    ap.add_argument("--mass-link", action="store_true", help="let mean mass depend on the same-year fluctuation of the taxon's density")
    ap.add_argument("--lag", type=int, default=1, choices=[0, 1, 2])
    ap.add_argument("--beta-scale", type=float, default=0.5, help="prior SD of the path coefficients")
    ap.add_argument("--curve", default="rw", choices=["piecewise", "rw"])
    ap.add_argument("--cluster-drift", action="store_true", help="local random-walk trend shared by each MPA's protected and reference stations")
    ap.add_argument("--station-trend", action="store_true", help="a slope of its own for every station")
    ap.add_argument("--sim-station-trend", type=float, default=0.0, help="simulate: SD of station slopes (log units per decade)")
    ap.add_argument("--sim-pretrend", type=float, nargs=2, default=[0.0, 0.0], metavar=("MEAN", "SD"), help="simulate: protected-minus-reference trend per decade at before-and-after MPAs, extended through the after years")
    ap.add_argument("--sim-cluster-trend", type=float, default=0.0, help="simulate: SD of yearly steps of a local random walk per MPA cluster")
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
    s25 = None
    if a.size25:
        s25 = find_beside(a.uses, "size25_share.csv")
        if s25 is None:
            raise SystemExit("size25_share.csv (from prep_mass.py) not found")
    if a.evidence != "all":
        if len(tags) != 1:
            raise SystemExit("--evidence ba/post is for one taxon at a time")
        D0 = build_multi(a.uses, tags, ed, lag=a.lag, regions=a.regions, omit_mpas=a.omit_mpa, omit_stations=a.omit_station, size25=s25)[0][tags[0]]
        drop = [m for m, f in zip(D0["mpas"], D0["mpa_ba"]) if (f == 1) != (a.evidence == "ba")]
        a.omit_mpa = list(a.omit_mpa) + drop
    data, edges = build_multi(a.uses, tags, ed, lag=a.lag, regions=a.regions, omit_mpas=a.omit_mpa, omit_stations=a.omit_station, size25=s25,
                              placebo_year=a.placebo_year, space_placebo_seed=a.space_placebo_seed, permute_seed=a.permute_seed)
    mass = {g: build_mass(a.uses, data[g], TAXA[g], variant=a.mass_variant, drop=a.mass_drop, holdout_seed=a.mass_holdout) for g in tags} if a.mass else None
    for D in data.values():
        D.pop("cells")
    config = {g: CONFIG[g] for g in tags}
    opts = dict(curve=a.curve, lag=a.lag, beta_scale=a.beta_scale, paths=not a.no_paths, split=a.split, warmup=a.warmup, samples=a.samples,
                chains=a.chains, seed=a.seed, depth=a.depth, accept=a.accept, step0=a.step0, thin=a.thin, drift=a.drift, inc_scale=a.inc_scale, cluster_drift=a.cluster_drift, station_trend=a.station_trend, pretrend=a.pretrend)
    if mass is not None:
        opts.update(mass_scale=a.mass_scale, mass_link=a.mass_link)
    meta = dict(mode=a.mode, evidence=a.evidence, permute_seed=a.permute_seed, placebo_year=a.placebo_year, space_placebo_seed=a.space_placebo_seed, relabelled={g: D["relabelled"] for g, D in data.items()},
                mass_holdout_seed=a.mass_holdout, omit_station=a.omit_station, size25=bool(a.size25), stations_removed={g: D["stations_removed"] for g, D in data.items()},
                units_size25_adjusted={g: D.get("units_size25_adjusted", 0) for g, D in data.items()}, chain_processes=(a.chains > 1 and not a.one_process), taxa={g: TAXA[g] for g in tags}, regions=a.regions, omit_mpa=a.omit_mpa,
                links={f"{k}<-{j}": {x: lk[x] for x in ("n_target_station_years", "n_linked", "n_stations_linked")} for k, j, lk in edges},
                n_units={g: int(len(D["y"])) for g, D in data.items()}, stations={g: len(D["stations"]) for g, D in data.items()}, **opts)
    if mass is not None:
        meta["mass_data"] = {g: (None if M is None else {k: M[k] for k in ("n_groups", "animals", "station_years_with_mass", "density_station_years", "kinds", "source", "variant", "dropped", "by_program")}
                                 | dict(stations=len(M["mass_stations"]), mpas_before_after=int(M["mpa_ba"].sum()), n_heldout=M["n_heldout"])) for g, M in mass.items()}
    if a.mode == "simulate":
        name = a.name or f"msim_{a.scenario}_{a.curve}_seed{a.seed}"
        ys, truth = simulate(data, edges, a.scenario, a.seed, lag=a.lag, direct=not a.no_direct, common_sd=a.common_sd,
                             station_trend_sd=a.sim_station_trend, cluster_trend_sd=a.sim_cluster_trend, pretrend=tuple(a.sim_pretrend))
        if mass is not None:
            mass = simulate_mass(mass, truth, a.seed, a.mass_scenario)
    else:
        name = a.name or f"mfit_{a.curve}{'_nopaths' if a.no_paths else ''}"
        ys = {g: D["y"] for g, D in data.items()}
    want_checks = mass is not None and (a.mass_ppc > 0 or a.mass_holdout is not None)
    nd = a.mass_ppc or 200
    if a.worker is not None:          # one chain, own seed; hand the kept draws back to the parent process
        opts.update(chains=1, seed=a.seed * 1000 + a.worker)
        mc, secs, data2, kw = fit(data, edges, ys, config, mass=mass, **opts)
        keep, ex = keep_sites(mc)
        np.savez(Path(a.worker_dir) / f"chain{a.worker}.npz", **keep, **{f"__{k}": v for k, v in ex.items()}, __seconds=np.array([secs]))
        if want_checks:
            np.savez(Path(a.worker_dir) / f"checks{a.worker}.npz", **mass_checks(mc, data2, edges, config, kw, nd, a.seed * 1000 + a.worker))
        return
    if a.chains > 1 and not a.one_process:
        keep, ex, secs, checks = run_chains(a, data, edges, ys, config, opts)
    else:
        mc, secs, data2, kw = fit(data, edges, ys, config, mass=mass, **opts)
        keep, ex = keep_sites(mc)
        checks = [mass_checks(mc, data2, edges, config, kw, nd, a.seed)] if want_checks else []
    res = summarize(keep, ex, edges)
    if checks and checks[0]:
        res["mass_ppc"], res["mass_holdout"] = merge_checks(checks, {g: D["programs"] for g, D in data.items()})
    if a.placebo_year is not None or a.space_placebo_seed is not None or a.permute_seed is not None:
        res["placebo"] = True
    if a.mode == "simulate":
        res["truth"] = truth; res["scenario"] = a.scenario
    else:
        res["mpas"] = {g: D["mpas"] for g, D in data.items()}
        if a.save_draws:
            try:
                np.savez_compressed(Path(a.out) / f"draws_{name}.npz", **{k.replace(":", "__"): v for k, v in keep.items()
                                    if k in ("beta", "beta_prot") or k.split(":")[-1].startswith("E10") or k.split(":")[-1] in ("g_knots", "increments")})
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
          "min_neff", res["min_neff"], "beta", res.get("beta", {}).get("mean"), "beta_prot", res.get("beta_prot", {}).get("mean"))


if __name__ == "__main__":
    main()
