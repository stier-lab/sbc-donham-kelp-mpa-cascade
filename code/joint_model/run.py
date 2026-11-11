"""Fit, simulate and check the joint density model.

Examples (from the repository root):
  python3 code/joint_model/run.py fit --taxon "S. purpuratus" --curve rw
  python3 code/joint_model/run.py simulate --scenario effect --seed 1 --curve rw
  python3 code/joint_model/run.py fit --taxon "S. pulcher" --curve rw --noncentered --level-dev --ppc 200
"""
import argparse, json, os, time
from pathlib import Path
import numpy as np
import numpyro


def write_retry(path, text, tries=60, wait=20):
    """Write a result file, retrying when the file system refuses (GRIT scratch was 100% full on 2026-11-10 and
    writes failed intermittently). Each attempt uses a NEW temporary file: on BeeGFS a file is pinned to a storage
    target when it is created, so rewriting the same file retries the same full target."""
    path = Path(path); n = len(text.encode())
    for i in range(tries):
        tmp = path.with_name(f"{path.name}.{os.getpid()}.{i}.part")
        try:
            tmp.write_text(text)
            if tmp.stat().st_size == n:
                os.replace(tmp, path)
                return
            raise OSError("short write")
        except OSError as e:
            print(f"write failed ({e}); retry {i + 1}/{tries}", flush=True)
            try:
                tmp.unlink()
            except OSError:
                pass
        fb = os.environ.get("MPA_JOINT_FALLBACK")
        if fb and i >= 4:   # scratch keeps refusing: small JSON results may go to the fallback directory (size-capped)
            try:
                fbd = Path(fb); fbd.mkdir(parents=True, exist_ok=True)
                used = sum(f.stat().st_size for f in fbd.iterdir() if f.is_file())
                if used + n < 28e6:
                    (fbd / path.name).write_text(text)
                    if (fbd / path.name).stat().st_size == n:
                        print(f"wrote fallback {fbd / path.name}", flush=True)
                        return
            except OSError as e:
                print(f"fallback write failed ({e})", flush=True)
        time.sleep(wait)
    raise OSError(f"could not write {path}")


def _early_cpus():
    """Parallel chains need the host device count set before JAX starts."""
    import sys
    if "--cpus" in sys.argv:
        n = int(sys.argv[sys.argv.index("--cpus") + 1])
        if n > 1:
            numpyro.set_host_device_count(n)
            return n
    return 1


CPUS = _early_cpus()
numpyro.enable_x64()
import jax
from numpyro.infer import MCMC, NUTS, Predictive, init_to_uniform, init_to_value
import sys
sys.path.insert(0, str(Path(__file__).parent))
from data import build, age_basis
from model import model

USES = "outputs/integrated_review_20261103/private/supported_measurement_uses.parquet"
KEYS = ["sigma_drift", "sigma_drift_region", "drift", "zi_prob", "E10", "increments", "tau_L", "L_mpa", "g_knots", "tau_D", "sigma_shock", "rho", "sigma_process", "sigma_region_year",
        "sigma_visit", "sigma_b", "inv_sqrt_phi", "mu_lambda", "sigma_lambda", "sigma_curve", "zone", "D_mpa", "shock", "lambda"]
HYPERS = ["mu_b", "sigma_b", "mu_lambda", "sigma_lambda", "sigma_region_year", "rho_raw", "sigma_process", "sigma_curve",
          "inc_first", "tau_D", "tau_L", "sigma_shock", "sigma_visit", "inv_sqrt_phi"]

# Simulation truths. "urchin": variance values from the first purple-urchin fit. "lobster": from the lobster fit
# (zero-heavy counts, step-like response with MPA-specific level shifts). "level" shifts the mean log density from the
# log of the observed mean so that the simulated share of zero counts matches the real data (0.29 urchin, 0.79 lobster).
TRUTHS = {
    "urchin": dict(sigma_b=1.0, sigma_region_year=0.6, rho=0.90, sigma_process=0.65, sigma_visit=0.4, inv_sqrt_phi=1.2,
                   mu_lambda=0.5, sigma_lambda=0.7, tau_D=1.5, tau_L=0.0, sigma_shock=0.15, zeta=[0.0, -0.3, -1.0, -3.0], level=-2.0),
    "lobster": dict(sigma_b=1.2, sigma_region_year=0.4, rho=0.86, sigma_process=0.35, sigma_visit=0.45, inv_sqrt_phi=1.25,
                    mu_lambda=-0.5, sigma_lambda=1.0, tau_D=0.6, tau_L=0.8, sigma_shock=0.4, zeta=[0.0, -0.9, -1.3, -1.8], level=-3.4),
}
SCENARIOS = {"null": [0, 0, 0, 0, 0],
             "effect": [-0.15, -0.35, -0.30, -0.10, 0.0],     # E10 = -0.65 (x0.52 per decade)
             "small": [-0.10, -0.10, -0.10, -0.10, -0.10],    # E10 = -0.25 (x0.78 per decade)
             "large": [1.2, 0.6, 0.4, 0.0, 0.0],              # E10 = +2.00 (x7.4 per decade)
             "lobster_effect": [1.0, 0.7, 0.45, 0.4, 0.4]}    # E10 = +1.925 (x6.9 per decade)
E10_WEIGHTS = np.array([1, 1, 0.5, 0, 0])


def simulate(D, B, scenario, seed, truth="urchin", trend=0.0, trend_sd=0.0):
    """Independent NumPy simulator of the model on the real survey design, at known parameter values."""
    rng = np.random.default_rng(seed)
    S, M, R, T = len(D["stations"]), len(D["mpas"]), len(D["regions"]), len(D["years"])
    tr = TRUTHS[truth]; inc = np.array(SCENARIOS[scenario], dtype=float)
    lam = np.where(D["mpa_ba"] == 1, tr["mu_lambda"] + tr["sigma_lambda"] * rng.normal(size=M), rng.normal(0, 1.0, size=M))
    b = D["log_mean_density"] + tr["level"] + tr["sigma_b"] * rng.normal(size=S) + D["st_prot"] * lam[D["st_mpa"]]
    gam = tr["sigma_region_year"] * rng.normal(size=(R, T))
    # optional secular trend shared by every station (log units per year). The fitted model has no trend term:
    # region-year effects are independent across years. Used to test whether a regional trend leaks into E10.
    slopes = trend + trend_sd * rng.normal(size=R)      # trend_sd > 0: each region has its own trend
    gam = gam + slopes[:, None] * (np.asarray(D["years"], dtype=float) - float(np.mean(D["years"])))[None, :]
    delta = np.zeros((T, S)); delta[0] = rng.normal(0, tr["sigma_process"] / np.sqrt(1 - tr["rho"] ** 2), size=S)
    for t in range(1, T):
        delta[t] = tr["rho"] * delta[t - 1] + tr["sigma_process"] * rng.normal(size=S)
    Dm = tr["tau_D"] * rng.normal(size=M); omega = tr["sigma_shock"] * rng.normal(size=T)
    Lm = tr["tau_L"] * rng.normal(size=M)
    nu = tr["sigma_visit"] * rng.normal(size=len(D["visit_prog"]))
    zeta = np.array(tr["zeta"])
    s, t = D["s"], D["t"]; m = D["st_mpa"][s]
    F = D["after"] * (B @ inc + Lm[m] + Dm[m] * D["age"] / 10.0 + omega[t])
    eta = b[s] + gam[D["st_region"][s], t] + delta[t, s] + F + nu[D["v"]] + zeta[D["zone"]]
    mu = np.exp(eta + D["log_area"])[D["unit_cell"]]; phi = 1.0 / tr["inv_sqrt_phi"] ** 2
    y = rng.negative_binomial(phi, phi / (phi + mu))
    out = {k: v for k, v in tr.items() if k not in ("zeta", "level")}
    out.update(trend=trend, trend_sd=trend_sd, E10=float(inc @ E10_WEIGHTS), g_knots=np.cumsum(inc).tolist(), zero_fraction=float((y == 0).mean()))
    return y, out


def attach_counts(D, y):
    """Sufficient statistics of unit counts y: cell totals, and distinct (count, program) pairs with multiplicity."""
    y = np.asarray(y)
    uq, mult = np.unique(np.stack([y, D["unit_prog"]], 1), axis=0, return_counts=True)
    D.update(cell_y=np.bincount(D["unit_cell"], weights=y, minlength=len(D["log_area"])),
             cell_nzero=np.bincount(D["unit_cell"], weights=(y == 0).astype(float), minlength=len(D["log_area"])),
             cell_npos=np.bincount(D["unit_cell"], weights=(y > 0).astype(float), minlength=len(D["log_area"])),
             uniq_y=uq[:, 0].astype(float), uniq_prog=uq[:, 1], uniq_mult=mult.astype(float),
             prog_units=np.bincount(D["unit_prog"], minlength=len(D["programs"])).astype(float))
    return D


def fit(D, B, y, curve="rw", shock=True, warmup=600, samples=600, chains=4, seed=0, depth=10, centered=True,
        level_dev=False, tau_d_scale=1.0, center_b=False, init="tight", dense=False, accept=0.95, marg_b=False, step0=1.0,
        drift="none", zi=False, inc_scale=1.0):
    """init "tight": every unconstrained parameter starts uniformly within +-0.5 (scales start between 0.6 and 1.65).
    init "legacy": NumPyro's default +-2 start (scales between 0.14 and 7.4) with station states at zero; it can
    start a chain in a degenerate everything-is-count-noise mode from which it does not escape."""
    D = attach_counts(dict(D), y)
    kw = dict(curve=curve, shock=shock, centered=centered, level_dev=level_dev, tau_d_scale=tau_d_scale, center_b=center_b, marg_b=marg_b,
              drift=drift, zi=zi, inc_scale=inc_scale)
    strat = init_to_uniform(radius=0.5) if init == "tight" else \
        init_to_value(values={"delta": np.zeros((len(D["years"]), len(D["stations"])))})
    dm = False
    if dense:  # dense mass matrix for the hyperparameter block only
        present = {"piecewise": {"inc_first"}, "rw": set(), "linear": {"inc_first", "sigma_curve"}}[curve]
        drop = set(present) | (set() if level_dev else {"tau_L"}) | (set() if shock else {"sigma_shock"})
        dm = [tuple(h for h in HYPERS if h not in drop)]
    kern = NUTS(model, init_strategy=strat, target_accept_prob=accept, max_tree_depth=depth, dense_mass=dm, step_size=step0)
    mc = MCMC(kern, num_warmup=warmup, num_samples=samples, num_chains=chains,
              chain_method="parallel" if CPUS > 1 else "sequential", progress_bar=False)
    t0 = time.time()
    mc.run(jax.random.PRNGKey(seed), D, B, y=y, extra_fields=("diverging", "num_steps", "potential_energy", "adapt_state.step_size"), **kw)
    ex = mc.get_extra_fields(group_by_chain=True)
    steps = np.asarray(ex["num_steps"])
    return mc, time.time() - t0, D, kw, ex


def summarize(mc, ex, keys=KEYS):
    from numpyro.diagnostics import summary
    grouped = mc.get_samples(group_by_chain=True)
    s = summary({k: v for k, v in grouped.items() if k in keys}, prob=0.9, group_by_chain=True)
    out = {}; worst, least = 1.0, np.inf
    for k in keys:
        if k not in s:
            continue
        v = s[k]
        out[k] = {f: np.asarray(v[f]).round(4).tolist() for f in ("mean", "5.0%", "95.0%", "n_eff", "r_hat")}
        worst = max(worst, float(np.nanmax(v["r_hat"]))); least = min(least, float(np.nanmin(v["n_eff"])))
    out["worst_rhat"] = round(worst, 4); out["min_neff"] = round(least, 1)
    out["divergences"] = int(np.asarray(ex["diverging"]).sum())
    pc = {k: np.asarray(grouped[k]).reshape(grouped[k].shape[0], grouped[k].shape[1], -1).mean(1).round(3).tolist()
          for k in ("E10", "inv_sqrt_phi", "sigma_b", "sigma_region_year", "sigma_process", "sigma_visit", "rho", "tau_D") if k in grouped}
    pc["potential_energy"] = np.asarray(ex["potential_energy"]).mean(1).round(1).tolist()
    pc["step_size"] = np.asarray(ex["adapt_state.step_size"])[:, -1].round(5).tolist()
    pc["mean_leapfrog_steps"] = np.asarray(ex["num_steps"]).mean(1).round(0).tolist()
    pc["divergences"] = np.asarray(ex["diverging"]).sum(1).tolist()
    out["per_chain"] = pc
    return out


def posterior_predictive(mc, D, B, kw, n_draws, seed):
    """Replicate unit counts from posterior draws, with new visit noise and new count noise, holding each
    station-year's latent density at its posterior draw. Tests the count model (zeros, dispersion); the
    protected/reference contrasts are conditional on the latent states, so they are a weak check."""
    post = mc.get_samples(); n = len(post["E10"]); idx = np.linspace(0, n - 1, min(n_draws, n)).astype(int)
    thin = {k: v[idx] for k, v in post.items()}
    eta = np.asarray(Predictive(model, posterior_samples=thin, return_sites=["eta_novisit"])(
        jax.random.PRNGKey(seed + 1), D, B, y=1, expose=True, **kw)["eta_novisit"])
    rng = np.random.default_rng(seed); P = len(D["programs"])
    sv = np.asarray(thin["sigma_visit"]); phi = 1.0 / np.asarray(thin["inv_sqrt_phi"]) ** 2
    uc, up = D["unit_cell"], D["unit_prog"]; y = np.asarray(D["y"])
    cell_T = D["mpa_T"][D["st_mpa"][D["s"]]]; cell_year = D["years"][D["t"]]
    u_prot = D["prot"][uc] > 0; u_age = (cell_year - cell_T)[uc]; u_area = D["unit_area"]

    def stats(yy):
        o = {}
        for p in range(P):
            m = up == p; nm = D["programs"][p]
            o[f"zero_fraction|{nm}"] = float((yy[m] == 0).mean())
            o[f"sd_log1p|{nm}"] = float(np.log1p(yy[m]).std())
            o[f"q99|{nm}"] = float(np.quantile(yy[m], 0.99))
        def lr(mask):
            a, b = mask & u_prot, mask & ~u_prot
            if a.sum() == 0 or b.sum() == 0:
                return np.nan
            return float(np.log((yy[a].sum() + 0.5) / u_area[a].sum()) - np.log((yy[b].sum() + 0.5) / u_area[b].sum()))
        o["logratio_after_minus_before"] = lr(u_age >= 0) - lr(u_age < 0)
        o["logratio_age10plus_minus_age0to5"] = lr(u_age >= 10) - lr((u_age >= 0) & (u_age <= 5))
        return o

    obs = stats(y); reps = []
    for d in range(len(idx)):
        nu = sv[d][D["visit_prog"]] * rng.normal(size=len(D["visit_prog"]))
        mu = np.exp(eta[d] + nu[D["v"]] + D["log_area"])[uc]; ph = phi[d][up]
        yr = rng.negative_binomial(ph, ph / (ph + mu))
        if "zi_prob" in thin:
            yr = np.where(rng.random(len(yr)) < np.asarray(thin["zi_prob"])[d][up], 0, yr)
        reps.append(stats(yr))
    out = {}
    for k, v in obs.items():
        r = np.array([x[k] for x in reps])
        out[k] = dict(observed=round(v, 4), rep_q05=round(float(np.nanquantile(r, .05)), 4), rep_median=round(float(np.nanmedian(r)), 4),
                      rep_q95=round(float(np.nanquantile(r, .95)), 4), p_rep_ge_obs=round(float(np.nanmean(r >= v)), 3))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["simulate", "fit"])
    ap.add_argument("--taxon", default="S. purpuratus")
    ap.add_argument("--curve", default="rw", choices=["piecewise", "rw", "linear"],
                    help="rw: random-walk increments (main). piecewise: independent increments (conservative bound).")
    ap.add_argument("--no-shock", action="store_true")
    ap.add_argument("--level-dev", action="store_true", help="add an MPA-specific level shift after designation")
    ap.add_argument("--noncentered", action="store_true", help="innovations parameterization of the station process")
    ap.add_argument("--center-b", action="store_true", help="sample station baselines directly")
    ap.add_argument("--marg-b", action="store_true", help="integrate station baselines out analytically (centered process only)")
    ap.add_argument("--step0", type=float, default=1.0, help="initial step size for warm-up adaptation")
    ap.add_argument("--accept", type=float, default=0.95, help="target acceptance probability")
    ap.add_argument("--dense", action="store_true", help="dense mass matrix for the hyperparameter block")
    ap.add_argument("--init", default="tight", choices=["tight", "legacy"])
    ap.add_argument("--depth", type=int, default=10)
    ap.add_argument("--tau-d-scale", type=float, default=1.0)
    ap.add_argument("--scenario", default="effect", choices=list(SCENARIOS))
    ap.add_argument("--truth", default="urchin", choices=list(TRUTHS))
    ap.add_argument("--trend-sd", type=float, default=0.0, help="simulate: SD of region-specific trends around --trend")
    ap.add_argument("--drift", default="none", choices=["none", "global", "regional"], help="experimental: random-walk year effect")
    ap.add_argument("--zi", action="store_true", help="experimental: zero-inflated counts")
    ap.add_argument("--inc-scale", type=float, default=1.0, help="experimental: multiplier on the response-increment prior scale")
    ap.add_argument("--regions", default="seven", choices=["seven", "three", "one"])
    ap.add_argument("--trend", type=float, default=0.0, help="simulate: secular trend in log density per year, all stations")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--warmup", type=int, default=600)
    ap.add_argument("--samples", type=int, default=600)
    ap.add_argument("--chains", type=int, default=4)
    ap.add_argument("--cpus", type=int, default=1, help="run chains in parallel on this many cores")
    ap.add_argument("--placebo-year", type=int, default=None)
    ap.add_argument("--space-placebo-seed", type=int, default=None)
    ap.add_argument("--omit-mpa", action="append", default=[])
    ap.add_argument("--omit-station", action="append", default=[])
    ap.add_argument("--pisco-cap", type=float, default=None)
    ap.add_argument("--ppc", type=int, default=0, help="number of posterior draws for predictive checks")
    ap.add_argument("--save-draws", action="store_true", help="also save posterior draws (.npz) next to the result")
    ap.add_argument("--name", default=None, help="result file name (without .json)")
    ap.add_argument("--uses", default=USES)
    ap.add_argument("--out", default="outputs/joint_model_20261109")
    a = ap.parse_args()
    D = build(a.uses, a.taxon, placebo_year=a.placebo_year, omit_mpas=a.omit_mpa, omit_stations=a.omit_station,
              space_placebo_seed=a.space_placebo_seed, pisco_cap=a.pisco_cap, regions=a.regions)
    B = age_basis(D["age"]); D.pop("cells")
    out = Path(a.out)
    try:
        out.mkdir(parents=True, exist_ok=True)
    except OSError:
        pass
    tag = a.taxon.replace(". ", "_").replace(" ", "_")
    opts = dict(curve=a.curve, shock=not a.no_shock, warmup=a.warmup, samples=a.samples, chains=a.chains, seed=a.seed,
                depth=a.depth, centered=not a.noncentered, level_dev=a.level_dev, tau_d_scale=a.tau_d_scale,
                center_b=a.center_b, init=a.init, dense=a.dense, marg_b=a.marg_b, step0=a.step0, accept=a.accept,
                drift=a.drift, zi=a.zi, inc_scale=a.inc_scale)
    meta = dict(taxon=a.taxon, mode=a.mode, n_obs=int(len(D["y"])), cells=int(len(D["log_area"])), stations=len(D["stations"]),
                n_mpas=len(D["mpas"]), programs=D["programs"], stations_removed=D["stations_removed"],
                relabelled=D["relabelled"], units_capped=D["units_capped"], placebo_year=a.placebo_year,
                omit_mpa=a.omit_mpa, omit_station=a.omit_station, pisco_cap=a.pisco_cap, regions=a.regions, n_regions=len(D["regions"]), **opts)
    if a.mode == "simulate":
        name = a.name or f"sim_{tag}_{a.scenario}_{a.curve}_seed{a.seed}"
        y, truth = simulate(D, B, a.scenario, a.seed, a.truth, a.trend, a.trend_sd)
        mc, secs, D2, kw, ex = fit(D, B, y, **opts)
        full = summarize(mc, ex)
        res = {k: full[k] for k in ("E10", "sigma_drift", "sigma_drift_region", "g_knots", "tau_D", "tau_L", "sigma_shock", "rho", "sigma_process", "sigma_region_year",
                                    "sigma_visit", "inv_sqrt_phi", "mu_lambda", "divergences", "worst_rhat", "min_neff") if k in full}
        res["truth"] = truth; res["scenario"] = a.scenario; res["truth_set"] = a.truth
    else:
        if (a.omit_mpa or a.omit_station) and D["stations_removed"] == 0:
            name = a.name or f"fit_{tag}_{a.curve}"
            write_retry(out / f"{name}.json", json.dumps(dict(ok=True, not_applicable=True, reason="station or MPA not in this taxon's data", **meta), indent=1))
            print(name, "not applicable"); return
        kind = "placebo" if (a.placebo_year or a.space_placebo_seed is not None) else "fit"
        name = a.name or f"{kind}_{tag}_{a.curve}{'_noshock' if a.no_shock else ''}{'_leveldev' if a.level_dev else ''}"
        mc, secs, D2, kw, ex = fit(D, B, D["y"], **opts)
        res = summarize(mc, ex)
        res["mpas"] = D["mpas"]; res["mpa_before_after"] = D["mpa_ba"].tolist(); res["years"] = D["years"].tolist()
        post = mc.get_samples(); e10 = np.asarray(post["E10"])
        res["E10_prob_negative"] = float((e10 < 0).mean()); res["E10_prob_beyond_0.22"] = float((np.abs(e10) > 0.223).mean())
        res["E10_prob_within_0.22"] = float((np.abs(e10) < 0.223).mean())
        if a.ppc:
            res["ppc"] = posterior_predictive(mc, D2, B, kw, a.ppc, a.seed)
        if a.save_draws:
            try:
                np.savez_compressed(out / f"draws_{name}.npz", **{k: np.asarray(v) for k, v in post.items() if k in KEYS})
            except OSError as e:
                print(f"draws not saved ({e})", flush=True)
    res.update(meta, seconds=round(secs), station_process="centered" if not a.noncentered else "innovations")
    write_retry(out / f"{name}.json", json.dumps(res, indent=1))
    print(name, "seconds", round(secs), "div", res["divergences"], "worst_rhat", res.get("worst_rhat"), "min_neff", res.get("min_neff"), "E10", res.get("E10"))


if __name__ == "__main__":
    main()
