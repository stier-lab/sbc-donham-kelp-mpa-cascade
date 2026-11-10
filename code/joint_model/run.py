"""Fit, simulate and summarize the joint density model.

Examples (from the repository root):
  python3 code/joint_model/run.py simulate --taxon "S. purpuratus" --scenario effect --seed 1
  python3 code/joint_model/run.py fit --taxon "S. purpuratus" --curve piecewise
"""
import argparse, json, time
from pathlib import Path
import numpy as np
import numpyro
numpyro.enable_x64()
import jax
from numpyro import handlers
from numpyro.infer import MCMC, NUTS, init_to_value
import sys
sys.path.insert(0, str(Path(__file__).parent))
from data import build, age_basis
from model import model

USES = "outputs/integrated_review_20261103/private/supported_measurement_uses.parquet"
KEYS = ["E10", "increments", "tau_L", "L_mpa", "g_knots", "tau_D", "sigma_shock", "rho", "sigma_process", "sigma_region_year", "sigma_visit",
        "sigma_b", "inv_sqrt_phi", "mu_lambda", "sigma_lambda", "sigma_curve", "zone", "D_mpa", "shock", "lambda"]

TRUTH = dict(  # scales chosen to resemble the data profile; effects chosen by scenario
    sigma_b=1.0, sigma_region_year=0.6, rho=0.90, sigma_process=0.65, sigma_visit=0.4, inv_sqrt_phi=1.2,
    mu_lambda=0.5, sigma_lambda=0.7, tau_D=1.5, sigma_shock=0.15, sigma_curve=1.0)
SCENARIOS = {"null": [0, 0, 0, 0, 0], "effect": [-0.15, -0.35, -0.30, -0.10, 0.0],  # E10 = -0.65 (x0.52 per decade)
             "small": [-0.10, -0.10, -0.10, -0.10, -0.10],
             "large": [1.2, 0.6, 0.4, 0.0, 0.0]}                                       # E10 = +2.0 (x7.4 per decade)                           # E10 = -0.25 (x0.78 per decade)


def simulate(D, B, scenario, seed):
    """Independent NumPy simulator of the model on the real survey design, at known parameter values."""
    rng = np.random.default_rng(seed)
    S, M, P, R, T = len(D["stations"]), len(D["mpas"]), len(D["programs"]), len(D["regions"]), len(D["years"])
    tr = TRUTH; inc = np.array(SCENARIOS[scenario], dtype=float)
    lam = np.where(D["mpa_ba"] == 1, tr["mu_lambda"] + tr["sigma_lambda"] * rng.normal(size=M), rng.normal(0, 1.0, size=M))
    b = D["log_mean_density"] + tr["sigma_b"] * rng.normal(size=S) + D["st_prot"] * lam[D["st_mpa"]]
    gam = tr["sigma_region_year"] * rng.normal(size=(R, T))
    delta = np.zeros((T, S)); delta[0] = rng.normal(0, tr["sigma_process"] / np.sqrt(1 - tr["rho"] ** 2), size=S)
    for t in range(1, T):
        delta[t] = tr["rho"] * delta[t - 1] + tr["sigma_process"] * rng.normal(size=S)
    Dm = tr["tau_D"] * rng.normal(size=M); omega = tr["sigma_shock"] * rng.normal(size=T)
    nu = tr["sigma_visit"] * rng.normal(size=len(D["visit_prog"]))
    zeta = np.array([0.0, -0.3, -1.0, -3.0]); mon = np.zeros(12)
    s, t = D["s"], D["t"]; m = D["st_mpa"][s]
    F = D["after"] * (B @ inc + Dm[m] * D["age"] / 10.0 + omega[t])
    eta = b[s] + gam[D["st_region"][s], t] + delta[t, s] + F + nu[D["v"]] + zeta[D["zone"]] + mon[D["month"]]
    mu = np.exp(eta + D["log_area"])[D["unit_cell"]]; phi = 1.0 / tr["inv_sqrt_phi"] ** 2
    y = rng.negative_binomial(phi, phi / (phi + mu))
    truth = dict(tr, increments=inc.tolist(), D_mpa=Dm.tolist(), shock=omega.tolist(), **{"lambda": lam.tolist()})
    return y, truth


def attach_counts(D, y):
    """Sufficient statistics of unit counts y: cell totals, and distinct (count, program) pairs with multiplicity."""
    y = np.asarray(y)
    uq, mult = np.unique(np.stack([y, D["unit_prog"]], 1), axis=0, return_counts=True)
    D.update(cell_y=np.bincount(D["unit_cell"], weights=y, minlength=len(D["log_area"])),
             uniq_y=uq[:, 0].astype(float), uniq_prog=uq[:, 1], uniq_mult=mult.astype(float),
             prog_units=np.bincount(D["unit_prog"], minlength=len(D["programs"])).astype(float))
    return D


def fit(D, B, y, curve="piecewise", shock=True, warmup=600, samples=600, chains=2, seed=0, depth=9, centered=True, level_dev=False):
    kern = NUTS(model, init_strategy=init_to_value(values={"delta": np.zeros((len(D["years"]), len(D["stations"])))}), target_accept_prob=0.95, max_tree_depth=depth)
    mc = MCMC(kern, num_warmup=warmup, num_samples=samples, num_chains=chains, chain_method="sequential", progress_bar=False)
    D = attach_counts(dict(D), y)
    t0 = time.time(); mc.run(jax.random.PRNGKey(seed), D, B, curve=curve, shock=shock, y=y, centered=centered, level_dev=level_dev, extra_fields=("diverging", "num_steps"))
    steps = float(np.asarray(mc.get_extra_fields()["num_steps"]).mean())   # forces computation to finish before timing
    print("mean leapfrog steps", round(steps))
    return mc, time.time() - t0


def summarize(mc, keys=KEYS):
    from numpyro.diagnostics import summary
    s = summary(mc.get_samples(group_by_chain=True), prob=0.9, group_by_chain=True)
    out = {}
    for k in keys:
        if k not in s:
            continue
        v = s[k]
        out[k] = {f: np.asarray(v[f]).round(4).tolist() for f in ("mean", "5.0%", "95.0%", "n_eff", "r_hat")}
    out["divergences"] = int(np.asarray(mc.get_extra_fields()["diverging"]).sum())
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["simulate", "fit"])
    ap.add_argument("--taxon", default="S. purpuratus")
    ap.add_argument("--curve", default="piecewise", choices=["piecewise", "rw", "linear"])
    ap.add_argument("--no-shock", action="store_true")
    ap.add_argument("--level-dev", action="store_true", help="add an MPA-specific level shift after designation")
    ap.add_argument("--noncentered", action="store_true", help="innovations parameterization of the station process")
    ap.add_argument("--scenario", default="effect", choices=list(SCENARIOS))
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--warmup", type=int, default=600)
    ap.add_argument("--samples", type=int, default=600)
    ap.add_argument("--chains", type=int, default=2)
    ap.add_argument("--placebo-year", type=int, default=None)
    ap.add_argument("--uses", default=USES)
    ap.add_argument("--out", default="outputs/joint_model_20261109")
    a = ap.parse_args()
    D = build(a.uses, a.taxon, placebo_year=a.placebo_year); B = age_basis(D["age"])
    cells = D.pop("cells")
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    tag = a.taxon.replace(". ", "_").replace(" ", "_")
    if a.mode == "simulate":
        y, truth = simulate(D, B, a.scenario, a.seed)
        mc, secs = fit(D, B, y, a.curve, not a.no_shock, a.warmup, a.samples, a.chains, a.seed)
        res = summarize(mc)
        res["truth"] = {k: (np.asarray(v).round(4).tolist()) for k, v in truth.items()}
        res["truth"]["E10"] = float(np.dot(SCENARIOS[a.scenario], [1, 1, .5, 0, 0]))
        res["truth"]["g_knots"] = np.cumsum(SCENARIOS[a.scenario]).tolist()
        name = f"sim_{tag}_{a.scenario}_{a.curve}_seed{a.seed}.json"
    else:
        mc, secs = fit(D, B, D["y"], a.curve, not a.no_shock, a.warmup, a.samples, a.chains, a.seed, centered=not a.noncentered, level_dev=a.level_dev)
        res = summarize(mc)
        res["mpas"] = D["mpas"]; res["mpa_before_after"] = D["mpa_ba"].tolist(); res["years"] = D["years"].tolist()
        post = mc.get_samples()
        e10 = np.asarray(post["E10"])
        res["E10_prob_negative"] = float((e10 < 0).mean()); res["E10_prob_beyond_0.22"] = float((np.abs(e10) > 0.223).mean())
        res["E10_prob_within_0.22"] = float((np.abs(e10) < 0.223).mean())
        np.savez_compressed(out / f"draws_{'placebo_' if a.placebo_year else ''}{tag}_{a.curve}{'_noshock' if a.no_shock else ''}{'_leveldev' if a.level_dev else ''}.npz",
                            **{k: np.asarray(v) for k, v in post.items() if k in KEYS})
        name = f"{'placebo' if a.placebo_year else 'fit'}_{tag}_{a.curve}{'_noshock' if a.no_shock else ''}{'_leveldev' if a.level_dev else ''}.json"
    res.update(seconds=round(secs), n_obs=int(len(D["y"])), stations=len(D["stations"]), n_mpas=len(D["mpas"]),
               station_process="innovations" if a.noncentered else "centered", warmup=a.warmup, samples=a.samples, chains=a.chains, taxon=a.taxon, curve=a.curve, shock=not a.no_shock, level_dev=a.level_dev)
    (out / name).write_text(json.dumps(res, indent=1))
    print(name, "seconds", round(secs), "div", res["divergences"], "E10", res.get("E10"))


if __name__ == "__main__":
    main()
