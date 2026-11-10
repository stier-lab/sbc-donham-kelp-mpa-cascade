"""Profile survey-unit counts, sizes and station histories to set priors for the joint hierarchical model.

Reads the integrated review (outputs/integrated_review_20261103) and writes aggregate tables only
(no unit- or record-level data) to outputs/joint_model_design_20261109/. Review-only; changes no inputs.
Run from the repository root: python3 code/tools/profile_joint_model_data.py
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
IR = ROOT / "outputs/integrated_review_20261103"
OUT = ROOT / "outputs/joint_model_design_20261109"
OUT.mkdir(parents=True, exist_ok=True)
FOCAL = ["M. franciscanus", "S. purpuratus", "P. interruptus", "S. pulcher", "M. pyrifera"]


def split_series(d):
    d[["MPA", "Taxa", "Source", "Resp"]] = d.series_id.str.split("|", expand=True)
    return d


def unit_profile():
    u = pd.read_parquet(IR / "private/measurement_units.parquet")
    u = u[u.Source != "Landsat"].copy()
    u["count"] = u.density * u.area
    u["integer"] = np.isclose(u["count"], np.round(u["count"]), atol=1e-6)
    u["count"] = np.round(u["count"])
    rows = []
    for (s, t), d in u.groupby(["Source", "Taxa"]):
        v = d.groupby(["site", "date"])["count"].agg(["size", "mean", "var"]).reset_index()
        vv = v[(v["size"] >= 3) & (v["mean"] > 0)]
        x, y = vv["mean"] ** 2, vv["var"] - vv["mean"]
        inv_phi = float((x * y).sum() / (x * x).sum()) if len(vv) else np.nan
        vy = d.groupby(["site", "year"]).date.nunique()
        rows.append(dict(Source=s, Taxa=t, units=len(d), stations=d.site.nunique(),
                         first_year=int(d.year.min()), last_year=int(d.year.max()),
                         unit_area_m2=f"{d.area.min():g}-{d.area.max():g}",
                         integer_counts=round(d.integer.mean(), 3),
                         zero_fraction=round((d["count"] == 0).mean(), 3),
                         missing_fraction=round(d.density.isna().mean(), 3),
                         mean_count_per_unit=round(d["count"].mean(), 2),
                         units_per_visit_median=float(v["size"].median()),
                         multi_visit_station_years=round((vy > 1).mean(), 3),
                         nb_inverse_dispersion=round(inv_phi, 3),
                         visits_for_dispersion=len(vv)))
    return pd.DataFrame(rows)


def variance_scales():
    s = split_series(pd.read_csv(IR / "station_annual_values.csv"))
    s = s[s.Source != "Landsat"].dropna(subset=["value"])
    rows = []
    for (src, tax, resp), d in s.groupby(["Source", "Taxa", "Resp"]):
        pos = d.value[d.value > 0]
        if pos.empty:
            continue
        d = d.assign(lv=np.log(d.value + pos.min() / 2))
        sm = d.groupby(["series_id", "status", "site"]).lv.mean().reset_index()
        d["res"] = d.lv - d.groupby(["series_id", "status", "site"]).lv.transform("mean")
        d["res2"] = d.res - d.groupby(["series_id", "status", "year"]).res.transform("mean")
        rm = d.groupby(["series_id", "status", "year"]).lv.mean().reset_index()
        rm["dev"] = rm.lv - rm.groupby(["series_id", "status"]).lv.transform("mean")
        p = rm.pivot_table(index=["series_id", "year"], columns="status", values="lv")
        res, ac, cor = [], [], []
        if {"mpa", "reference"} <= set(p.columns):
            p = p.dropna()
            lr = p["mpa"] - p["reference"]
            for sid, g in lr.groupby(level=0):
                g = g.droplevel(0).sort_index()
                if len(g) < 6:
                    continue
                t = g.index.values.astype(float)
                X = np.c_[np.ones_like(t), t]
                e = g.values - X @ np.linalg.lstsq(X, g.values, rcond=None)[0]
                res.append(e.std(ddof=2))
                cons = np.diff(t) == 1
                if cons.sum() > 3 and e[:-1][cons].std() > 0 and e[1:][cons].std() > 0:
                    ac.append(np.corrcoef(e[:-1][cons], e[1:][cons])[0, 1])
                pp = p.loc[sid]
                if len(pp) > 5 and pp["mpa"].std() > 0 and pp["reference"].std() > 0:
                    cor.append(np.corrcoef(pp["mpa"], pp["reference"])[0, 1])
        med = lambda a: float(np.median(a)) if a else np.nan
        rows.append(dict(Source=src, Taxa=tax, Resp=resp, series=d.series_id.nunique(), stations=d.site.nunique(),
                         between_station_sd=sm.groupby(["series_id", "status"]).lv.std().median(),
                         within_station_residual_sd=d.res2.std(), role_year_sd=rm.dev.std(),
                         log_ratio_detrended_sd=med(res), log_ratio_lag1=med(ac),
                         protected_reference_year_correlation=med(cor)))
    return pd.DataFrame(rows).round(3)


def effect_scales():
    m = split_series(pd.read_csv(IR / "model_results.csv"))
    m = m[(m.rho == 0) & m.term.isin(["recovery_slope", "level_change"])]
    e = m.groupby(["model", "term", "Taxa", "Resp"]).agg(
        histories=("estimate", "size"), median=("estimate", "median"), sd=("estimate", "std"),
        q05=("estimate", lambda x: x.quantile(.05)), q95=("estimate", lambda x: x.quantile(.95)),
        median_se=("se", "median")).reset_index()
    p = pd.read_csv(IR / "pooled_results.csv")
    p = p[p.Source == "All programs"].assign(between_MPA_sd=lambda x: np.sqrt(x.between_MPA_variance))
    p = p[["Taxa", "Resp", "model", "term", "comparisons", "MPAs", "estimate", "between_MPA_sd", "boundary"]]
    el = split_series(pd.read_csv(IR / "model_eligibility.csv"))
    el = el.groupby(["Taxa", "Resp"]).agg(histories=("series_id", "size"), MPAs=("MPA", "nunique"),
                                          median_before_years=("before_years", "median"),
                                          histories_with_3plus_before=("before_years", lambda x: int((x >= 3).sum())),
                                          median_after_years=("after_years", "median")).reset_index()
    return e.round(3), p.round(3), el


def size_profile():
    inv = pd.read_csv(IR / "size_pool_inventory.csv")
    inv = inv[inv.Source.isin(["KFM", "PISCO", "LTER"])]
    inv["log_mean_mass"] = np.log(inv.mean_mass.where(inv.mean_mass > 0))
    a = inv.groupby(["Source", "Taxa"]).agg(
        pools=("pool_id", "size"), median_n=("n_sized", "median"),
        q10_n=("n_sized", lambda x: x.quantile(.1)), singleton_fraction=("singleton", "mean"),
        under10_fraction=("n_sized", lambda x: (x < 10).mean()),
        between_pool_sd_log_mean_mass=("log_mean_mass", "std")).reset_index()
    mem = pd.read_parquet(IR / "private/size_pool_members.parquet")
    mem = mem[mem.mass > 0].assign(lm=lambda x: np.log(x.mass))
    mem["n"] = mem.groupby("pool_id").frequency.transform("sum")
    mem = mem[mem.n >= 5]
    mem["mu"] = (mem.lm * mem.frequency).groupby(mem.pool_id).transform("sum") / mem.n
    w = mem.assign(sq=mem.frequency * (mem.lm - mem.mu) ** 2).groupby("pool_id").agg(sq=("sq", "sum"), n=("n", "first"))
    w["sd"] = np.sqrt(w.sq / (w.n - 1))
    w = w.reset_index()
    w["Source"] = w.pool_id.str.split("|").str[1]
    w["Taxa"] = w.pool_id.str.split("|").str[-1]
    wsd = w.groupby(["Source", "Taxa"]).sd.median().rename("median_within_pool_sd_log_mass").reset_index()
    a = a.merge(wsd, on=["Source", "Taxa"], how="left")
    u = pd.read_parquet(IR / "private/supported_measurement_uses.parquet")
    b = u[u.Resp == "Bio"]
    basis = pd.crosstab([b.Source, b.Taxa], b.size_basis).reset_index()
    return a.round(3), basis


def prior_predictive(seed=20261109, sims=4000):
    """Simulate log density for one MPA (1-6 protected and 1-6 reference stations, 20 years, designation at year 8)."""
    rng = np.random.default_rng(seed)
    hn = lambda s: np.abs(rng.normal(0, s))
    out = []
    years = np.arange(20)
    T0 = 8
    after = years >= T0
    for _ in range(sims):
        s_stat, s_reg, s_proc, s_obs = hn(1.0), hn(1.0), hn(1.0), hn(0.5)
        rho = np.clip(rng.normal(0.3, 0.3), -0.9, 0.95)
        mu_theta, tau_theta = rng.normal(0, 0.15), hn(0.1)
        mu_step, tau_step = rng.normal(0, 1.0), hn(0.75)
        theta = rng.normal(mu_theta, tau_theta)
        step = rng.normal(mu_step, tau_step)
        lam = rng.normal(0, hn(1.0))
        reg = rng.normal(0, s_reg, size=20)
        def station(protected):
            z = np.empty(20)
            z[0] = rng.normal(0, s_proc / np.sqrt(max(1 - rho ** 2, 0.05)))
            for t in range(1, 20):
                z[t] = rho * z[t - 1] + rng.normal(0, s_proc)
            base = rng.normal(lam if protected else 0, s_stat)
            eff = protected * after * (step + theta * (years - T0))
            return base + reg + z + eff + rng.normal(0, s_obs, size=20)
        n_p, n_r = rng.integers(1, 7), rng.integers(1, 7)
        P = np.mean([station(1) for _ in range(n_p)], axis=0)
        R = np.mean([station(0) for _ in range(n_r)], axis=0)
        lr = P - R
        X = np.c_[np.ones(20), years]
        e = lr - X @ np.linalg.lstsq(X, lr, rcond=None)[0]
        post = years >= T0
        Xp = np.c_[np.ones(post.sum()), years[post]]
        slope = np.linalg.lstsq(Xp, lr[post], rcond=None)[0][1]
        out.append(dict(log_ratio_detrended_sd=e.std(ddof=2), log_ratio_lag1=np.corrcoef(e[:-1], e[1:])[0, 1],
                        protected_reference_year_correlation=np.corrcoef(P - P.mean(), R - R.mean())[0, 1],
                        post_slope=slope, step=step, theta=theta))
    return pd.DataFrame(out)


def main():
    units = unit_profile(); units.to_csv(OUT / "unit_profile.csv", index=False)
    var = variance_scales(); var.to_csv(OUT / "variance_scales.csv", index=False)
    eff, pooled, elig = effect_scales()
    eff.to_csv(OUT / "effect_scales.csv", index=False)
    pooled.to_csv(OUT / "pooled_between_mpa_sd.csv", index=False)
    elig.to_csv(OUT / "history_support.csv", index=False)
    sizes, basis = size_profile()
    sizes.to_csv(OUT / "size_profile.csv", index=False)
    basis.to_csv(OUT / "biomass_size_basis_units.csv", index=False)
    pp = prior_predictive()
    q = pp.quantile([.05, .25, .5, .75, .95]).T.round(3)
    q.to_csv(OUT / "prior_predictive_quantiles.csv")
    obs = var[var.Resp == "Den"][["log_ratio_detrended_sd", "log_ratio_lag1", "protected_reference_year_correlation"]]
    cover = {c: float(((obs[c] >= pp[c].quantile(.05)) & (obs[c] <= pp[c].quantile(.95))).mean()) for c in obs}
    (OUT / "prior_predictive_check.json").write_text(json.dumps(
        {"observed_density_series_within_prior_90pct": cover,
         "observed_post_slope_q05_q95": [float(eff[eff.term == "recovery_slope"].q05.min()), float(eff[eff.term == "recovery_slope"].q95.max())],
         "prior_post_slope_q05_q95": [float(pp.post_slope.quantile(.05)), float(pp.post_slope.quantile(.95))]}, indent=2))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
