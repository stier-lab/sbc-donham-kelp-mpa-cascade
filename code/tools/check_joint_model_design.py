"""Design checks for the joint hierarchical model (v2), from committed integrated-review annual values.

Run from the analysis repository root: python3 code/tools/check_joint_model_design.py
Reads outputs/integrated_review_20261103/annual_values.csv; writes outputs/joint_model_design_20261109/design_checks.json.
Descriptive only. Detrending biases autocorrelations toward zero in short series.
"""
import itertools, json, warnings
from pathlib import Path
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = Path.cwd()
a = pd.read_csv(ROOT / "outputs/integrated_review_20261103/annual_values.csv")
a[["MPA", "Taxa", "Source", "Resp"]] = a.series_id.str.split("|", expand=True)
a = a[(a.Source != "Landsat") & a.ln_ratio.notna() & (a.Resp == "Den")]
out = {}

# A. lag-1..3 autocorrelation of detrended protected/reference log ratios (series with >=15 years)
rows = []
for sid, g in a.groupby("series_id"):
    g = g.sort_values("year")
    if len(g) < 15:
        continue
    t = g.year.values.astype(float)
    X = np.c_[np.ones_like(t), t, (t >= g.MPA_Start.iloc[0]).astype(float)]
    s = pd.Series(g.ln_ratio.values - X @ np.linalg.lstsq(X, g.ln_ratio.values, rcond=None)[0], index=g.year.values)
    def ac(k):
        x = s.reindex(s.index + k); m = x.notna().values
        return np.corrcoef(s.values[m], x.values[m])[0, 1] if m.sum() > 8 else np.nan
    rows.append(dict(Taxa=g.Taxa.iloc[0], Source=g.Source.iloc[0], r1=ac(1), r2=ac(2), r3=ac(3)))
r = pd.DataFrame(rows)
out["autocorrelation_by_taxon"] = r.groupby("Taxa")[["r1", "r2", "r3"]].median().round(2).to_dict("index")
out["autocorrelation_by_program"] = r.groupby("Source")[["r1", "r2", "r3"]].median().round(2).to_dict("index")
out["series_used"] = int(len(r))

# B. correlation between MPAs of detrended post-protection residuals (same taxon, different MPA)
res = {}
for sid, g in a.groupby("series_id"):
    g = g[g.year >= g.MPA_Start].sort_values("year")
    if len(g) < 8:
        continue
    t = g.year.values.astype(float); X = np.c_[np.ones_like(t), t]
    res[sid] = pd.Series(g.ln_ratio.values - X @ np.linalg.lstsq(X, g.ln_ratio.values, rcond=None)[0], index=g.year.values)
pairs = []
for s1, s2 in itertools.combinations(res, 2):
    m1, t1, _, _ = s1.split("|"); m2, t2, _, _ = s2.split("|")
    if t1 != t2 or m1 == m2:
        continue
    j = pd.concat([res[s1], res[s2]], axis=1).dropna()
    if len(j) >= 8:
        pairs.append(dict(Taxa=t1, r=np.corrcoef(j.iloc[:, 0], j.iloc[:, 1])[0, 1]))
p = pd.DataFrame(pairs)
out["cross_mpa_residual_correlation"] = p.groupby("Taxa").r.agg(["size", "mean"]).round(3).to_dict("index")

# C. shape: mean change in log ratio from ages 0-3, by age since designation
rows = []
for sid, g in a.groupby("series_id"):
    g = g[g.year >= g.MPA_Start].assign(age=lambda x: x.year - x.MPA_Start)
    base = g[g.age <= 3].ln_ratio.mean()
    if len(g) < 8 or np.isnan(base):
        continue
    rows.append(g.assign(d=g.ln_ratio - base, bin=pd.cut(g.age, [-1, 3, 7, 11, 15, 21],
                labels=["0-3", "4-7", "8-11", "12-15", "16-20"]))[["Taxa", "bin", "d", "series_id"]])
d = pd.concat(rows)
shape = d.groupby(["Taxa", "bin"], observed=True).agg(change=("d", "mean"), histories=("series_id", "nunique")).round(2)
out["change_by_age"] = {t: shape.loc[t].to_dict("index") for t in shape.index.get_level_values(0).unique()}

# D. placebo feasibility: KFM before-protection years per density history
k = a[(a.Source == "KFM") & (a.year < a.MPA_Start)]
out["kfm_before_years"] = k.groupby("series_id").year.agg(["min", "max", "size"]).astype(int).to_dict("index")

dest = ROOT / "outputs/joint_model_design_20261109"
dest.mkdir(parents=True, exist_ok=True)
(dest / "design_checks.json").write_text(json.dumps(out, indent=2))
print(json.dumps({k: out[k] for k in ("autocorrelation_by_taxon", "cross_mpa_residual_correlation")}, indent=1))
