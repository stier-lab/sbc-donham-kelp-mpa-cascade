"""Build the mean-mass inputs for the joint model from the integrated review's registered sources (run on the Mac).

  python3 code/joint_model/prep_mass.py

Writes (private, git-ignored; may go to GRIT scratch only):
  outputs/joint_model_20261109/private/mass_groups.parquet   one row per mass group (see data_mass.py)
  outputs/joint_model_20261109/private/size25_share.csv      share of measured urchins >= 25 mm by station-year
and a public-safe summary of the checks: outputs/joint_model_20261109/validation/prep_mass_checks.json

Each group carries its mean mass under the adopted allometry (mean), under exponent - 0.1 and + 0.1 (mean_blo,
mean_bhi), and for urchins the number and mean mass of animals >= 25 mm test diameter (n25, mean25), which is the size
range PISCO measures. Everything is recomputed from source rows and checked against the integrated review's own values.
"""
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "code/integrated_review"))
import common as C

OUTP = ROOT / "outputs/joint_model_20261109/private"; OUTP.mkdir(parents=True, exist_ok=True)
REV = ROOT / "outputs/integrated_review_20261103"
COEF = {"M. franciscanus": (.00059, 2.917), "S. purpuratus": (.00059, 2.870), "P. interruptus": (.001352821, 2.913963113), "S. pulcher": (.0144, 3.04)}
URCHINS = ["M. franciscanus", "S. purpuratus"]
checks = {}


def mvar(tx, size, db):
    a, b = COEF[tx]
    return a * np.asarray(size, dtype=float) ** (b + db)


for tx, (a, b) in COEF.items():          # the coefficients here must be the review's
    assert np.allclose(C.mass(tx, np.array([10., 50., 120.])), a * np.array([10., 50., 120.]) ** b)

rows = []

def add(tx, source, site, year, date, zone, kind, size, count, thr=None):
    """One group from individual sizes (same units as the review passes to mass()) and their frequencies."""
    size = np.asarray(size, dtype=float); count = np.asarray(count, dtype=float); n = count.sum()
    r = dict(Taxa=tx, Source=source, site=site, year=int(year), date=date, zone=zone, kind=kind, n=float(n),
             mean=float(np.average(mvar(tx, size, 0), weights=count)), mean_blo=float(np.average(mvar(tx, size, -.1), weights=count)),
             mean_bhi=float(np.average(mvar(tx, size, .1), weights=count)), n25=np.nan, mean25=np.nan)
    if thr is not None:
        k = size >= thr; r["n25"] = float(count[k].sum())
        r["mean25"] = float(np.average(mvar(tx, size[k], 0), weights=count[k])) if k.any() else np.nan
    rows.append(r)

# 1. size pools: same construction as code/integrated_review/biomass.py
s = C.read_tracked(C.DATA / "data/processed/source_build/analysis_sizes.csv")
s["Taxa"] = s.classcode.map(C.TAXA); s["Source"] = s.campus.map({"UCSB": "PISCO", "VRG": "PISCO", "KFM": "KFM", "LTER": "LTER"})
s = s[s.Source.isin(["KFM", "PISCO"]) & s.Taxa.notna() & s["size"].gt(0) & s["count"].gt(0)].copy()
s = s[~((s.Source == "PISCO") & s.Taxa.isin(URCHINS) & s["size"].lt(25))]
s["zone"] = s.zone.fillna("not_recorded")
s["date"] = pd.to_datetime(s.date, errors="coerce").fillna(pd.to_datetime(s[["year", "month", "day"]], errors="coerce")).dt.strftime("%Y-%m-%d")
s = s[s.date.notna()]
s = s.loc[~pd.Series([C.blocked(x, y, t) for x, y, t in zip(s.Source, s.site, s.Taxa)], index=s.index)]
for (source, site, year, date, zone, tx), g in s.groupby(["Source", "site", "year", "date", "zone", "Taxa"]):
    add(tx, source, site, year, date, zone, "pool", g["size"], g["count"], thr=25 if tx in URCHINS else None)
inv = pd.read_csv(REV / "size_pool_inventory.csv")
p = pd.DataFrame(rows); p["pool_id"] = ["|".join(map(str, ("size", r.Source, r.site, r.year, r.date, r.zone, r.Taxa))) for r in p.itertuples()]
m = p.merge(inv[["pool_id", "mean_mass", "n_sized"]], on="pool_id", how="outer", indicator=True)
anim = inv[inv.Taxa.isin(COEF) & inv.Source.isin(["KFM", "PISCO"])]
mm = m[m.pool_id.isin(anim.pool_id) | m._merge.ne("right_only")]
checks["pools"] = dict(pools=int(len(p)), inventory_animal_pools=int(len(anim)), unmatched=int((mm._merge != "both").sum()),
                       max_rel_diff_mean=float((mm["mean"] / mm.mean_mass - 1).abs().max()), max_diff_n=float((mm.n - mm.n_sized).abs().max()))
assert checks["pools"]["unmatched"] == 0 and checks["pools"]["max_rel_diff_mean"] < 1e-9 and checks["pools"]["max_diff_n"] < 1e-9, checks
# KFM kelp pools have no allometry (stipes x monthly conversion): copy from the inventory
kk = inv[(inv.Taxa == "M. pyrifera") & (inv.Source == "KFM") & (inv.n_sized > 0) & (inv.mean_mass > 0)]
for r in kk.itertuples():
    rows.append(dict(Taxa=r.Taxa, Source="KFM", site=r.site, year=int(r.year), date=r.date, zone=r.zone, kind="pool", n=float(r.n_sized),
                     mean=r.mean_mass, mean_blo=r.mean_mass, mean_bhi=r.mean_mass, n25=np.nan, mean25=np.nan))

# 2. survey units whose counted animals carry sizes on the same source rows
u = pd.read_parquet(REV / "private/supported_measurement_uses.parquet")
u = u[(u.Resp == "Bio") & u.biomass.notna() & (u.density > 0)].drop_duplicates(["unit_id", "Taxa"])
SPEC = {"pisco_fish": ("fish_tl", "count", 1.0, dict(encoding="latin1")), "lter_fish": ("SIZE", "COUNT", 1.0, {}),
        "lter_lobster": ("SIZE_MM", "COUNT", 1.0, {}), "lter_benthic_19_33": ("SIZE", "COUNT", 10.0, {})}
BASIS = {"same_record_fish_length", "same_record_lobster_length", "same_quadrat_size_class_mean"}
worst = {}
for sid_, (scol, ccol, mult, kw) in SPEC.items():
    src = C.read_source(sid_, **kw); src["_row"] = np.arange(1, len(src) + 1); src = src.set_index("_row")
    uu = u[(u.source_id == sid_) & u.size_basis.isin(BASIS)]
    d = []
    for r in uu.itertuples():
        g = src.loc[[int(x) for x in str(r.source_rows).split(";")]]
        if sid_ == "lter_benthic_19_33":
            g = g[g.SP_CODE.map(C.TAXA) == r.Taxa]
        c = pd.to_numeric(g[ccol], errors="coerce").to_numpy(); z = pd.to_numeric(g[scol], errors="coerce").to_numpy() * mult
        ok = (c > 0) & (z > 0)
        if not ok.any():
            continue
        add(r.Taxa, r.Source, r.site, r.year, r.date, r.zone, "unit", z[ok], c[ok], thr=25 if r.Taxa in URCHINS else None)
        d.append((abs(rows[-1]["mean"] / (r.biomass / r.density) - 1), abs(rows[-1]["n"] - r.density * r.area)))
    d = np.array(d); worst[sid_] = dict(units=int(len(d)), max_rel_diff_mean=float(d[:, 0].max()), max_diff_n=float(d[:, 1].max()))
    assert d[:, 0].max() < 1e-8 and d[:, 1].max() < 1e-6, (sid_, worst[sid_])
checks["same_record_units"] = worst
# PISCO kelp plants (stipes x monthly conversion; no allometry): unit mean = biomass / density
k = u[(u.Taxa == "M. pyrifera") & (u.Source == "PISCO") & (u.size_basis == "same_record_stipes")]
for r in k.itertuples():
    mu = r.biomass / r.density
    rows.append(dict(Taxa=r.Taxa, Source="PISCO", site=r.site, year=int(r.year), date=r.date, zone=r.zone, kind="unit",
                     n=float(np.round(r.density * r.area)), mean=mu, mean_blo=mu, mean_bhi=mu, n25=np.nan, mean25=np.nan))

# 3. KFM sheephead: roving-diver size surveys (Scott Gabara's 2026-11-02 export), all KFM sites, years through 2023
paths = [Path(x["path"]) for x in json.loads((ROOT / "outputs/kfm_sheephead_biomass_20261103/inputs.json").read_text())]
for x in json.loads((ROOT / "outputs/kfm_sheephead_biomass_20261103/inputs.json").read_text()):
    assert C.sha(x["path"]) == x["sha256"], "changed input " + x["path"]
f = pd.read_csv(paths[1]); f = f[f.scientificName.str.contains("pulcher", na=False)].copy()
f["site"] = f.siteName.map(pd.read_csv(paths[3]).set_index("source_name").site)
f = f[f.site.notna() & f.totalLengthInCentimeters.gt(0) & f.speciesCount.gt(0)]
f["year"] = f.startDate.str[:4].astype(int); f = f[f.year <= 2023]
n0 = len(rows)
for (site, date), g in f.groupby(["site", "startDate"]):
    add("S. pulcher", "KFM", site, int(date[:4]), date[:10], "not_recorded", "roving", g.totalLengthInCentimeters, g.speciesCount)
kf = pd.DataFrame(rows[n0:])
vs = pd.read_csv(ROOT / "outputs/kfm_sheephead_biomass_20261103/visit_support.csv"); vs = vs[vs["mode"] == "exact"]
j = vs.merge(kf, left_on=["site_id", "size_date"], right_on=["site", "date"], how="left")
checks["kfm_sheephead"] = dict(size_surveys=int(len(kf)), sites=int(kf.site.nunique()), fish=float(kf.n.sum()), review_exact_visits=int(len(j)),
                               unmatched=int(j["mean"].isna().sum()), max_rel_diff_mean=float((j["mean"] / j.mean_mass_g - 1).abs().max()),
                               max_diff_n=float((j.n - j.n_fish).abs().max()))
assert checks["kfm_sheephead"]["unmatched"] == 0 and checks["kfm_sheephead"]["max_rel_diff_mean"] < 1e-9, checks["kfm_sheephead"]

g = pd.DataFrame(rows)
assert g["mean"].gt(0).all() and g.n.gt(0).all()
assert not any(C.blocked(a, b, t) for a, b, t in zip(g.Source, g.site, g.Taxa))
g.to_parquet(OUTP / "mass_groups.parquet", index=False)

# 4. share of measured urchins >= 25 mm by station-year (KFM pools, LTER quadrat groups), for a common density definition
q = g[g.Taxa.isin(URCHINS) & g.Source.isin(["KFM", "LTER"])].groupby(["Taxa", "Source", "site", "year"], as_index=False).agg(n=("n", "sum"), n25=("n25", "sum"))
q["share25"] = q.n25 / q.n
q.to_csv(OUTP / "size25_share.csv", index=False)

checks["groups"] = {f"{t}|{s}|{k}": dict(groups=int(len(x)), animals=float(x.n.sum()), station_years=int(x.drop_duplicates(["site", "year"]).shape[0]))
                    for (t, s, k), x in g.groupby(["Taxa", "Source", "kind"])}
z = g[g.Taxa.isin(COEF)]
checks["exponent_shift_log_mean"] = {t: dict(minus=float(np.average(np.log(x.mean_blo / x["mean"]), weights=x.n)), plus=float(np.average(np.log(x.mean_bhi / x["mean"]), weights=x.n)))
                                     for t, x in z.groupby("Taxa")}
checks["share25"] = {f"{t}|{s}": dict(station_years=int(len(x)), animal_weighted_share=float(x.n25.sum() / x.n.sum()), q10=float(x.share25.quantile(.1)), median=float(x.share25.median()))
                     for (t, s), x in q.groupby(["Taxa", "Source"])}
(ROOT / "outputs/joint_model_20261109/validation/prep_mass_checks.json").write_text(json.dumps(checks, indent=1))
print(json.dumps(checks, indent=1))
