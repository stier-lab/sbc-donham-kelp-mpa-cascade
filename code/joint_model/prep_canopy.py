"""Quarterly Landsat kelp canopy by MPA and reference polygon, rebuilt from the pixel file with Tom Bell's corrected masks.

  python3 code/joint_model/prep_canopy.py            (analysis repo root; needs the local delivery and the cached NetCDF)

Inputs (all local, none tracked): data/external/tom_bell_2026-10-21/{MPA_Ref_Locations_v2.xlsx, MPA_Runs_new.xlsx},
data/external/tom_bell_2026-10-20/03_Site_List.csv, local_archive/metadata_audit_2026-10-16/CAkelpCanopyEnv_2021.nc
(LANDSAT_V2_DIR and LANDSAT_NETCDF override the first and last).

Writes outputs/joint_model_20261109/private/canopy_quarterly.parquet (one row per polygon, quarter and depth band:
pixels, pixels observed, pixels with canopy, canopy biomass kg wet, canopy area m2) and
outputs/joint_model_20261109/validation/prep_canopy_checks.json.

Check: with Tom's coverage rules (a quarter is missing if more than 25% of the polygon's pixels are missing; the annual
value is the maximum of the quarters and missing if more than one quarter is missing) the table must reproduce the
delivered annual arrays for every non-empty polygon. The script stops if it does not.
"""
import json, os
from pathlib import Path
import numpy as np
import pandas as pd
import netCDF4
import openpyxl

V2 = Path(os.environ.get("LANDSAT_V2_DIR", "data/external/tom_bell_2026-10-21"))
NC = os.environ.get("LANDSAT_NETCDF", "local_archive/metadata_audit_2026-10-16/CAkelpCanopyEnv_2021.nc")
SITES = Path("data/external/tom_bell_2026-10-20/03_Site_List.csv")
OUT = Path("outputs/joint_model_20261109")
(OUT / "private").mkdir(parents=True, exist_ok=True); (OUT / "validation").mkdir(parents=True, exist_ok=True)

nc = netCDF4.Dataset(NC); nc.set_auto_mask(False)
lat, lon, depth = nc["lat"][:], nc["lon"][:], nc["depth"][:]
qy, qq = nc["year"][:].astype(int), nc["quarter"][:].astype(int)
key = lambda a, b: [f"{x:.6f}|{y:.6f}" for x, y in zip(np.round(a, 6), np.round(b, 6))]
keys = {k: i for i, k in enumerate(key(lat, lon))}
assert len(keys) == len(lat), "duplicate pixel coordinates"
years = np.arange(1984, 2022); assert sorted(set(qy)) == list(years) and all((qy == y).sum() == 4 for y in years)

ws = openpyxl.load_workbook(V2 / "MPA_Ref_Locations_v2.xlsx", read_only=True, data_only=True)["Pixel Locations"]
rows = [list(r) for r in ws.iter_rows(values_only=True)]
ncol = max(len(r) for r in rows); rows = [r + [None] * (ncol - len(r)) for r in rows]
num = lambda x: float(x) if isinstance(x, (int, float)) and np.isfinite(x) else np.nan
groups = []
for j in range(1, ncol, 2):              # column pairs (lat, lon); row 1 = site number, row 2 = MPA / Ref_1 / Ref_2
    if rows[0][j] is None:
        continue
    a = np.array([num(r[j]) for r in rows[2:]]); b = np.array([num(r[j + 1]) for r in rows[2:]])
    assert (np.isfinite(a) == np.isfinite(b)).all()
    idx = np.array([keys[k] for k in key(a[np.isfinite(a)], b[np.isfinite(b)])], dtype=int)
    assert len(set(idx)) == len(idx), "repeated pixel in one polygon"
    groups.append(dict(site=int(rows[0][j]), group=str(rows[1][j]), idx=idx))
assert len(groups) == 72, len(groups)

names = pd.read_csv(SITES, quotechar="'"); names.columns = ["site", "name", "lat", "lon"]
name = dict(zip(names.site, names.name))
BANDS = {"all": (-np.inf, np.inf), "to10m": (-np.inf, 10.0), "over10m": (10.0, np.inf)}
dabs = np.abs(depth)
rec = []; q_all = np.full((72, len(qy)), np.nan)
for t in range(len(qy)):
    bio = nc["biomass"][t, :].astype(float); bio[bio < 0] = np.nan
    area = nc["area"][t, :].astype(float); area[area < 0] = np.nan
    for j, g in enumerate(groups):
        if not len(g["idx"]):
            continue
        for band, (lo, hi) in BANDS.items():
            ix = g["idx"][(dabs[g["idx"]] > lo) & (dabs[g["idx"]] <= hi)] if band != "all" else g["idx"]
            if not len(ix):
                continue
            x = bio[ix]; obs = np.isfinite(x)
            rec.append((g["site"], g["group"], int(qy[t]), int(qq[t]), band, len(ix), int(obs.sum()), int((x[obs] > 0).sum()),
                        float(np.nansum(x)), float(np.nansum(area[ix]))))
        x = bio[g["idx"]]
        if np.isnan(x).mean() <= 0.25:
            q_all[j, t] = np.nansum(x)
tab = pd.DataFrame(rec, columns=["site", "group", "year", "quarter", "depth_band", "pixels", "pixels_observed", "pixels_with_canopy", "biomass_kg", "area_m2"])
tab.insert(1, "name", tab.site.map(name))
tab.to_parquet(OUT / "private" / "canopy_quarterly.parquet", index=False)

# ---- check against the delivered annual arrays
wb = openpyxl.load_workbook(V2 / "MPA_Runs_new.xlsx", read_only=True, data_only=True)
deliv = {}
for g, s in (("MPA", "MPAs"), ("Ref_1", "MPAC1"), ("Ref_2", "MPAC2")):
    m = np.array([[num(c) for c in r] for r in wb[s].iter_rows(values_only=True)], dtype=float)
    assert m.shape == (24, 38), (s, m.shape); deliv[g] = m
annual = np.full((72, 38), np.nan)
for k, y in enumerate(years):
    x = q_all[:, qy == y]; okq = np.isfinite(x).sum(1) >= 3
    annual[okq, k] = np.nanmax(x[okq], axis=1)
worst, n_series, empty = 0.0, 0, 0
for j, g in enumerate(groups):
    b = deliv[g["group"]][g["site"] - 1]
    if not len(g["idx"]):
        empty += 1; continue
    a = annual[j]; assert (np.isfinite(a) == np.isfinite(b)).all(), ("missing pattern", g["site"], g["group"])
    f = np.isfinite(a); worst = max(worst, float(np.abs(a[f] - b[f]).max()) if f.any() else 0.0); n_series += 1
assert worst < 1e-6, worst

# ---- descriptive facts used in the plan (all polygons with pixels; whole polygon)
a = tab[tab.depth_band == "all"].copy(); a["covered"] = a.pixels_observed >= 0.75 * a.pixels
c = a[a.covered]
zero_by_q = c.groupby("quarter").apply(lambda d: float((d.biomass_kg == 0).mean())).round(3).to_dict()
miss_by_q = a.groupby("quarter").apply(lambda d: float((~d.covered).mean())).round(3).to_dict()
mx = c.loc[c.groupby(["site", "group", "year"]).biomass_kg.idxmax()]; mx = mx[mx.biomass_kg > 0]
max_q = (mx.quarter.value_counts(normalize=True).sort_index().round(3)).to_dict()
dep = {f"{g['site']}|{g['group']}": dict(pixels=int(len(g["idx"])), median_depth_m=float(np.median(dabs[g["idx"]])), share_to_10m=float((dabs[g["idx"]] <= 10).mean()))
       for g in groups if len(g["idx"])}
chk = dict(polygons=72, non_empty=n_series, empty=empty, annual_reconstruction_max_abs_difference=worst, quarters=int(len(qy)),
           share_of_covered_polygon_quarters_with_zero_canopy=zero_by_q, share_of_polygon_quarters_failing_25pct_coverage=miss_by_q,
           quarter_of_annual_maximum=max_q, depth_sign_note=f"depth min {float(depth.min()):.1f}, max {float(depth.max()):.1f} (absolute value used)",
           depth_by_polygon=dep, rows=int(len(tab)))
(OUT / "validation" / "prep_canopy_checks.json").write_text(json.dumps(chk, indent=1))
print(json.dumps({k: v for k, v in chk.items() if k != "depth_by_polygon"}, indent=1))
