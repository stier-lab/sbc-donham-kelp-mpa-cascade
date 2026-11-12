"""Mean individual mass inputs for the joint model (design section 4.2; stage 3, experimental).

A "mass group" is a set of n animals (or plants) measured together whose mean mass is known:
  - a size pool (KFM and PISCO urchins, PISCO lobster, KFM kelp): animals measured at one site, date and zone.
    Taken from the integrated review's size_pool_inventory.csv (pool mean and number measured; no individual records);
  - one survey unit whose counted animals were all sized on the same record (sheephead; LTER lobster; PISCO kelp
    plants with stipe counts) or carry a recorded size-class mean (LTER urchins). Mean mass = biomass / density.
Only animals measured at the station itself are used. Borrowed and calibrated sizes are left out: the model's station
and region-year terms do the borrowing instead. LTER kelp is left out: LTER counts stipes and each stipe has a fixed
converted mass, so it carries no information on mass per plant.

Under a Gamma model for individual mass with shape k and mean M, the mean of n animals is Gamma(n k, mean M) exactly,
so group means with n are all the likelihood needs.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from data import age_basis

SAME_RECORD = {"same_record_fish_length", "same_record_lobster_length", "same_record_stipes", "same_quadrat_size_class_mean"}
ZMAP = {"INNER": 0, "MID": 1, "OUTER": 2, "DEEP": 3, "OUTER MIDDLE": 2, "INMID": 1, "OUTMID": 2}


def inventory_path(uses_path):
    p = Path(uses_path)
    for c in (p.parent / "size_pool_inventory.csv", p.parent.parent / "size_pool_inventory.csv"):
        if c.exists():
            return c
    raise FileNotFoundError("size_pool_inventory.csv not found beside or above " + str(p))


def find_beside(uses_path, name):
    """Private prepared inputs sit beside the survey-unit table (GRIT) or in the joint model's private folder (Mac)."""
    p = Path(uses_path)
    for c in (p.parent / name, p.parents[2] / "joint_model_20261109/private" / name if len(p.parents) > 2 else p.parent / name):
        if c.exists():
            return c
    return None


def mass_groups(uses_path, taxon, inv_path=None, variant="main", drop=()):
    """All mass groups for one taxon: station, year, zone, n, mean mass, kind.

    Preferred source: mass_groups.parquet from prep_mass.py (recomputed from source rows and checked against the
    integrated review; adds KFM sheephead roving-diver sizes, allometry-exponent variants and the >= 25 mm subset).
    variant: "main"; "b_lo" / "b_hi" (allometric exponent - 0.1 / + 0.1); "ge25" (urchins: only animals >= 25 mm
    test diameter in KFM and LTER groups, the size range PISCO measures).
    drop: programs, or "KFM-roving", whose mass groups are left out (sensitivities).
    Fallback when that file is absent: pool means from size_pool_inventory.csv and unit means from the survey-unit
    table (same values for every group both sources hold; no KFM sheephead, no variants)."""
    f = find_beside(uses_path, "mass_groups.parquet")
    if f is not None:
        g = pd.read_parquet(f); g = g[g.Taxa == taxon].copy()
        col = {"main": "mean", "ge25": "mean", "b_lo": "mean_blo", "b_hi": "mean_bhi"}[variant]
        g["mean"] = g[col]
        if variant == "ge25" and taxon in ("S. purpuratus", "M. franciscanus"):
            k = g.Source.isin(["KFM", "LTER"])
            g.loc[k, "n"] = g.loc[k, "n25"]; g.loc[k, "mean"] = g.loc[k, "mean25"]
        g = g.assign(station=g.Source + "|" + g.site, year=g.year.astype(int))
        for d in drop:
            g = g[~((g.Source == "KFM") & (g.kind == "roving"))] if d == "KFM-roving" else g[g.Source != d]
        g = g[["station", "year", "zone", "n", "mean", "kind"]]
        return g[(g.n > 0) & (g["mean"] > 0)], str(f.name)
    if variant != "main":
        raise FileNotFoundError("mass_groups.parquet (from prep_mass.py) is needed for variant " + variant)
    u = pd.read_parquet(uses_path)
    b = u[(u.Resp == "Bio") & (u.Taxa == taxon) & (u.Source != "Landsat")].drop_duplicates("unit_id")
    b = b[b.size_basis.isin(SAME_RECORD) & b.biomass.notna() & (b.density > 0)]
    if taxon == "M. pyrifera":
        b = b[b.Source != "LTER"]
    a = pd.DataFrame(dict(station=b.Source + "|" + b.site, year=b.year.astype(int), zone=b.zone, n=np.round(b.density * b.area),
                          mean=b.biomass / b.density, kind="unit"))
    inv = pd.read_csv(inv_path or inventory_path(uses_path))
    inv = inv[(inv.Taxa == taxon) & inv.Source.isin(["KFM", "PISCO"]) & (inv.n_sized > 0) & (inv.mean_mass > 0)]
    p = pd.DataFrame(dict(station=inv.Source + "|" + inv.site, year=inv.year.astype(int), zone=inv.zone, n=inv.n_sized.astype(float),
                          mean=inv.mean_mass, kind="pool"))
    g = pd.concat([a, p], ignore_index=True)
    for d in drop:
        g = g[~g.station.str.startswith(d + "|")]
    return g[(g.n > 0) & (g["mean"] > 0)], "size_pool_inventory.csv + survey-unit table"


def build_mass(uses_path, D, taxon, min_before=3, inv_path=None, variant="main", drop=()):
    """Mass inputs aligned to the density inputs D of the same taxon. Returns None if there are no groups."""
    g, source = mass_groups(uses_path, taxon, inv_path, variant, drop)
    s_ix = {s: i for i, s in enumerate(D["stations"])}
    y0, y1 = int(D["years"][0]), int(D["years"][-1])
    g = g[g.station.isin(s_ix) & (g.year >= y0) & (g.year <= y1)].copy()
    if len(g) == 0:
        return None
    g["sd"] = g.station.map(s_ix)                       # index into the density stations
    keep = np.array(sorted(g.sd.unique()))              # density-station indices that have mass data
    m_ix = {int(s): i for i, s in enumerate(keep)}
    g["s"] = g.sd.map(m_ix); g["t"] = g.year - y0
    T = len(D["years"]); yy = np.asarray(D["years"])[:, None]
    span = g.groupby("s").year.agg(["min", "max"]).sort_index()
    first = (yy == span["min"].values[None, :]).astype(float)
    inside = ((yy > span["min"].values[None, :]) & (yy <= span["max"].values[None, :])).astype(float)
    st_mpa = np.asarray(D["st_mpa"])[keep]; st_prot = np.asarray(D["st_prot"])[keep]
    T_st = np.asarray(D["mpa_T"])[st_mpa]
    # before-after status for mass: an MPA needs min_before pre-designation years with mass data in both roles
    sy = g.drop_duplicates(["s", "year"]); sy = sy.assign(mpa=st_mpa[sy.s.values], prot=st_prot[sy.s.values], T=T_st[sy.s.values])
    ba = np.zeros(len(D["mpas"]), dtype=int)
    for m in range(len(D["mpas"])):
        d = sy[(sy.mpa == m) & (sy.year < sy["T"])]
        ba[m] = int((d.groupby("year").prot.nunique() == 2).sum() >= min_before)
    prog = np.asarray(D["st_prog"])[keep]
    return dict(
        mass_stations=[D["stations"][i] for i in keep], dens_index=keep, n_groups=int(len(g)),
        g_s=g.s.values.astype(int), g_t=g.t.values.astype(int), g_n=g.n.values.astype(float), g_mean=g["mean"].values.astype(float),
        g_prog=prog[g.s.values], g_zone=g.zone.map(ZMAP).fillna(0).values.astype(int),
        st_prog=prog, st_mpa=st_mpa, st_prot=st_prot, st_region=np.asarray(D["st_region"])[keep],
        span_first=first, span_inside=inside, span=first + inside,
        grid_after=np.asarray(D["grid_after"])[:, keep], grid_age=np.asarray(D["grid_age"])[:, keep],
        grid_B=np.asarray(D["grid_B"])[:, keep, :], mpa_ba=ba,
        log_mean_mass=float(np.log(np.average(g["mean"], weights=g.n))),
        station_years_with_mass=int(len(sy)), density_station_years=int(len(set(zip(D["s"], D["t"])))),
        animals=float(g.n.sum()), kinds=g.kind.value_counts().to_dict(),
        source=source, variant=variant, dropped=list(drop),
        by_program={k: int(v) for k, v in g.station.str.split("|").str[0].value_counts().items()},
    )
