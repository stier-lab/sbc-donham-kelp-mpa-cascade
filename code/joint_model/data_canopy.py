"""Landsat kelp canopy as a response of its own (tag "can"): inputs for the joint model.

Source: canopy_quarterly.parquet from prep_canopy.py (quarterly canopy biomass summed over Tom Bell's corrected pixel
masks for each MPA and its one or two reference areas, 1984-2021). Cohort: the 18 comparisons frozen on 2026-11-01
(outputs/satellite_primary_20261101/CRITERIA.md); designation years from its cohort.csv. No site is added or dropped here.

What the canopy can and cannot say (Reed, Rassweiler, Nelson & Kui 2026, Front. Mar. Sci. 13:1845447): surface canopy is
12-38% of standing biomass; it tracks subsurface biomass well in February-May and poorly in autumn; it says little about
plant density; and a forest can hold much kelp under little or no canopy. So here:
  - canopy is its own response with its own protection effect, never pooled with diver counts of plants or stipes;
  - the main season is spring (January-June quarters); annual maximum, summer-autumn and all quarters are sensitivities;
  - a quarter with no canopy is "no canopy detected", modelled as a detection event whose probability rises with the
    latent canopy level, not as a biomass of zero.

One observation = one polygon in one quarter with at least 75% of its pixels observed (Tom Bell's coverage rule):
canopy density = biomass / (observed pixels x 900 m2), kg wet per m2.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from data import REGION, COARSE, KNOTS

# site number in Tom Bell's Site_List -> (name used in this project, designation year); the frozen 18
COHORT = {2: ("Harris Point SMR", 2003), 3: ("Carrington Pt SMR", 2003), 4: ("Skunk Pt SMR", 2003), 5: ("South Point SMR", 2003),
          6: ("Gull Island SMR", 2003), 8: ("Santa Barbara Island SMR", 2003), 10: ("Farnsworth Onshore SMCA", 2012),
          11: ("Swamis SMCA", 2012), 12: ("South La Jolla SMR", 2012), 13: ("Cabrillo SMR", 2012), 14: ("Point Conception SMR", 2012),
          15: ("Naples SMCA", 2012), 16: ("Campus Point SMCA", 2012), 17: ("Judith Rk SMR", 2003), 18: ("Point Dume SMCA", 2012),
          19: ("Point Dume SMR", 2012), 20: ("Abalone Cove SMCA", 2012), 21: ("Point Vicente SMCA", 2012)}
REGION_CAN = dict(REGION, **{"Point Conception SMR": "SB mainland"})     # proposed, like the others; Adrian to confirm
SEASONS = {"spring": (1, 2), "autumn": (3, 4), "all": (1, 2, 3, 4), "q2": (2,), "max": None}


def find_canopy(uses_path, name="canopy_quarterly.parquet"):
    """The table sits beside the uses parquet (GRIT scratch data/) or in outputs/joint_model_20261109/private/ (Mac)."""
    for p in (Path(uses_path).parent / name, Path("outputs/joint_model_20261109/private") / name, Path("data") / name):
        if p.exists():
            return str(p)
    return None


def build_canopy(path, season="spring", after="designation", depth_band="all", regions="three", permute_seed=None, placebo_year=None,
                 omit_mpas=(), min_cover=0.75, **_):
    """after: "designation" = the designation year counts as after (the joint model's convention for diver surveys);
              "next" = after starts the following year (the audited convention of the annual-maximum satellite analysis).
    permute_seed: randomization placebo, the protected label is given at random to one polygon of each comparison.
    placebo_year: time placebo, years before the real designation only, with this fake designation year for every MPA."""
    q = pd.read_parquet(path); q = q[(q.depth_band == depth_band) & q.site.isin(list(COHORT)) & (q.pixels > 0)].copy()
    q["MPA"] = q.site.map(lambda s: COHORT[s][0]); q["T"] = q.site.map(lambda s: COHORT[s][1])
    q = q[~q.MPA.isin(list(omit_mpas))]
    q["station"] = q.MPA + "|" + q.group
    q = q[q.pixels_observed >= min_cover * q.pixels].copy()
    q["dens"] = q.biomass_kg / (q.pixels_observed * 900.0)
    if season == "max":      # Tom Bell's annual maximum: needs at least three covered quarters
        n = q.groupby(["station", "year"]).quarter.transform("size")
        q = q[n >= 3].sort_values("dens").groupby(["station", "year"], as_index=False).last(); q["quarter"] = 0
    else:
        q = q[q.quarter.isin(SEASONS[season])].copy()
    prot = {s: g == "MPA" for s, g in q.drop_duplicates("station").set_index("station").group.items()}
    relabelled = []
    if permute_seed is not None:
        rng = np.random.default_rng(permute_seed)
        for m, d in q.drop_duplicates("station").groupby("MPA"):
            st = sorted(d.station); pick = st[rng.integers(len(st))]
            for s in st:
                prot[s] = s == pick
            relabelled.append(pick)
    if placebo_year is not None:
        q = q[q.year < q["T"] + (1 if after == "next" else 0)].copy(); q["T"] = int(placebo_year)
    stations = sorted(q.station.unique()); s_ix = {s: i for i, s in enumerate(stations)}
    mpas = sorted(q.MPA.unique()); m_ix = {m: i for i, m in enumerate(mpas)}
    reg = {m: (REGION_CAN[m] if regions == "seven" else COARSE["three"][REGION_CAN[m]] if regions == "three" else "All") for m in mpas}
    regs = sorted(set(reg.values())); r_ix = {r: i for i, r in enumerate(regs)}
    y0, y1 = int(q.year.min()), int(q.year.max()); years = np.arange(y0, y1 + 1)
    st = q.drop_duplicates("station").set_index("station").loc[stations]
    span = q.groupby("station").year.agg(["min", "max"]).loc[stations]; yy = years[:, None]
    quarters = [int(k) for k in sorted(q.quarter.unique())]; k_ix = {k: i for i, k in enumerate(quarters)}
    pos = (q.dens > 0).values
    D = dict(canopy=True, taxon="Kelp canopy (Landsat)", season=season, after=after, depth_band=depth_band,
             stations=stations, mpas=mpas, programs=["Landsat"], regions=regs, years=years, knots=KNOTS, quarters=quarters,
             st_prog=np.zeros(len(stations), dtype=int), st_mpa=st.MPA.map(m_ix).values, st_region=st.MPA.map(reg).map(r_ix).values,
             st_prot=np.array([float(prot[s]) for s in stations]), st_pixels=st.pixels.values.astype(float),
             mpa_ba=np.ones(len(mpas), dtype=int), mpa_T=np.array([int(q[q.MPA == m]["T"].iloc[0]) for m in mpas]),
             after_shift=1 if after == "next" else 0,
             span_first=(yy == span["min"].values[None, :]).astype(float),
             span_inside=((yy > span["min"].values[None, :]) & (yy <= span["max"].values[None, :])).astype(float),
             c_s=q.station.map(s_ix).values, c_t=(q.year - y0).values.astype(int), c_q=q.quarter.map(k_ix).values,
             c_logpix=np.log(q.pixels_observed.values.astype(float)), y=q.dens.values.astype(float),
             log_mean_density=float(np.log(q.dens[pos].mean())),
             c_centre=float(np.log(q.dens[pos].mean()) + np.log(q.pixels_observed.values.astype(float)).mean()),
             relabelled=relabelled, stations_removed=0,
             cells=q[["station", "year", "quarter", "MPA", "group", "T", "dens", "pixels", "pixels_observed"]])
    return D


def attach_canopy(D, y):
    """Observed canopy densities (kg wet per m2; 0 = no canopy detected) as model inputs."""
    y = np.asarray(y, dtype=float); pos = y > 0
    D.update(c_pos=pos.astype(float), c_logy=np.where(pos, np.log(np.where(pos, y, 1.0)), 0.0))
    return D
