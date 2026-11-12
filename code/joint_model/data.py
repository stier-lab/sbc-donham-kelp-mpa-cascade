"""Build model inputs for the joint density model from the integrated review's survey-unit table.

Station selection, exclusions and roles come unchanged from the integrated review
(outputs/integrated_review_20261103/private/supported_measurement_uses.parquet).
"""
from pathlib import Path
import numpy as np
import pandas as pd

REGION = {  # proposed regions (docs/JOINT_MODEL_DESIGN.md section 9, decision 2; not yet confirmed by Adrian)
    "Anacapa Island SMR 2003": "N Channel Islands", "Scorpion SMR": "N Channel Islands", "Gull Island SMR": "N Channel Islands",
    "Harris Point SMR": "N Channel Islands", "South Point SMR": "N Channel Islands", "Carrington Pt SMR": "N Channel Islands",
    "Skunk Pt SMR": "N Channel Islands", "Judith Rk SMR": "N Channel Islands", "Santa Barbara Island SMR": "Santa Barbara Island",
    "Campus Point SMCA": "SB mainland", "Naples SMCA": "SB mainland", "Point Dume SMCA": "Malibu", "Point Dume SMR": "Malibu",
    "Abalone Cove SMCA": "Palos Verdes", "Point Vicente SMCA": "Palos Verdes", "Blue Cavern Onshore SMCA": "Catalina",
    "Long Point SMR": "Catalina", "Farnsworth Onshore SMCA": "Catalina", "Matlahuayl SMR": "San Diego",
    "South La Jolla SMR": "San Diego", "Cabrillo SMR": "San Diego", "Swamis SMCA": "San Diego"}
COARSE = {"three": {"N Channel Islands": "Northern islands", "Santa Barbara Island": "Southern islands", "Catalina": "Southern islands",
                    "SB mainland": "Mainland", "Malibu": "Mainland", "Palos Verdes": "Mainland", "San Diego": "Mainland"}}
KNOTS = np.array([0, 4, 8, 12, 16, 20], dtype=float)
ZONES = ["INNER", "MID", "OUTER", "DEEP"]


def build(uses_path, taxon, min_before=3, placebo_year=None, omit_mpas=(), omit_stations=(),
          space_placebo_seed=None, pisco_cap=None, regions="seven", size25=None, permute_seed=None):
    """Model inputs for one taxon.

    placebo_year: time placebo (KFM years before the real designation, fake designation year).
    omit_mpas / omit_stations: leave-out sensitivities ("KFM|a-k-01" style station names).
    space_placebo_seed: drop every protected station, then relabel one randomly chosen reference
        station (with at least five post-designation survey years) as protected in each MPA with at
        least three reference stations.
    pisco_cap: PISCO subsamples dense urchins in the field and extrapolates (delivered files do not
        record the area searched). With a cap, a PISCO unit count above it is treated as `cap`
        animals counted on a proportionally smaller area, which bounds the information in the count.
    """
    u = pd.read_parquet(uses_path)
    # Kelp has no density series in the review table: its diver counts sit on the biomass rows (KFM and PISCO count
    # plants, LTER counts stipes; program baselines absorb the difference in level). Exploratory use only.
    resp = "Bio" if taxon == "M. pyrifera" else "Den"
    u = u[(u.Resp == resp) & (u.Taxa == taxon) & (u.Source != "Landsat") & u.density.notna()].copy()
    u["MPA"] = u.series_id.str.split("|").str[0]
    if placebo_year is not None:
        # placebo in time: KFM years before the real designation only, with a fake designation year
        u = u[(u.Source == "KFM") & (u.year < u.MPA_Start)].copy()
        u = u[u.groupby("MPA").year.transform("min") <= placebo_year - min_before]
        u["MPA_Start"] = placebo_year
    u["zone"] = u.zone.replace({"OUTER MIDDLE": "OUTER", "INMID": "MID", "OUTMID": "OUTER"})
    u.loc[~u.zone.isin(ZONES), "zone"] = "none"
    u["station"] = u.Source + "|" + u.site
    n_before_omit = u.station.nunique()
    u = u[~u.MPA.isin(list(omit_mpas)) & ~u.station.isin(list(omit_stations))].copy()
    removed = n_before_omit - u.station.nunique()
    relabelled = []
    if permute_seed is not None:
        # Label permutation (randomization placebo): every station stays in the data. Within each MPA and the reference
        # stations assigned to it, the "protected" label is dealt out at random to the same number of stations as are
        # really protected. Unlike space_placebo_seed, successive seeds give nearly independent relabellings, so the
        # estimates from many seeds form a reference distribution for the real estimate under "labels do not matter".
        start = u.groupby("MPA").MPA_Start.first()
        hm = u.sort_values("MPA").groupby("station").MPA.first()
        real = u.groupby("station").status.first()
        u["MPA"] = u.station.map(hm); u["MPA_Start"] = u.MPA.map(start)
        rng = np.random.default_rng(permute_seed); new = {}
        for m in sorted(hm.unique()):
            sts = sorted(hm.index[hm == m]); k = int((real[sts] == "mpa").sum())
            if 0 < k < len(sts):
                pick = set(rng.choice(sts, size=k, replace=False).tolist())
                new.update({s_: ("mpa" if s_ in pick else "reference") for s_ in sts})
        u["status"] = u.station.map(new).fillna(u.status)
        relabelled = sorted(s_ for s_, v in new.items() if v == "mpa" and real[s_] != "mpa")
    if space_placebo_seed is not None:
        start = u.groupby("MPA").MPA_Start.first()
        u = u[u.status == "reference"].copy()
        hm = u.sort_values("MPA").groupby("station").MPA.first()
        u["MPA"] = u.station.map(hm); u["MPA_Start"] = u.MPA.map(start)
        rng = np.random.default_rng(space_placebo_seed)
        post_years = u[u.year >= u.MPA_Start].groupby("station").year.nunique()
        for m in sorted(hm.unique()):
            sts = sorted(hm.index[hm == m])
            ok = [s for s in sts if post_years.get(s, 0) >= 5]      # a fake protected station needs a usable history
            if len(sts) >= 3 and ok:
                pick = ok[int(rng.integers(len(ok)))]
                u.loc[u.station == pick, "status"] = "mpa"; relabelled.append(pick)
    # one MPA per station: protected stations have exactly one; shared references take the first alphabetically
    assert u[u.status == "mpa"].groupby("station").MPA.nunique().max() == 1
    assert u.groupby("station").status.nunique().max() == 1
    home = u.sort_values("MPA").groupby("station").agg(MPA=("MPA", "first"), status=("status", "first"),
                                                         Source=("Source", "first"), T=("MPA_Start", "first"))
    units = u.drop_duplicates("unit_id").copy()
    units["count"] = np.round(units.density * units.area)
    units["month"] = pd.to_datetime(units.date).dt.month
    capped = 0
    if pisco_cap:
        m = (units.Source == "PISCO") & (units["count"] > pisco_cap)
        units.loc[m, "area"] = units.loc[m, "area"] * pisco_cap / units.loc[m, "count"]
        units.loc[m, "count"] = pisco_cap; capped = int(m.sum())
    adjusted = 0
    if size25 is not None and taxon in ("S. purpuratus", "M. franciscanus"):
        # Common density definition (design 4.2, sensitivity): PISCO counts urchins >= 25 mm only, KFM and LTER count
        # all sizes. A KFM or LTER count of all sizes has mean (density >= 25 mm) * area / q, where q is the share of
        # measured animals >= 25 mm, so the unit's area is divided by q. q comes from the station-year when at least
        # 30 animals were measured there, otherwise from all years at that station, otherwise from the program.
        q = pd.read_csv(size25); q = q[q.Taxa == taxon].copy(); q["station"] = q.Source + "|" + q.site
        sy = q[q.n >= 30].set_index(["station", "year"]).share25
        stq = q.groupby("station")[["n25", "n"]].sum(); stq = (stq.n25 / stq.n)[stq.n >= 30]
        prq = q.groupby("Source")[["n25", "n"]].sum(); prq = prq.n25 / prq.n
        m = units.Source.isin(["KFM", "LTER"])
        share = pd.Series(pd.MultiIndex.from_frame(units[["station", "year"]]).map(sy), index=units.index)
        share = share.fillna(units.station.map(stq)).fillna(units.Source.map(prq)).clip(0.05, 1.0)
        units.loc[m, "area"] = units.loc[m, "area"] / share[m]; adjusted = int(m.sum())
    # cells = units sharing a visit, zone and unit area (identical expected count). The likelihood uses cell
    # totals for the mean and every unit count for the dispersion, which is exact for the negative binomial.
    key = ["station", "year", "date", "zone", "month", "area"]
    units["cell"] = units.groupby(key).ngroup()
    cells = units.groupby(key, as_index=False).agg(y=("count", "sum"), n=("unit_id", "size"), cell=("cell", "first"))
    assert (cells.cell.values == np.arange(len(cells))).all()
    cells = cells.join(home, on="station")
    stations = sorted(cells.station.unique()); s_ix = {s: i for i, s in enumerate(stations)}
    mpas = sorted(cells.MPA.unique()); m_ix = {m: i for i, m in enumerate(mpas)}
    programs = sorted(cells.Source.unique()); p_ix = {p: i for i, p in enumerate(programs)}
    # region sets: "seven" (proposed), "three" (northern islands, southern islands, mainland), "one"
    reg = {m: (REGION[m] if regions == "seven" else COARSE["three"][REGION[m]] if regions == "three" else "All") for m in REGION}
    regions = sorted({reg[m] for m in mpas}); r_ix = {r: i for i, r in enumerate(regions)}
    y0, y1 = int(cells.year.min()), int(cells.year.max()); years = np.arange(y0, y1 + 1)
    cells["visit"] = cells.station + "|" + cells.date.astype(str)
    visits = sorted(cells.visit.unique()); v_ix = {v: i for i, v in enumerate(visits)}
    st = home.loc[stations]
    # before-after status per MPA: >= min_before years before designation with both roles observed
    sy = cells.drop_duplicates(["station", "year"])
    ba = {}
    for m in mpas:
        d = sy[(sy.MPA == m) & (sy.year < sy["T"])]
        both = d.groupby("year").status.nunique()
        ba[m] = int((both == 2).sum() >= min_before)
    zone_ix = cells.zone.map({"none": 0, "INNER": 0, "MID": 1, "OUTER": 2, "DEEP": 3}).values
    age = (cells.year - cells["T"]).clip(lower=0).values.astype(float)
    Tn = len(years)
    span = cells.groupby("station").year.agg(["min", "max"]).loc[stations]
    yy = years[:, None]
    span_first = (yy == span["min"].values[None, :]).astype(float)
    span_inside = ((yy > span["min"].values[None, :]) & (yy <= span["max"].values[None, :])).astype(float)
    D = dict(
        span_first=span_first, span_inside=span_inside, stations_removed=int(removed), relabelled=relabelled, units_capped=capped, units_size25_adjusted=adjusted,
        unit_station=units.station.map(s_ix).values, unit_year=units.year.values.astype(int), unit_area=units.area.values.astype(float),
        taxon=taxon, stations=stations, mpas=mpas, programs=programs, regions=regions, years=years, knots=KNOTS,
        y=units["count"].values.astype(int), unit_cell=units.cell.values, unit_prog=units.Source.map(p_ix).values,
        log_area=np.log(cells.area.values), n_units=cells.n.values.astype(float),
        s=cells.station.map(s_ix).values, t=(cells.year - y0).values.astype(int), v=cells.visit.map(v_ix).values,
        zone=zone_ix, month=(cells.month - 1).values.astype(int),
        obs_prog=cells.Source.map(p_ix).values,
        st_prog=st.Source.map(p_ix).values, st_mpa=st.MPA.map(m_ix).values,
        st_region=st.MPA.map(reg).map(r_ix).values, st_prot=(st.status == "mpa").values.astype(float),
        after=(cells.year >= cells["T"]).values.astype(float) * (cells.status == "mpa").values.astype(float),
        prot=(cells.status == "mpa").values.astype(float), age=age,
        mpa_ba=np.array([ba[m] for m in mpas]), mpa_T=np.array([int(cells[cells.MPA == m]["T"].iloc[0]) for m in mpas]),
        visit_prog=np.array([p_ix[v.split("|")[0]] for v in visits]),
        log_mean_density=float(np.log(cells.y.sum() / cells.area.sum())),
        cells=cells[["station", "year", "date", "zone", "y", "area", "n", "MPA", "status", "Source", "T"]],
    )
    return D


def age_basis(age, knots=KNOTS):
    """Piecewise-linear basis: column j is the share of knot interval j completed at this age."""
    w = np.diff(knots)
    return np.clip((age[:, None] - knots[None, :-1]) / w[None, :], 0.0, 1.0)
