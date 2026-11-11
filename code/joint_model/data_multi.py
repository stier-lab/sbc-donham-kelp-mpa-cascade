"""Inputs for the multi-taxon model (stage 2): one single-taxon build per taxon, plus the index maps that link taxa.

Each taxon keeps its own stations, years and MPAs exactly as data.build() returns them. A link (edge) from a source
taxon to a target taxon is a map from every target station-year to the same station in the source taxon `lag` years
earlier, with a mask that is 1 only where the source taxon has a latent state there (inside its surveyed span).
"""
import numpy as np
from data import build, age_basis

TAXA = {"lob": "P. interruptus", "shp": "S. pulcher", "pur": "S. purpuratus", "red": "M. franciscanus", "kelp": "M. pyrifera"}
# target <- source. Predators last year -> urchins this year; urchins last year -> kelp this year.
EDGES = [("pur", "lob"), ("pur", "shp"), ("red", "lob"), ("red", "shp"), ("kelp", "pur"), ("kelp", "red")]
# Sampler parameterisation that converged for each taxon in the single-taxon checks (REPORT_recommended_structure.md).
CONFIG = {"lob": dict(centered=False, level_dev=True), "shp": dict(centered=False, level_dev=True),
          "pur": dict(centered=True, level_dev=False), "red": dict(centered=True, level_dev=False),
          "kelp": dict(centered=True, level_dev=False)}


def grid(D):
    """Protection function inputs on the full year x station grid (the single-taxon model has them per cell only)."""
    T_st = np.asarray(D["mpa_T"])[D["st_mpa"]]
    yy = np.asarray(D["years"])[:, None]
    D["grid_after"] = ((yy >= T_st[None, :]) * D["st_prot"][None, :]).astype(float)
    D["grid_age"] = np.clip(yy - T_st[None, :], 0, None).astype(float)
    D["grid_B"] = age_basis(D["grid_age"].ravel()).reshape(D["grid_age"].shape + (-1,))
    D["span"] = D["span_first"] + D["span_inside"]
    return D


def link(Dk, Dj, lag):
    """Target k <- source j: source year index, source station index and validity mask on the target's grid."""
    Tk, Sk = len(Dk["years"]), len(Dk["stations"])
    j_ix = {s: i for i, s in enumerate(Dj["stations"])}
    src_s = np.array([j_ix.get(s, -1) for s in Dk["stations"]])
    src_t = np.arange(Tk) + int(Dk["years"][0]) - int(Dj["years"][0]) - lag
    ok_t = (src_t >= 0) & (src_t < len(Dj["years"]))
    st, ss = np.clip(src_t, 0, len(Dj["years"]) - 1), np.clip(src_s, 0, None)
    valid = ok_t[:, None] * (src_s >= 0)[None, :] * Dj["span"][st][:, ss]
    return dict(src_t=st, src_s=ss, valid=valid.astype(float),
                n_target_station_years=int(Dk["span"].sum()), n_linked=int((valid * Dk["span"]).sum()),
                n_stations_linked=int(((valid * Dk["span"]).sum(0) > 0).sum()))


def build_multi(uses_path, tags=("lob", "shp", "pur", "red", "kelp"), edges=None, lag=1, regions="three", **kw):
    """Returns (data, edges): data[tag] is the single-taxon input plus grid arrays; edges is a list of
    (target, source, link) in an order in which every source is computed before its target."""
    data = {g: grid(build(uses_path, TAXA[g], regions=regions, **kw)) for g in tags}
    edges = [e for e in (EDGES if edges is None else edges) if e[0] in tags and e[1] in tags]
    order = list(tags)
    assert all(order.index(j) < order.index(k) for k, j in edges), "list source taxa before their targets"
    return data, [(k, j, link(data[k], data[j], lag)) for k, j in edges]
