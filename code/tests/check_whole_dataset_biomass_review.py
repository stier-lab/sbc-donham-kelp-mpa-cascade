"""Independent gates for the review-only whole-dataset biomass outputs."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

root = Path(__file__).resolve().parents[2]
folder = root / "outputs/whole_dataset_biomass_update_20261105"
integrated = root / "outputs/integrated_review_20261103"
read = lambda name: pd.read_csv(folder / name)
annual = pd.read_csv(integrated / "annual_values.csv")
pooled = pd.read_csv(integrated / "pooled_results.csv")
scope = pd.read_csv(root / "config/integrated_comparison_scope.csv")
checks = json.loads((folder / "verification.json").read_text())
assert checks["paper_updated"] is False
assert checks["retained_histories"] == 175
assert checks["supported_comparison_years"] == int(annual.response_supported.sum())
assert int(scope.exclude.sum()) == 31
assert not set(annual.series_id) & set(scope.loc[scope.exclude, "series_id"])
assert not len(annual[annual.Source.eq("KFM") & annual.Taxa.eq("P. interruptus") & annual.Resp.eq("Bio")])
assert not len(annual[annual.Source.eq("KFM") & annual.Taxa.eq("S. pulcher") & annual.Resp.eq("Bio")])

effect = read("pooled_effects.csv")
for row in effect.itertuples():
    expected = pooled[(pooled.Taxa == row.Taxa) & (pooled.Resp == row.Resp)
                      & (pooled.Source == "All programs") & (pooled.model == row.model)
                      & (pooled.term == row.term)]
    assert len(expected) == 1 and np.isclose(row.estimate, expected.estimate.iloc[0])
    assert np.isclose(row.estimate_percent, 100 * np.expm1(row.estimate))

steps = read("series_change_steps.csv")
for before, delta, after in [
    ("old_all_slope", "old_coverage_step", "old_shared_slope"),
    ("old_shared_slope", "annual_input_step", "new_shared_slope"),
    ("new_shared_slope", "new_coverage_step", "new_all_slope"),
    ("new_all_slope", "joint_model_step", "joint_slope"),
]:
    z = steps[[before, delta, after]].dropna()
    assert np.allclose(z[before] + z[delta], z[after], atol=1e-10)
assert steps.annual_input_step_scope.str.contains("not biomass alone").all()

pairs = read("matched_year_pathway_annual.csv")
assert not pairs.duplicated(["panel", "MPA", "Source", "year"]).any()
links = {
    "a": ("P. interruptus", "Bio", "S. purpuratus", "Bio"),
    "b": ("P. interruptus", "Den", "S. purpuratus", "Den"),
    "c": ("S. pulcher", "Bio", "S. purpuratus", "Bio"),
    "d": ("S. pulcher", "Den", "S. purpuratus", "Den"),
    "e": ("S. purpuratus", "Den", "M. pyrifera", "Bio"),
    "f": ("M. franciscanus", "Den", "M. pyrifera", "Bio"),
}
for panel, (xt, xr, yt, yr) in links.items():
    p = pairs[pairs.panel.eq(panel)]
    x = annual[annual.Taxa.eq(xt) & annual.Resp.eq(xr)]
    y = annual[annual.Taxa.eq(yt) & annual.Resp.eq(yr)]
    z = p.merge(x, on=["MPA", "Source", "year"], validate="one_to_one")
    assert len(z) == len(p) and z.response_supported.all()
    assert np.allclose(z.x, z.ln_ratio)
    z = p.merge(y, on=["MPA", "Source", "year"], validate="one_to_one")
    assert len(z) == len(p) and z.response_supported.all()
    assert np.allclose(z.y, z.ln_ratio)
    q = x.merge(y, on=["MPA", "Source", "year"], suffixes=("_x", "_y"))
    q = q.merge(p[["MPA", "Source", "year"]], on=["MPA", "Source", "year"])
    assert q.protected_stations_x.eq(q.protected_stations_y).all()
    assert q.reference_stations_x.eq(q.reference_stations_y).all()
    assert q.shared_zones_x.eq(q.shared_zones_y).all()

fit = read("matched_year_pathway_fits.csv")
assert len(fit) == 18 and fit.status.eq("fitted").sum() == checks["conditional_pathway_fits"]
assert (fit.loc[fit.status.eq("fitted"), "fitted_MPAs"] >= 4).all()
assert not (fit.loc[fit.status.ne("fitted"), "p"].notna()).any()
heat = read("heatwave_state_pairs.csv")
assert not heat.duplicated(["MPA", "Source", "Taxa", "period"]).any()
assert np.allclose(heat.inside_ratio / heat.outside_ratio, np.exp(heat.log_inside_over_outside))
heat_fit = read("heatwave_period_refits.csv")
assert len(heat_fit) == 15 and heat_fit.status.eq("fitted").all()
assert (heat_fit.se > 0).all()

for relative, digest in json.loads((folder / "input_hashes.json").read_text()).items():
    assert hashlib.sha256((root / relative).read_bytes()).hexdigest() == digest, relative
for name in ["matched_year_pathway_screen.png", "matched_year_pathway_screen.pdf", "index.html", "README.md"]:
    assert (folder / name).is_file() and (folder / name).stat().st_size > 0
print("PASS whole-dataset biomass review gates")
