"""Review-only downstream summaries from the integrated source-to-model build.

This module never reads or writes the manuscript's frozen effect sizes. It keeps
post-protection slopes, before-after changes, pathway screens, and heatwave
state comparisons separate because they answer different questions.
"""
from __future__ import annotations

from pathlib import Path
import hashlib
import html
import json
import os
import subprocess

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
IN = ROOT / "outputs/integrated_review_20261103"
OUT = ROOT / "outputs/whole_dataset_biomass_update_20261105"
KFM = ROOT / "outputs/kfm_sheephead_biomass_20261103"
OUT.mkdir(parents=True, exist_ok=True)


def read(name: str) -> pd.DataFrame:
    return pd.read_csv(IN / name)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def saved(frame: pd.DataFrame, name: str) -> None:
    frame.to_csv(OUT / name, index=False)


complete = json.loads((IN / "completed.json").read_text())
verify = json.loads((IN / "verification.json").read_text())
assert verify["passed"] and complete["paper_updated"] is False
assert all(sha(ROOT / path) == digest for path, digest in complete["source_code_hashes"].items())
assert sha(ROOT / "config/integrated_review_policy.json") == complete["policy_sha256"]
kfm_check = json.loads((KFM / "verification.json").read_text())
assert kfm_check["input_hashes_verified"] and kfm_check["eligible_comparison_years"] == 52
for record in json.loads((KFM / "inputs.json").read_text()):
    path = Path(record["path"])
    assert path.is_file() and sha(path) == record["sha256"], path

annual = read("annual_values.csv")
joint = read("joint_model_results.csv")
pooled = read("pooled_results.csv")
decomp = read("input_change_decomposition.csv")
scope = pd.read_csv(ROOT / "config/integrated_comparison_scope.csv")
assert not annual.duplicated(["series_id", "year"]).any()
assert annual.series_id.nunique() == verify["candidate_histories"]
assert int(scope.exclude.sum()) == 31

# Model-based effects for one clearly defined estimand at a time.
targets = {
    ("post_trend", "recovery_slope"): "annual protected/reference change",
    ("flat_before_after", "average_observed_post_change"): "average observed before-after change",
}
effects = pooled[pooled.Source.eq("All programs")].copy()
effects["estimand"] = [targets.get((m, t), "") for m, t in zip(effects.model, effects.term)]
effects = effects[effects.estimand.ne("")].copy()
for column in ("estimate", "low", "high"):
    effects[column + "_percent"] = 100 * np.expm1(effects[column])
effects = effects.sort_values(["estimand", "Taxa", "Resp"])
saved(effects, "pooled_effects.csv")

# The first three differences use the same simple slope. They are data/coverage
# diagnostics, not a decomposition of biomass alone and not the joint model.
primary_slopes = joint[
    joint.analysis.eq("estimated_AR1")
    & joint.model.eq("post_trend")
    & joint.term.eq("recovery_slope")
][["series_id", "estimate", "se", "low", "high", "n_years"]].rename(
    columns={"estimate": "joint_slope", "se": "joint_se", "low": "joint_low", "high": "joint_high", "n_years": "joint_years"}
)
stages = decomp.merge(primary_slopes, on="series_id", how="left", validate="one_to_one")
stages["old_coverage_step"] = stages.old_shared_slope - stages.old_all_slope
stages["annual_input_step"] = stages.new_shared_slope - stages.old_shared_slope
stages["new_coverage_step"] = stages.new_all_slope - stages.new_shared_slope
stages["joint_model_step"] = stages.joint_slope - stages.new_all_slope
stages["annual_input_step_scope"] = "All annual-input changes on matched years; not biomass alone"
stages["joint_model_step_scope"] = "Same new years; conditional joint model versus unweighted diagnostic slope"
saved(stages, "series_change_steps.csv")

# Directional pathway screen on identical observed post-protection years within
# one program and reserve. No cross-taxon covariance is available for a new
# inferential pathway model, so no p-values are computed here.
links = [
    ("a", "P. interruptus", "Bio", "S. purpuratus", "Bio"),
    ("b", "P. interruptus", "Den", "S. purpuratus", "Den"),
    ("c", "S. pulcher", "Bio", "S. purpuratus", "Bio"),
    ("d", "S. pulcher", "Den", "S. purpuratus", "Den"),
    ("e", "S. purpuratus", "Den", "M. pyrifera", "Bio"),
    ("f", "M. franciscanus", "Den", "M. pyrifera", "Bio"),
]
post = annual[annual.response_supported & annual.t.ge(0) & annual.ln_ratio.notna()].copy()
pair_rows = []
annual_pair_rows = []
path_coverage = []
for letter, xt, xr, yt, yr in links:
    x = post[post.Taxa.eq(xt) & post.Resp.eq(xr)]
    y = post[post.Taxa.eq(yt) & post.Resp.eq(yr)]
    xy = x.merge(y, on=["MPA", "Source", "year"], suffixes=("_x", "_y"), validate="one_to_one")
    same_sites = (xy.protected_stations_x.eq(xy.protected_stations_y)
                  & xy.reference_stations_x.eq(xy.reference_stations_y)
                  & xy.shared_zones_x.eq(xy.shared_zones_y))
    path_coverage.append(dict(panel=letter, same_MPA_source_years=len(xy),
                              same_station_zone_years=int(same_sites.sum()),
                              different_footprint_years=int((~same_sites).sum())))
    xy = xy[same_sites]
    for (mpa, source), group in xy.groupby(["MPA", "Source"]):
        if len(group) < 3:
            continue
        years = group.year.to_numpy(float)
        for row in group.itertuples():
            annual_pair_rows.append(dict(panel=letter, MPA=mpa, Source=source,
                                         year=int(row.year), x=float(row.ln_ratio_x),
                                         y=float(row.ln_ratio_y)))
        pair_rows.append(dict(
            panel=letter, MPA=mpa, Source=source, x_taxon=xt, x_response=xr,
            y_taxon=yt, y_response=yr, years=len(group), first_year=int(years.min()),
            last_year=int(years.max()), x_slope=float(np.polyfit(years, group.ln_ratio_x, 1)[0]),
            y_slope=float(np.polyfit(years, group.ln_ratio_y, 1)[0]),
            shared_calendar_years=";".join(map(str, sorted(group.year.astype(int)))),
            protected_stations=";".join(sorted(set(";".join(group.protected_stations_x).split(";")))),
            reference_stations=";".join(sorted(set(";".join(group.reference_stations_x).split(";")))),
        ))
pathway = pd.DataFrame(pair_rows)
assert not pathway.duplicated(["panel", "MPA", "Source"]).any()
saved(pathway, "matched_year_pathway_screen.csv")
saved(pd.DataFrame(annual_pair_rows), "matched_year_pathway_annual.csv")
saved(pd.DataFrame(path_coverage), "matched_year_pathway_coverage.csv")
path_sum = []
for letter, g in pathway.groupby("panel"):
    path_sum.append(dict(panel=letter, pairs=len(g), MPAs=g.MPA.nunique(),
                         sources=";".join(sorted(g.Source.unique())),
                         positive_x=int((g.x_slope > 0).sum()),
                         negative_y=int((g.y_slope < 0).sum()),
                         opposite_direction=int(((g.x_slope * g.y_slope) < 0).sum()),
                         scope="descriptive; identical post-protection years; no cross-taxon uncertainty"))
saved(pd.DataFrame(path_sum), "matched_year_pathway_summary.csv")

# State-based heatwave check. Keep source-specific pairs so different survey
# scales are never averaged together. A pair needs observed baseline and period
# values on both sides; observed years and zero/undefined ratios are explicit.
windows = {"baseline": range(2010, 2014), "during": range(2014, 2017),
           "after": range(2017, 2020), "recent": range(2020, 2024)}
taxa_response = {"M. pyrifera": "Bio", "M. franciscanus": "Den",
                 "S. purpuratus": "Den", "P. interruptus": "Den", "S. pulcher": "Den"}
state = annual[annual.response_supported & annual.protected.notna() & annual.reference.notna()].copy()
state = state[state.apply(lambda row: taxa_response.get(row.Taxa) == row.Resp, axis=1)]
heat_rows = []
for (mpa, source, taxon), group in state.groupby(["MPA", "Source", "Taxa"]):
    w = {}
    for name, years in windows.items():
        z = group[group.year.isin(years)]
        w[name] = z
    base = w["baseline"]
    for period in ("during", "after", "recent"):
        z = w[period]
        if base.empty or z.empty:
            continue
        p0, r0 = base.protected.mean(), base.reference.mean()
        p1, r1 = z.protected.mean(), z.reference.mean()
        if min(p0, r0, p1, r1) <= 0:
            continue
        inside, outside = p1 / p0, r1 / r0
        footprint = pd.concat([base, z])[["protected_stations", "reference_stations", "shared_zones"]].drop_duplicates()
        heat_rows.append(dict(MPA=mpa, Source=source, Taxa=taxon, period=period,
                              baseline_years=len(base), period_years=len(z),
                              baseline_year_list=";".join(map(str, sorted(base.year.astype(int)))),
                              period_year_list=";".join(map(str, sorted(z.year.astype(int)))),
                              inside_ratio=inside, outside_ratio=outside,
                              log_inside_over_outside=float(np.log(inside / outside)),
                              constant_station_zone_footprint=len(footprint) == 1))
heat = pd.DataFrame(heat_rows)
assert not heat.duplicated(["MPA", "Source", "Taxa", "period"]).any()
saved(heat, "heatwave_state_pairs.csv")
heat_sum = []
for (source, taxon, period), group in heat.groupby(["Source", "Taxa", "period"]):
    values = group.log_inside_over_outside.to_numpy(float)
    p = float(wilcoxon(values).pvalue) if len(values) >= 5 and np.any(values != 0) else np.nan
    heat_sum.append(dict(Source=source, Taxa=taxon, period=period, pairs=len(group),
                         constant_footprint_pairs=int(group.constant_station_zone_footprint.sum()),
                         median_inside=float(group.inside_ratio.median()),
                         median_outside=float(group.outside_ratio.median()),
                         median_inside_over_outside=float(np.exp(np.median(values))),
                         paired_wilcoxon_p=p,
                         scope="source-specific state ratio; observed windows; conditional paired check"))
saved(pd.DataFrame(heat_sum), "heatwave_state_summary.csv")

# Plot six pathway panels without embedded title, subtitle, key or caption.
fig, axes = plt.subplots(2, 3, figsize=(11.5, 6.5))
plain_names = {"P. interruptus": "Lobster", "S. pulcher": "Sheephead",
               "S. purpuratus": "Purple urchin", "M. franciscanus": "Red urchin",
               "M. pyrifera": "Kelp"}
for axis, (letter, xt, xr, yt, yr) in zip(axes.flat, links):
    data = pathway[pathway.panel.eq(letter)]
    axis.axhline(0, color="#aaa", lw=.6)
    axis.axvline(0, color="#aaa", lw=.6)
    axis.scatter(data.x_slope, data.y_slope, s=22, alpha=.75, color="#246a73")
    axis.text(0, 1.02, letter, transform=axis.transAxes, va="bottom", fontweight="bold")
    axis.set_xlabel(plain_names[xt] + " " + ("biomass" if xr == "Bio" else "density") + " slope")
    axis.set_ylabel(plain_names[yt] + " " + ("biomass" if yr == "Bio" else "density") + " slope")
    axis.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(OUT / "matched_year_pathway_screen.png", dpi=180)
fig.savefig(OUT / "matched_year_pathway_screen.pdf")
plt.close(fig)

subprocess.run(["Rscript", "code/R/review_biomass_pathway_refits.R"], cwd=ROOT,
               env=dict(os.environ, OPENBLAS_NUM_THREADS="1"), check=True)
pathway_fits = pd.read_csv(OUT / "matched_year_pathway_fits.csv")
assert len(pathway_fits) == 6 * pathway.Source.nunique()
subprocess.run(["Rscript", "code/R/review_biomass_heatwave_refits.R"], cwd=ROOT,
               env=dict(os.environ, OPENBLAS_NUM_THREADS="1"), check=True)
heatwave_fits = pd.read_csv(OUT / "heatwave_period_refits.csv")
assert len(heatwave_fits) == 15 and heatwave_fits.status.eq("fitted").all()

source_hashes = {str(p.relative_to(ROOT)): sha(p) for p in [IN / "completed.json", IN / "verification.json",
    IN / "annual_values.csv", IN / "joint_model_results.csv", IN / "pooled_results.csv",
    IN / "input_change_decomposition.csv", ROOT / "config/integrated_comparison_scope.csv",
    KFM / "verification.json", KFM / "trends.csv"]}
(OUT / "input_hashes.json").write_text(json.dumps(source_hashes, indent=2))
checks = {"integrated_completed": True, "integrated_checks": verify["checks"],
          "retained_histories": annual.series_id.nunique(), "supported_comparison_years": int(annual.response_supported.sum()),
          "fitted_post_slopes": len(primary_slopes), "matched_pathway_pairs": len(pathway),
          "conditional_pathway_fits": int(pathway_fits.status.eq("fitted").sum()),
          "conditional_heatwave_taxa_fitted": int(heatwave_fits.Taxa.nunique()),
          "heatwave_state_pairs": len(heat), "KFM_sheephead_supported_years": 52,
          "paper_updated": False}
(OUT / "verification.json").write_text(json.dumps(checks, indent=2))

labels = {"M. pyrifera": "Kelp", "M. franciscanus": "Red urchin",
          "S. purpuratus": "Purple urchin", "P. interruptus": "Lobster",
          "S. pulcher": "Sheephead"}
post_effects = effects[effects.model.eq("post_trend") & effects.term.eq("recovery_slope")]
effect_table = "| Response | Comparisons | Annual change, 95% interval |\n|---|---:|---:|\n"
for row in post_effects.itertuples():
    effect_table += (f"| {labels[row.Taxa]} {'biomass' if row.Resp == 'Bio' else 'density'} | "
                     f"{row.comparisons} | {row.estimate_percent:+.1f}% "
                     f"({row.low_percent:+.1f}% to {row.high_percent:+.1f}%) |\n")

readme = f"""# Whole-dataset biomass update — review build

This build uses the corrected integrated source-to-model run completed {complete['completed']}. It is a review of the current working cohort. The manuscript and its frozen output files have not changed.

## Pooled post-protection trends

{effect_table}
These are conditional changes in the protected/reference ratio per observed calendar year. They are distinct from the flat-baseline before–after estimates in the linked effects table.

## Results to use for review

- [Pooled effects](pooled_effects.csv) separate post-protection slopes from average observed before–after changes. The [integrated review figures](../integrated_review_20261103/index.html) show conditional intervals.
- [Series change steps](series_change_steps.csv) compare the older annual inputs with the corrected inputs using one simple slope. Coverage changes, matched-year annual-input changes, and the change from that simple slope to the joint model are separate columns. The annual-input step includes all source, survey and biomass changes; it is not a biomass-only estimate. These diagnostic slopes do not replace the fitted joint model.
- [Matched-year pathway pairs](matched_year_pathway_screen.csv), [annual inputs](matched_year_pathway_annual.csv), [coverage losses](matched_year_pathway_coverage.csv), [conditional fits](matched_year_pathway_fits.csv), and [figure](matched_year_pathway_screen.png) use identical observed post-protection years, protected/reference stations and recorded zones within one source and MPA. The conditional fit uses calendar-gap outcome uncertainty, but does not estimate predictor error or cross-taxon covariance. It does not replace the paper's Table 3 without a model decision.
- [Heatwave state pairs](heatwave_state_pairs.csv) and [source-specific summary](heatwave_state_summary.csv) compare each source's protected and reference abundance in observed 2010–2013, 2014–2016, 2017–2019 and 2020–2023 windows. The pair table gives exact years and flags station/zone changes across windows. The [conditional period refit](heatwave_period_refits.csv) uses corrected annual ratios, source effects, an MPA intercept, calendar-gap AR1 and a linear recovery term. It is a review model, not an exact rerun of the old heatwave model or a causal test of protection.

## Scope and holds

The integrated run contains {annual.series_id.nunique()} response histories and {int(annual.response_supported.sum()):,} supported comparison-years. All {int(scope.exclude.sum())} researcher/ecological exclusions remain in force. Missing biomass never removes valid density. The separately checked KFM sheephead biomass has 52 supported comparison-years; it is not pooled into the full-history predator primary result. KFM lobster has no supported lengths. Scott's pending reference-site confirmation may change the KFM site roster.

The primary model carries conditional count, observed-size, shared-donor, shared-reference and calendar-process uncertainty. LTER group-mean variation, conversion-coefficient covariance, detection, habitat equivalence, predictor error and cross-taxon covariance in the pathway refit remain unresolved. The heatwave refit uses the corrected annual ratios but is a new conditional model; it has not replaced the paper's resilience suite. Final figure and manuscript concordance require a reviewed cohort and a decision about these downstream estimands.

## Rebuild

Run `python3 code/tools/run_whole_dataset_biomass_review.py` from the analysis repository to rebuild the integrated cohort, separate KFM fish review, and this summary in order. For an already completed integrated build, the driver accepts `--reuse-integrated`; use `--reuse-kfm` only when its source and review record have been checked. The report verifies the integrated build's completion and code/policy hashes. [Input hashes](input_hashes.json) and [checks](verification.json) record this review run.
"""
(OUT / "README.md").write_text(readme)
links = ["README.md", "pooled_effects.csv", "series_change_steps.csv", "matched_year_pathway_screen.csv", "matched_year_pathway_annual.csv", "matched_year_pathway_coverage.csv", "matched_year_pathway_effects.csv", "matched_year_pathway_fits.csv",
         "matched_year_pathway_summary.csv", "heatwave_state_pairs.csv", "heatwave_state_summary.csv", "heatwave_period_refits.csv", "input_hashes.json", "verification.json"]
body = "<h1>Whole-dataset biomass review</h1><p>Corrected cohort and models; manuscript unchanged.</p>"
body += f"<p>{checks['retained_histories']} retained histories; {checks['supported_comparison_years']:,} supported comparison-years; {checks['fitted_post_slopes']} fitted post-protection slopes.</p>"
body += '<p><a href="../integrated_review_20261103/index.html">Integrated model figures and diagnostics</a></p>'
body += '<img src="matched_year_pathway_screen.png" alt="Six matched-year pathway scatter panels">'
body += "<ul>" + "".join(f'<li><a href="{n}">{html.escape(n.replace("_", " "))}</a></li>' for n in links) + "</ul>"
(OUT / "index.html").write_text('<!doctype html><html><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Whole-dataset biomass review</title><style>body{font:16px/1.5 system-ui;max-width:1050px;margin:32px auto;padding:16px}img{max-width:100%;height:auto}a{color:#176371}</style>' + body + '</html>')
print(json.dumps(checks))
