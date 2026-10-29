# Primary Processing Labels

`AnalysisType` in the primary exports is a retained processing label, not a
complete before/after eligibility audit. Do not infer that all CI series failed
a pretrend screen or lacked before-period observations.

Run `python3 code/tools/check_processing_labels.py` after the primary and matched
coverage exports. It joins the exact retained identities and writes
`tables/table_s_processing_routes.csv` without changing any fit, estimate,
variance or primary membership.

The current 145 time series carry 108 CI and 37 pBACIPS labels. Of the CI series,
84 have fewer than three paired before-years and 24 have at least three. All
89 PISCO series use the direct CI route, including 22 with at least three
before-years. The stored CI rows do not provide failed-screen P values.

Campus Point LTER lobster biomass is directly constructed as a calendar-time
linear CI trend in module 08, but its helper exports the pBACIPS label. The
coverage table records one before-labelled year. The 37-fit display therefore
includes this exception and must not be called 37 independently verified
before/after-eligible fits. The matched comparison separately contains 60 time
series with at least three paired years on each side of designation.

The V30 manuscript and essential supplement report construction design:
109 CI and 36 pBACIPS contrasts. Table S2 joins this route metadata to the
model-selection table and assigns the Campus Point LTER lobster-biomass trend
to CI for display (41 CI linear contrasts, 27 pBACIPS linear contrasts).
The CI coverage summary consequently has 85 series with fewer than three
paired before-years and 24 with at least three. This presentation change does
not alter stored fit labels, estimates, variances or primary membership.

These facts clarify the analysis actually performed. They do not approve the
inherited routes, reconstruct every pretrend decision, or substitute a new
before/after refit for the primary synthesis. Any future route change requires
new primary outputs and their dependent uncertainty, display and manuscript
checks.
