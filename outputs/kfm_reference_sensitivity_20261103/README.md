# KFM reference-site sensitivity

November 3, 2026. These tests answer Scott Gabara’s concerns about the current integrated review. They do not change the manuscript, primary roster or researcher exclusions.

## Findings

- **Harris Point:** Hare Rock is the protected station. Wyckoff Ledge is the only current reference for six responses; sheephead also uses Miracle Mile. Omitting Miracle Mile changes the sheephead post-protection ratio trend from −20.8% to −20.4% per year. Omitting Wyckoff leaves no reference for the six other responses. Sheephead retains 15 post-protection years, with no before baseline. On those same 15 years, the current-reference trend is −16.5% (95% interval −25.6 to −6.4), compared with −13.8% (−24.2 to −2.0) using Miracle Mile alone. The larger full-record change from −20.8% partly reflects lost years.
- **Santa Barbara Island:** Webster’s Arch affects sheephead density only; its benthic exclusion already applies. Removing it retains Arch Point and Cat Canyon as references and changes the local sheephead trend from +3.4% (−3.8 to +11.2) to −0.4% (−7.9 to +7.6) per year. Protected stations are SE Sea Lion Rookery, Graveyard Canyon and Southeast Reef, as available each year. No years are lost.
- **Rodes Reef:** no current measurement uses; omitting it adds no change. Gull KFM remains excluded under the current researcher-scope rules; Gull PISCO is separate and unchanged.
- **Combined omission:** removing Wyckoff, Miracle Mile and Webster’s Arch removes all seven Harris Point histories, loses 217 supported comparison-years, and reduces fitted histories from 162 to 155. Before–after eligible histories decline from 62 to 57. These are response histories and uses, not independent sites or unique calendar years.
- **Pooled post-protection trends:** kelp changes from +16.9% to +16.5% per year; purple-urchin density from −9.7% to −10.3%; lobster density from +7.1% to +10.4%; sheephead density from −3.1% to −1.6%. Predator trend intervals still include increases and decreases. Red-urchin density changes only modestly (−9.6% to −10.2%), although its interval endpoint crosses zero; do not turn that boundary change into a different biological conclusion.
- **Before–after lobster:** the flat-baseline average observed post contrast changes materially, from +229% (−40 to +1695) to +424% (+77 to +1451). The matched-cohort calculation attributes this to removing Harris Point, not to changing the remaining reference averages. Allowing the pre-protection trend to continue gives +158% (−97 to +19050) after omission. The evidence remains sensitive to the counterfactual model.
- **Before–after kelp:** the flat-baseline contrast changes from +280% (+56 to +827) to +292% (+63 to +843). The pre-trend alternative remains uncertain: about +10% before omission and −9% afterward, with very wide intervals in both cases.

## What the tests isolate

Four omission cases reuse the integrated estimator and joint uncertainty model. Dates, retained within-station observations, zone support, frozen ratio offsets, nonlinear size draws and within-station aggregation stay fixed. Equal weights among remaining stations are renormalized. A year without both roles cannot produce a paired estimate. No alternative station is introduced and no outcome determines the omission.

For Wyckoff and the combined case, a matched-main refit uses the original reference definition on exactly the alternative’s supported series/year keys. Full main → matched main measures the consequence of changed coverage. Matched main → omission measures the reference-definition change on that coverage, including reestimated shared process covariance. These differences add on the log-estimate scale; they are not separately identified causal mechanisms. Pooled weights and nuisance covariance can change after histories are removed.

Annual values and measurement draws for unaffected histories are unchanged. Their fitted coefficients can still move through reestimated covariance within a program/species/response group. This explains why the SBI estimate after combined omission (+3.2%, interval −19.1 to +31.8) differs from omitting Webster alone. It is not a second change to SBI’s observations. In the combined case, the KFM sheephead post-trend AR estimate reaches its upper bound of 0.95; fixed-AR alternatives also span zero at SBI. Treat this local uncertainty as dependence-sensitive.

Do not interpret the local trend as density change at the protected station alone: it is change relative to the reference. These are not habitat-validation tests. Miracle Mile is also on the southern side that Scott questioned, so its survival as a numerical comparison does not make it a preferred ecological replacement.

## Recommendation and remaining question

Keep these as explicit sensitivities while asking Scott whether he recommends retaining or omitting the Hare Rock versus Wyckoff/Miracle Mile comparison, and whether the Webster’s Arch habitat concern also warrants excluding its fish records. Do not choose by which result better supports a cascade. Rodes and the existing benthic exclusions need no further numerical repair. No new reference omission has been adopted here.

## Figure captions

**Local sheephead trends.** Points and lines show estimates and conditional 95% intervals for annual change in the protected/reference sheephead density ratio, using the unchanged zero-handling offsets. Panel A separates full-record, common-year and reference changes at Harris Point. Panel B shows the Webster-only and combined omissions at SBI. Current means the integrated review, not the unchanged manuscript. The intervals are for individual estimates, not tests of the difference between scenarios. No before–after result is plotted.

**Pooled trends.** Gray circles: current full cohort. Ochre squares: current references after removing the histories/years absent from the combined test. Green triangles: combined reference omission. Intervals are conditional 95% intervals. Overlapping ochre/green estimates show where the change is explained by loss of Harris Point rather than changes to remaining references. Positive values mean an increasing protected/reference ratio. No before–after result is plotted.

## Verification and reproducibility

All four omissions and two matched-main joint refits completed. Checks confirm exact retained measurement records, unchanged unaffected annual values and draws, fixed offsets, matched support keys, absence of omitted references, original ecological exclusion reasons and unchanged primary numerical files. All coefficient covariance matrices are symmetric and nonnegative within numerical tolerance; none of the 60 model groups per scenario reaches an upper variance bound. Depending on the scenario, 11–13 groups reach the AR search boundary and 36–38 have a zero variance component. These remain conditional fits with the same diagnostic qualifications as the current review; they are not estimates of total uncertainty.

The initial numerical run was followed by metadata-only ledger cleanup: existing exclusion reasons were preserved and dropped-row station metadata labeled as original candidate coverage. Neither changes model values or draws. The rerun script includes that cleanup from the start.

Run `python3 code/tools/run_kfm_reference_sensitivity.py`, then `check_kfm_reference_sensitivity.py`, `report_kfm_reference_sensitivity.py` and `plot_kfm_reference_sensitivity.py` from the analysis root. Individual source/draw files stay in the existing private sensitivity folders. `verification.json`, `checks.csv`, `diagnostics.csv` and the input records accompany this report. Exact scenarios and matching definitions are in `design.json`.
