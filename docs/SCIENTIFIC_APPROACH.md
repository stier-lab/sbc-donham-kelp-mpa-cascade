# Scientific approach and interpretation

Current as of 20 October 2026. This is analysis-facing guidance; manuscript files and private correspondence remain in the sibling manuscript repository. Dated audits record earlier states, not the current headline findings.

## Three questions

1. Do predator, grazer and kelp populations show directions compatible with a cascade? Use the primary population contrasts and elapsed-year diver trajectories. The synthesis has 145 source-specific histories from 21 MPAs, not 145 independent MPAs. Positive sheephead biomass and kelp contrasts and negative purple-urchin density contrasts match some predictions. Positive red-urchin biomass can reflect direct protection of a harvested grazer. Other primary population responses are inconclusive.
2. Do those population responses covary across MPAs? Use the six primary meta-regressions, exact permutation check, common-calendar programme comparisons and matched size analyses. The four-MPA lobster–purple-urchin biomass association is sensitive to the test and shared references. Neither a complete causal cascade nor absence of predator control is established.
3. Does kelp have a disturbance advantage inside MPAs, and is there evidence that designation improved that advantage? Keep the modern diver comparison separate from the historical satellite comparison. Modern baseline-relative biomass is higher inside MPAs, but is not necessarily regrowth after a loss. The same-pair historical comparison has intervals permitting poorer or better relative performance after designation. It does not establish either a designation benefit or inherent site resilience as the explanation.

## Quantities are not interchangeable

The primary synthesis combines before/after step contrasts, fitted model-time changes and mean log ratios over retained input years. Most diver trend models use survey indices; Landsat and one lobster history use elapsed years. The exponentiated pooled contrast is not a common biomass ratio or percentage increase. Table 2 and Figure 3 describe population evidence, not a uniform treatment effect.

The matched-design bridge instead compares arithmetic before/after changes and post-only levels on identical histories. For kelp, the annual-mean before/after estimate is 0.693 log units (95% interval 0.168–1.218; adjusted P = 0.046; 17 MPAs). Eleven MPAs contribute satellite observations only. The six-MPA diver-only interval includes zero. The descriptive calendar-AR1 GLS estimate is 0.072 (−0.683–0.827; 17 MPAs). GLS re-estimates period differences; it is not merely a new standard error for the arithmetic contrast. Both uncertainty and model specification limit a general designation-effect claim.

Sources: `tables/table_s_design_matched_summary.csv`, filtered by taxon M. pyrifera, response Bio, Region All and Matched before/after histories; distinguish Annual means from AR1. Construction is in module 27.

## Reef comparability

Module 45 provides a descriptive KFM-only diagnostic on fixed legacy reefs. Four core and three secondary non-target sessile taxa use three kelp-matched or four urchin-matched MPA calendars. Sparse records, heterogeneous starting differences and imprecise intervals do not validate equivalent reefs. Positive full-preperiod kelp point slopes are window-sensitive; they do not identify a stable untreated trajectory. These comparators cannot estimate the fraction of the kelp result caused by siting, and are not validated unaffected negative controls.

Repository-only as of 21 October 2026: Adrian removed this check from V27 and its supplement pending discussion with Bob. Retain its code, outputs and limitations here; do not automatically reinsert it into the paper. It does not correct the primary response. Full results and scope: `NONFISHED_COMPARISON_PLAN.md`, `NONFISHED_COMPARISON_RESULTS.md` and `tables/table_s_nonfished_*`.

## Satellite source decisions remain open

The corrected October delivery resolves the earlier mask discrepancies. All 24 MPA and 34 nonempty reference arrays reconstruct exactly from the corrected masks and cached NetCDF, including missingness. All 60 existing input rows match; the six previously removed reference slots remain empty. No primary values changed. Some reference splits were recovered against the delivered arrays, so reproduction does not establish a unique historical clustering solution.

Retain the current equal-available-reference procedure for this release; do not silently substitute preferred-reference selection, remove further observations or assign new historical footprints. Selection using a fixed pre-2012 correlation is not uniformly pre-designation for island MPAs. A preferred-reference sensitivity requires a documented source-specific selection rule. Inclusion and historical-footprint decisions require contributor confirmation. `LANDSAT_V2_VERIFICATION_2026-10-21.md` records the completed reconstruction and its limits; the 20 October audit is historical.

## Evidence hierarchy and completion

Main displays answer the three ecological questions. The essential supplement defines measurements, estimates and consequential sensitivities; repository outputs retain full diagnostics and exploratory checks. Descriptive and inconclusive results must not be promoted to causal identification. Agreement with predicted directions is not evidence that each trophic link operated.

Before final submission: resolve historical eligibility; document any approved reference-selection change and rerun all dependent outputs; confirm author roles; and run current manuscript, numerical, source and clean-reproduction checks. Completed satellite reconstruction does not resolve these separate decisions.
