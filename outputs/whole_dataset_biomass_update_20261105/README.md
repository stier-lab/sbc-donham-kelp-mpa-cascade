# Whole-dataset biomass update — review build

This build uses the corrected integrated source-to-model run completed 2026-11-05T19:58:21.917334+00:00. It is a review of the current working cohort. The manuscript and its frozen output files have not changed.

## Pooled post-protection trends

| Response | Comparisons | Annual change, 95% interval |
|---|---:|---:|
| Red urchin biomass | 11 | -13.2% (-26.0% to +1.8%) |
| Red urchin density | 19 | -9.6% (-18.3% to +0.0%) |
| Kelp biomass | 36 | +16.9% (+4.9% to +30.3%) |
| Lobster biomass | 5 | +6.2% (-19.7% to +40.6%) |
| Lobster density | 19 | +7.1% (-7.6% to +24.0%) |
| Sheephead biomass | 19 | -13.5% (-30.3% to +7.4%) |
| Sheephead density | 24 | -3.1% (-10.5% to +4.9%) |
| Purple urchin biomass | 10 | -6.9% (-18.7% to +6.6%) |
| Purple urchin density | 19 | -9.7% (-16.9% to -1.7%) |

These are conditional changes in the protected/reference ratio per observed calendar year. They are distinct from the flat-baseline before–after estimates in the linked effects table.

## Results to use for review

- [Pooled effects](pooled_effects.csv) separate post-protection slopes from average observed before–after changes. The [integrated review figures](../integrated_review_20261103/index.html) show conditional intervals.
- [Series change steps](series_change_steps.csv) compare the older annual inputs with the corrected inputs using one simple slope. Coverage changes, matched-year annual-input changes, and the change from that simple slope to the joint model are separate columns. The annual-input step includes all source, survey and biomass changes; it is not a biomass-only estimate. These diagnostic slopes do not replace the fitted joint model.
- [Matched-year pathway pairs](matched_year_pathway_screen.csv), [annual inputs](matched_year_pathway_annual.csv), [coverage losses](matched_year_pathway_coverage.csv), [conditional fits](matched_year_pathway_fits.csv), and [figure](matched_year_pathway_screen.png) use identical observed post-protection years, protected/reference stations and recorded zones within one source and MPA. The conditional fit uses calendar-gap outcome uncertainty, but does not estimate predictor error or cross-taxon covariance. It does not replace the paper's Table 3 without a model decision.
- [Heatwave state pairs](heatwave_state_pairs.csv) and [source-specific summary](heatwave_state_summary.csv) compare each source's protected and reference abundance in observed 2010–2013, 2014–2016, 2017–2019 and 2020–2023 windows. The pair table gives exact years and flags station/zone changes across windows. The [conditional period refit](heatwave_period_refits.csv) uses corrected annual ratios, source effects, an MPA intercept, calendar-gap AR1 and a linear recovery term. It is a review model, not an exact rerun of the old heatwave model or a causal test of protection.

## Scope and holds

The integrated run contains 175 response histories and 2,939 supported comparison-years. All 31 researcher/ecological exclusions remain in force. Missing biomass never removes valid density. The separately checked KFM sheephead biomass has 52 supported comparison-years; it is not pooled into the full-history predator primary result. KFM lobster has no supported lengths. Scott's pending reference-site confirmation may change the KFM site roster.

The primary model carries conditional count, observed-size, shared-donor, shared-reference and calendar-process uncertainty. LTER group-mean variation, conversion-coefficient covariance, detection, habitat equivalence, predictor error and cross-taxon covariance in the pathway refit remain unresolved. The heatwave refit uses the corrected annual ratios but is a new conditional model; it has not replaced the paper's resilience suite. Final figure and manuscript concordance require a reviewed cohort and a decision about these downstream estimands.

## Rebuild

Run `python3 code/tools/run_whole_dataset_biomass_review.py` from the analysis repository to rebuild the integrated cohort, separate KFM fish review, and this summary in order. For an already completed integrated build, the driver accepts `--reuse-integrated`; use `--reuse-kfm` only when its source and review record have been checked. The report verifies the integrated build's completion and code/policy hashes. [Input hashes](input_hashes.json) and [checks](verification.json) record this review run.
