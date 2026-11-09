# Lobster and purple-urchin size uncertainty

November 1 update: this sensitivity was rerun on the 143-history satellite-screened cohort. All nine directions persist, but the previous statement that all adjusted decisions persist is superseded. Purple-urchin biomass meets the threshold in 30.6% of individual-weight and 50.3% of survey-group draws; red-urchin biomass and lobster density in 39.0% and 58.6%. Independent second-seed detection fractions differ by at most 0.062 and lie within the conservative four-SE Monte Carlo check. The dated account below describes the preceding cohort; current tables and [satellite policy](SATELLITE_BASELINE_POLICY_2026-11-01.md) control current interpretation.

The pooled population results withstand the tested changes in lobster and purple-urchin sizes. Evidence for the four-MPA lobster–urchin relationship is less stable: its slope remains negative, but its nominal interval often includes zero. This supports retaining the paper's cautious interpretation of trophic links.

## What was tested

The analysis retains all 145 primary histories, their fitted years and model forms, using the current calendar clock and surveyed endpoints. It changes the size distributions used by 13 biomass histories: three PISCO lobster histories and ten purple-urchin histories. The two LTER lobster histories use directly measured lengths and retain their corrected primary values.

The shared source-record ledger contains 27,211 original frequency records feeding 414 positive-biomass pools: 353 paired-year pools and 61 additional pools that contribute to the original zero adjustment. Every occurrence of a source record receives the same perturbation wherever it is reused. South Point's partially overlapping KFM/PISCO pools retain their different diameter filters. Campus Point and Naples retain their shared reference measurements.

Two sets of 1,000 positive-weight perturbations examine different assumptions:

- **Individual measurements:** a frequency record representing `n` animals receives a Gamma(n, 1) total weight. Mean mass is calculated from the weighted masses of individuals, preserving the nonlinear length–mass relationship.
- **Recorded survey groups:** observations sharing program, method, site, year, recorded date, location, zone and transect receive one exponential weight. Animal frequencies still determine the original estimator's weighting. Missing identifiers remain missing; these groups are not claimed to be independently verified dives.

These are conditional reweighting diagnostics, motivated by the [Bayesian bootstrap](https://doi.org/10.1214/aos/1176345338). Neither scheme is inherently conservative. Positive weights retain every observed size and prevent an empty pool; they cannot invent unobserved sizes.

Each annual biomass is multiplied by the perturbed mean individual mass divided by its original mean mass. This anchors the calculation to the released biomass, including its original simulation approximation, while holding count surveys and coefficients fixed. The adaptive zero adjustment is recalculated over the complete original annual domain before selecting fitted years. Gull Island's purple-urchin adjustment depends on an unpaired reference observation in 2002 (4.859944 g/m²), absent from the paired-only raw export. Its sizes are included in the perturbations.

The retained contrasts and their conditional standard errors are recalculated, followed by the original joint meta-analysis and four-MPA lobster–purple-urchin regression. The nine-test adjustment is retained. Failure stops the run; failed fits are not silently removed.

Separate stress tests halve or double the assumed mean mass in each of the 17 sparse lobster pools, one at a time (34 scenarios). Four additional scenarios change protected and reference masses in opposite directions, either for the six singleton pools or for all 17 sparse pools. These factors are imposed assumptions, not estimated confidence limits. Only retained paired-year lobster pools receive these stress factors.

## Findings

All nine pooled directions and the same four adjusted statistical classifications persist across the 2,000 perturbations and all 38 stress scenarios. The second seed, with another 250 perturbations per scheme, gives the same pooled conclusions.

| Result | Primary estimate | Individual-weight central 95% range | Survey-group-weight central 95% range |
|---|---:|---:|---:|
| Pooled lobster biomass contrast | 0.346 | 0.295 to 0.389 | 0.302 to 0.397 |
| Pooled purple-urchin biomass contrast | -0.739 | -0.761 to -0.718 | -0.792 to -0.701 |
| Four-MPA lobster–purple-urchin slope | -0.555 | -0.616 to -0.508 | -0.631 to -0.490 |

These ranges describe variation in **point estimates under the specified perturbations**. They are not replacement confidence intervals, and mixed mean/trend contrasts are not common fold changes.

The negative purple-urchin response persists, while the pooled lobster biomass estimate remains too uncertain to establish a clear overall increase. Under the deliberately imposed sparse-pool stress scenarios, the pooled lobster estimate ranges from 0.140 to 0.567; its adjusted P remains above 0.05 (range 0.077–0.667).

The four-MPA slope is negative in every perturbation. However, its nominal 95% interval includes zero in **38.9%** of individual-weight perturbations and **27.3%** of survey-group perturbations. The second seed gives 37.6% and 30.0%, respectively. Thus the sign is stable under these checks, while the apparent precision is sensitive to the observed sizes. These frequencies are not probabilities that a trophic relationship exists or is absent. The exact-permutation test was not rerun for each perturbation.

Across seeds, the largest pooled median difference is 0.0017 and the largest central-range endpoint difference is 0.0047. No pooled direction or adjusted classification changes.

## What remains unresolved

- Six retained lobster pools contain one measured animal; 142 of the 353 retained positive pools contain only one recorded survey group. Their lack of within-pool variation under the corresponding scheme does not establish precise biomass estimates.
- Shared weights preserve dependence in size perturbations. They do not repair all covariance assumptions inside each original fitted model, shared count surveys or the treatment of the pathway predictor as fixed.
- Count sampling, visibility/detection, size representativeness, measurement error, allometric coefficient uncertainty and other taxa's conversion uncertainty are held fixed. The LTER direct-size histories receive no new sampling-uncertainty analysis.
- Figure 4's mixed model, calendar-AR1 fits and seasonal sensitivities were not refit in this targeted check. Existing temporal checks address different uncertainty questions.

The pooled interpretation can remain. The trophic interpretation should remain cautious, and these results should not be used to narrow the published intervals or claim complete uncertainty propagation. A short SI description can report this targeted check; the detailed ledgers and diagnostics belong in the repository. No primary numerical result or manuscript/SI prose was changed by this analysis.

## Reproduction and verification

From the analysis repository, with the source-data repository available through `DONHAM_DATA_REPO` or its default sibling path:

```sh
Rscript code/R/run_size_uncertainty.R
```

The runner refreshes module 55's inventory, reconstructs the full annual domain with existing source caches, runs module 57, checks independent reconstruction, and repeats the diagnostic with a second seed. Source cache replacement is prohibited. Raw membership and draw-level records remain ignored under `local_archive/primary_size_uncertainty_2026-10-30/`.

Compact outputs are `tables/table_s_size_uncertainty_{pools,contrast_summary,pooled_summary,pathway_summary,covariance,contrasts,pooled,pathway,inputs,status,seed_comparison}.csv`. The contrasts, pooled and pathway tables contain the baseline and explicit stress cases; draw-level results remain local. `inputs` pins the full-domain manifest, raw-size file, source membership, code and protected primary files.

Executed checks reproduce all 145 primary estimates/SEs, all nine pooled estimates/adjusted P values and the four-MPA baseline slope/P. Independent calculations verify frequency weighting, partial overlap, shared survey weights, singleton limits, every first-draw retained contrast, and the unpaired-year zero-adjustment dependency. Primary input, result-table and temporal-cache hashes are unchanged. Source cache hashes are unchanged. The inventory check and independent second-seed check pass. All conditional size fits completed; rejected nonlinear initialization attempts elsewhere in the pipeline do not replace the converged retained fits.

The earlier paired-only exploratory run is archived locally and is superseded by the complete-domain results above.
