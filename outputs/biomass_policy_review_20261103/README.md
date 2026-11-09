# Biomass support and uncertainty: current review

November 3, 2026. These are review calculations. They have not replaced manuscript results.

## Coverage

The [comparison table](comparison_support.csv) lists each protected area, program and response with its candidate, supported and unsupported biomass years, and years that involve borrowed sizes. The [year table](comparison_year_support.csv) shows the underlying visits and support state. Of 1,990 candidate biomass comparison-years across 94 histories, 1,449 support a biomass ratio. Counts and density remain available where their surveys support them; a missing size is not a site exclusion. PISCO has 446 supported of 713 candidate biomass years; 20 of those supported years involve approved borrowed sizes.

## What borrowing changes

A [direct-size-only refit](../integrated_review_20261103/sensitivities/direct_sizes_only/coverage_changes.csv) withholds the 137 calibrated borrowed-size visits and keeps the source roster, counts, exclusions, sites, zones and model code fixed. Twenty previously supported PISCO biomass comparison-years become unavailable: eight red urchin, seven lobster and five purple urchin. Those years occur in six histories. The same 2,919 other supported annual estimates and measurement draws are numerically identical to the integrated review. All density pooled results are unchanged to numerical precision.

The pooled all-program red-urchin biomass recovery slope changes from -13.19% per year (conditional 95% interval -25.97% to +1.79%; 11 comparisons) to -18.38% (-31.53% to -2.70%; nine comparisons). The purple-urchin slope changes from -6.92% (-18.74% to +6.62%; ten comparisons) to -7.43% (-16.40% to +2.50%; nine comparisons). Lobster biomass remains uncertain. **These changes reflect lost years and comparison membership; the retained yearly biomass values did not change.** They are not evidence that borrowing creates or hides a protection effect. The [pooled comparison](pooled_sensitivity_results.csv) lists all responses.

The earlier six-history held-out donor check is [available separately](../biomass_next_steps_20261102/borrowing/heldout_validation.csv). At zero extra timing tolerance, 420 eligible other-site predictions had 93.6% interval coverage and median absolute log prediction error 0.313. The eight eligible other-visit predictions had 100% coverage and median absolute log error 0.474. These repeated calibration records are not 428 independent sites; coverage and point error both matter. The integrated donor-attempt and calibration ledgers remain the authority for which borrowed visits were actually used.

## LTER group means

The source records give group mean urchin diameter, not individual diameters. A mean-preserving lognormal illustration with an assumed within-group diameter coefficient of variation (CV) of 10%, 20% or 30% changes red-urchin mean mass by factors of 1.028, 1.116 or 1.272; purple-urchin mean mass by 1.027, 1.111 or 1.260. [The formula table](lter_group_mean_scenarios.csv) records these assumptions. They are **not estimates of actual LTER variation**.

I refit the integrated model under the 20% and 30% cases, applying each taxon's same factor to all positive LTER group-mean biomass values while preserving observed zeros and coverage. Pooled red-urchin biomass changes from -13.19% per year to -13.26% and -13.33%; purple changes from -6.92% to -6.91% and -6.89%. The near cancellation occurs because the same factor applies to protected and reference observations. These fits do not test unequal size variation between roles, years or size classes, and they do not recover individual-size uncertainty. Do not describe these factors as confidence bounds.

## Conversion equations

The [coefficient-status table](conversion_uncertainty_status.csv) records the formulas used and the absence of coefficient covariance in the current input files. The point equations are retained. A probabilistic coefficient resampling analysis needs the underlying calibration observations or a published coefficient variance-covariance matrix. Scaling each coefficient independently by an arbitrary percentage would not provide a defensible confidence interval. Kelp's month-specific coefficients are observed calibration inputs, but their sampling covariance is likewise unavailable here.

## Emily's archived size pool

The [common-year archive comparison](../emily_biomass_bridge_20261103/README.md) now reads Emily's archived combined size file and physical-site crosswalk, with her listed site exclusions and 25 mm urchin threshold. On the exact current PISCO count roster, 30 comparison-years across six urchin histories support both size approaches. Descriptive red-urchin slopes average -21.65% using current visit-supported sizes and -22.66% with Emily's broader pool; purple slopes average -4.78% and -11.40%. These are selected matched histories and unweighted descriptive averages, not the pooled paper model. The matched joint-model refits separate size-pool point estimates from their induced measurement covariance on those 30 years. Purple-urchin biomass moves by -4.92 percentage points per year when the size pool changes with current covariance fixed, versus +0.42 points when only the covariance changes with current biomass fixed. This is a limited common-cohort attribution, not a replay of Emily's 1,000 draws or her downstream pBACIPS model.

## Decision supported by this review

Retain the documented same-visit preference and narrow, calibrated reference borrowing for the review cohort. Keep density where biomass is unsupported. Show borrowing dependence and the smaller direct-only cohort alongside any eventual biomass result. The matched model-level split is complete for the 30 supported PISCO urchin years; those selected years cannot stand in for the full biomass cohort. Emily's original R-model uncertainty has not been replicated. LTER within-group and allometric-coefficient uncertainty remain unmeasured by the source files. Reviewer decisions should turn on donor comparability and source evidence, not whether an interval crosses zero.

Rebuild the tables with `python3 code/tools/report_biomass_support.py` and `python3 code/tools/compare_emily_size_pool.py`. The direct-only and LTER scenarios use `code/tools/run_direct_size_only.py` and `code/tools/run_lter_group_mean_scenario.py` with `.2` or `.3`. [Validation](verification.json) checks fixed rosters, matched annual values and support.
