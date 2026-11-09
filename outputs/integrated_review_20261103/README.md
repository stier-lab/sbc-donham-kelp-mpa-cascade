# Integrated pipeline review

Updated November 3, 2026. The paper and frozen outputs are unchanged.

The integrated build starts with registered contributor/provider files, applies source repairs and response-specific exclusions, rebuilds survey units and biomass, and fits before–after changes separately from post-protection trends. Main coverage uses documented stations available each year and shared program zones. Calendar years determine time; the designation year is after protection. Three matched before years permit before–after estimation; four is a sensitivity. Flat-baseline recovery is the main before–after model, with a pre-trend sensitivity.

## What the build contains

- 3,483 comparison-years with both roles and shared recorded zones; 2,939 have a supported response under the stated measurement rules.
- 175 histories reach the response/eligibility table; 162 have sufficient years for at least one fitted model.
- 31 candidate histories retain their recorded researcher/ecological exclusions. Another 0 retained candidates lack a usable comparison in this build; see the roster rather than interpreting them as ecological exclusions.
- 10,457 of 13,969 biomass visits have measurement support. Missing biomass never defines the density roster.
- 18 integration checks passed; 487 frozen files remain unchanged.

## How to read the results

[Coverage](coverage_summary.csv), [station assignments](station_assignments.csv), [comparison-year roster](comparison_year_roster.csv) and [model eligibility](model_eligibility.csv) show exactly which data enter each stage. [Visit biomass support](visit_biomass_support.csv), [donor attempts](donor_attempts.csv), [size pools](size_pool_inventory.csv) and [borrowing variance](borrowing_variance.csv) show how sizes are used.

[Joint models](joint_model_results.csv) preserve documented shared-reference process dependence as well as shared measurements. [Pooled results](pooled_results.csv) keep before–after changes, recovery slopes and post-only trends separate, and add MPA-level heterogeneity. Across-program summaries include source fixed effects and average the represented program means equally. [Coverage and sparse-size sensitivities](pooled_coverage_sensitivities.csv) use the same fitting procedure. These are conditional estimates, not total uncertainty.

[The paper-input comparison](paper_input_comparison.csv) compares annual responses directly with the frozen paper inputs. [The decomposition](input_change_decomposition.csv) uses one simple unweighted trend on old/all, old/matched, new/matched and new/all years to isolate data and coverage changes. Those diagnostic slopes are not the joint model estimates and do not reproduce the paper's former mixture of estimands.

## Interpretation limits

Of 60 main joint model groups, 36 estimate a variance component at zero and 13 reach the AR search boundary. These flags and the dependence sensitivities must accompany interpretation; they are not reasons to exclude a site. The [uncertainty ledger](uncertainty_limits.csv) distinguishes estimable variation from missing information.

The KFM sheephead restricted-period analysis remains a separately checked, 52-comparison-year secondary branch. KFM lobster contributes density; missing lengths have not been invented. Historical site-pool biomass is reproduced in the isolated source build but is not automatically adopted as supported biomass here.

## Rebuild

Run `python3 code/integrated_review/run.py` from the analysis repository. `--analysis-only` skips the inherited R reconstruction but still rebuilds and tests the registered source layer. The final completion record records code and policy hashes. Individual measurements remain in ignored private/staging directories; the source manifest retains access conditions.

Review these results and diagnostics before regenerating pathway/resilience analyses, final paper tables/figures, or manuscript text. Those paper outputs are not certified by this review.
