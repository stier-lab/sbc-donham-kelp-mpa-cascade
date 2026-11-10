# Code index

Generated 2026-11-09. Paths are stable because dated reports cite scripts by path; do not move a script without updating every reference.

## Current entry points

| Script | Purpose |
|---|---|
| `code/tools/run_whole_dataset_biomass_review.py` | Rebuild registered sources, comparisons, supported biomass, joint uncertainty, pathway and heatwave checks (whole-dataset review). |
| `code/integrated_review/run.py` | Integrated source-to-recovery review; policy in config/integrated_review_policy.json. |
| `code/integrated_review/check.py` | Integrated review invariants. |
| `code/tools/report_dataset_inventory.py` | Regenerate the current decision inventory and viewer (no model reruns). |
| `code/tools/check_current_decisions.py` | Verify current decision tables against their sources. |
| `code/tools/check_documentation_consolidation.py` | Check guide links, archive manifests and status assertions. |
| `code/tools/profile_joint_model_data.py` | Data profile and prior predictive check for the joint hierarchical model design. |
| `code/R/run_all.R` | Frozen manuscript pipeline (current paper results). |
| `code/R/run_figures_only.R` | Regenerate paper figures from the cached snapshot. |
| `code/R/run_resilience.R` | Resilience suite; membership in code/R/resilience_modules.R. |

## Packages and pipelines

- `code/R/` — frozen manuscript pipeline modules (00–58). Module map: `CLAUDE.md`.
- `code/integrated_review/` — integrated review package (`run.py`, `biomass.py`, `coverage.py`, `estimate.py`, `sensitivities.py`, `report.py`, `check.py`).
- `code/replications/` — external-method replications (Malakhoff 2021).
- `code/tests/` — unit and regression tests.
- `code/tools/archive/` — scripts no document or other script references (moved 2026-11-09).

## Review and audit scripts, by the report they produce or check

Status of each report is in `outputs/README.md`. These scripts reproduce dated evidence; they are not the current pipeline unless listed above.

### `outputs/all_program_station_refits_20261102/`

- `code/tools/refit_all_program_stations.py`
- `code/tools/report_all_program_stations.py` — Display the largest omission shifts per program among newly examined responses.
- `code/tools/stage_all_program_station_refits.py`
- `code/tools/verify_all_program_omissions.R`
- `code/tools/verify_all_program_stations.R`

### `outputs/anacapa_comparability_20261101/`

- `code/tools/audit_anacapa_comparability.py` — Source-record audit; isolated diagnostics, never alters primary inputs.
- `code/tools/plot_anacapa_comparability.py` — Date-matched year illustration: zones retained as essential identifiers.

### `outputs/archive/`

- `code/tools/resolve_biomass_support.py` — Review-only size-support sensitivities. Never writes manuscript inputs.

### `outputs/biomass_exclusion_reconciliation_20261102/`

- `code/tools/verify_kfm_excluded_donor_refits.R`

### `outputs/biomass_joint_checks_20261102/`

- `code/tools/export_biomass_joint_sources.R` — Read only: extract the previously reconstructed processor objects to private audit files.
- `code/tools/report_biomass_joint.py` — Companion artwork contains identifiers and axes only; explanation is in the report.

### `outputs/biomass_next_steps_20261102/`

- `code/tools/biomass_observation_model.py` — Conditional joint size/count observation model and temporal-process diagnostic.
- `code/tools/rebuild_borrowing_support.py` — Read-only input, diagnostic-only borrowing rebuild. Requires numpy/pandas/scipy.
- `code/tools/rebuild_density_baseline.py` — Read-only six-history PISCO density rebuild. Requires Python 3, pandas, numpy.
- `code/tools/run_biomass_next_steps.py` — Rebuild scoped review evidence; canonical analysis and manuscript are untouched.
- `code/tools/test_density_baseline.py` — Independent csv/stdlib source reconstruction plus numerical edge tests.

### `outputs/biomass_policy_review_20261103/`

- `code/tools/run_direct_size_only.py` — Review-only direct-size sensitivity, preserving existing density and survey roster.
- `code/tools/run_lter_group_mean_scenario.py` — Illustrative group-size CV sensitivity; assumed CV is not observed LTER data.

### `outputs/biomass_protocols_20261102/`

- `code/tools/check_kfm_date_timing.py` — Only site-years actually supplying retained donor pools are relevant to the review cohort.
- `code/tools/check_kfm_size_site_labels.py`
- `code/tools/reconcile_kfm_urchin_dates.py` — Recover dates only from exact site/year/species size-frequency agreement.
- `code/tools/report_biomass_protocols.py`
- `code/tools/verify_sbi_size_reconciliation.R`

### `outputs/biomass_resolution_20261102/`

- `code/tools/report_biomass_resolution.py` — Distinguish loss of years from changing donors on those same years.
- `code/tools/verify_biomass_resolution.R`
- `code/tools/verify_missing_stipe_guard.R` — Exercise exact source expressions against the real reconstructed station table.

### `outputs/biomass_support_20261102/`

- `code/tools/audit_direct_biomass_support.py`
- `code/tools/audit_expanded_size_pools.R` — Expanded inventory; retain raw size records and primary calculations.
- `code/tools/report_biomass_support.py` — Current review biomass support by physical comparison/year; no paper edits.
- `code/tools/verify_biomass_missing_stipes.R`
- `code/tools/write_biomass_support_report.py` — What is supported

### `outputs/change_attribution_20261101/`

- `code/tools/attribute_analysis_changes.R`
- `code/tools/check_change_attribution.py`
- `code/tools/report_change_attribution.py`

### `outputs/coverage_decisions_20261101/`

- `code/tools/audit_coverage_decisions.py` — Coverage diagnostics only; preserve method identities and outcome-independent choices.
- `code/tools/report_coverage_decisions.py` — Decisions recorded

### `outputs/cross_site_seasons_20261101/`

- `code/tools/audit_cross_site_seasons.py` — Isolated, source-based seasonal and timing coverage audit. See output PLAN.md.
- `code/tools/refit_cross_site_seasons.R`
- `code/tools/report_cross_site_seasons.py` — Fixed physical station pairs across every displayed timing-rule contrast.
- `code/tools/verify_cross_site_seasons.py`

### `outputs/dataset_inventory_20261101/`

- `code/tools/audit_all_dataset_sources.py` — Read-only source inventory. No abundance estimation, exclusions or model changes.
- `code/tools/audit_dataset_measurement_design.py` — Read-only measurement/design screen; raw scope, never automatic exclusions.
- `code/tools/extend_dataset_audit.py` — Resolve cross-program mapping labels and cover all24 delivered satellite polygons.
- `code/tools/inventory_measurement_review.py` — Interpret source screens without changing scientific eligibility.
- `code/tools/inventory_site_resolution.py` — Apply reviewed site hierarchy to inventory only; never change analytic membership.
- `code/tools/verify_consolidated_inventory.py` — Check consolidated values, missingness, keyed joins and viewer payload against dated evidence.
- `code/tools/verify_consolidated_models.R` — Independent R refit of every displayed diagnostic coefficient and uncertainty interval.
- `code/tools/verify_dataset_inventory.py` — Independent reconciliation of coverage tables, keys and unchanged source artifacts.

### `outputs/depth_reference_policy_20261101/`

- `code/tools/report_depth_reference_policy.py` — Distinct pair summaries, not independent replicates or uncertainty intervals.
- `code/tools/verify_depth_reference_policy.py` — Raw depths independently re-read. Dates/transects, not taxon records, define units.

### `outputs/emily_biomass_bridge_20261103/`

- `code/tools/compare_emily_size_pool.py` — Review-only common-cohort reconstruction of the archived PISCO urchin size pool.

### `outputs/fish_lter_reconciliation_20261102/`

- `code/tools/decompose_lter_reconciliation.py` — Source quadrats: unavailable counts remain unavailable; no biological absence inferred.
- `code/tools/reconcile_lter_annual_biomass.py` — Compare current raw kelp reconstruction to pinned SBC50.21; no promotion.
- `code/tools/reconcile_scott_fish.py` — Read-only source reconciliation; current count roster, no biomass promotion.

### `outputs/full_followup_recovery_20261101/`

- `code/tools/check_full_followup_recovery.py`
- `code/tools/check_recovery_influence.R`
- `code/tools/refit_full_followup_recovery.R`
- `code/tools/report_full_followup_recovery.py`

### `outputs/joint_borrowed_biomass_20261102/`

- `code/tools/joint_borrowed_biomass.py` — Joint working Gaussian likelihood for raw size means and annual count contributions.
- `code/tools/run_joint_biomass_review.py` — Rebuild the isolated joint-model review, checks and consolidated guide.

### `outputs/kfm_git_origin_20261102/`

- `code/tools/audit_kfm_site_identity_history.py`

### `outputs/kfm_missing_size_recovery_20261103/`

- `code/tools/audit_kfm_missing_size_recovery.py` — Trace seven missing visits and quantify review-only visit retention alternatives.
- `code/tools/check_kfm_missing_size_recovery.py`

### `outputs/kfm_pipeline_guard_20261102/`

- `code/tools/audit_older_kfm_identity.py`
- `code/tools/build_kfm_identity_bundle.py`
- `code/tools/check_kfm_pipeline_guard.R` — A diagnostic subset tests repaired joins; it is NOT adopted as eligibility.

### `outputs/kfm_reference_sensitivity_20261103/`

- `code/tools/check_kfm_reference_sensitivity.py` — Check reference omissions, matched-year decomposition and unchanged primary files.
- `code/tools/plot_kfm_reference_sensitivity.py`
- `code/tools/report_kfm_reference_sensitivity.py`
- `code/tools/run_kfm_reference_sensitivity.py` — KFM reference sensitivity: same pipeline, explicit omissions, no primary promotion.

### `outputs/kfm_remaining_sources_20261102/`

- `code/tools/check_remaining_kfm_sources.py`
- `code/tools/verify_hare_1987_source.py`
- `code/tools/verify_remaining_kfm_effects.R`

### `outputs/kfm_scott_reconciliation_20261102/`

- `code/tools/reconcile_scott_kfm_sizes.py` — Reconcile original Scott size frequencies and refit isolated review histories.
- `code/tools/verify_scott_kfm_refits.R`

### `outputs/kfm_sheephead_biomass_20261103/`

- `code/tools/kfm_sheephead_biomass_review.py` — Restricted KFM review: nonlinear joint count/size propagation; no primary edits.
- `code/tools/report_kfm_sheephead_biomass.py` — Figures, review report and cross-scenario summary for restricted KFM biomass.
- `code/tools/run_kfm_sheephead_biomass_review.py` — Rebuild the restricted KFM review, its sensitivities, report and guide tables.

### `outputs/kfm_source_comparison_20261102/`

- `code/tools/compare_kfm_original_measurements.py` — Independent, count-preserving comparison of raw and verified KFM urchin sizes.

### `outputs/lter_integration_pisco_check_20261102/`

- `code/tools/check_cached_pisco_sizes.py`
- `code/tools/report_lter_integration.py`

### `outputs/lter_original_sizes_20261102/`

- `code/tools/check_lter_original_sizes.py` — Site names independently decoded in MBON's provider-integrated count table.

### `outputs/main_figure_preview_20261106/`

- `code/tools/preview_main_figures.py` — Review-only previews of five main figure roles from the integrated cohort.

### `outputs/model_eligibility_20261101/`

- `code/tools/audit_baseline_model_eligibility.py`
- `code/tools/check_baseline_eligibility.R`
- `code/tools/check_common_estimand.py`
- `code/tools/plot_model_eligibility.py`
- `code/tools/pool_model_eligibility.R` — Reuse the same pooling specification; input effects are unchanged primary contrasts.
- `code/tools/report_model_eligibility.py`

### `outputs/nonlinear_biomass_calibration_20261102/`

- `code/tools/nonlinear_biomass_calibration.py` — Paired conditional nonlinear calibration of the joint working likelihood.
- `code/tools/run_nonlinear_biomass_review.py` — Rebuild conditional nonlinear calibration and documentation; no paper changes.

### `outputs/parallel_resolution_20261102/`

- `code/tools/audit_biomass_joint.py` — Conditional joint measurement resampling; immutable primary inputs.
- `code/tools/audit_biomass_matching.py` — Source-record count/size compatibility checks; no new cutoffs or primary edits.
- `code/tools/audit_depth_reference_policy.py` — Outcome-blind coverage/depth audit; no outcome model fitted here.
- `code/tools/audit_kfm_2015_aliases.py` — Capture every retained history using the suspect delivered labels. Post2015 identity
- `code/tools/audit_other_program_footprints.py` — Match known retained site/month membership, not every station in a broad MPA map.
- `code/tools/audit_pisco_size_zones.py` — Trace 23 size-zone flags to raw records; isolated matched-year biomass sensitivities.
- `code/tools/audit_station_coverage_refits.py`
- `code/tools/biomass_source_guard.py` — Admission guard for review biomass draws; never changes density or raw data.
- `code/tools/build_dataset_inventory.py` — Build linked, review-only dataset/site/comparison inventory from raw scan and evidence.
- `code/tools/check_pisco_timing_donors.py` — Review-only date pairing and held-out size-donor validation. No primary changes.
- `code/tools/collect_parallel_resolution.py` — Collect isolated agent evidence without adopting recommendations as policy.
- `code/tools/consolidate_current_decisions.py` — Current decision views; overlays documentation, never scientific estimates.
- `code/tools/inventory_latest_audits.py` — Consolidate dated review evidence without promoting it to the manuscript.
- `code/tools/prepare_all_program_stations.py`
- `code/tools/propagate_pisco_borrowing.py` — Review pipeline: all-month density and empirical predictive biomass draws.
- `code/tools/refit_kfm_excluded_donors.py` — Isolated correction of proven excluded-donor leakage, using original <=2011 sizes.
- `code/tools/refit_pisco_all_months.py` — Review-only all-month PISCO swath refits; six flagged animal comparisons.
- `code/tools/refit_sbi_size_reconciliation.py`
- `code/tools/report_kfm_exclusion_reconciliation.py` — Why this matters
- `code/tools/report_pisco_borrowing.py`
- `code/tools/run_pisco_biomass_review.py` — Single entry point for the isolated, non-production PISCO biomass review.
- `code/tools/verify_biomass_joint.R` — Verify independent arithmetic for KFM expected mass using original plant records.
- `code/tools/verify_pisco_borrowing.R`

### `outputs/pelican_pisco_recovery_20261102/`

- `code/tools/review_pelican_pisco_recovery.py` — Strict sums/means retain an explicit unknown measurement throughout this site-year.

### `outputs/pisco_all_months_20261102/`

- `code/tools/report_pisco_all_months.py`
- `code/tools/verify_pisco_all_months.R`

### `outputs/pisco_fish_effort_20261101/`

- `code/tools/audit_pisco_fish_effort.py` — Retain all-taxon survey evidence, even when sheephead was absent.
- `code/tools/refit_pisco_fish_effort.R`
- `code/tools/report_pisco_fish_effort.py`
- `code/tools/verify_pisco_fish_effort.py`

### `outputs/pisco_fish_rules_20261101/`

- `code/tools/build_pisco_fish_rules.py` — Authorized PISCO measurement policy, isolated from manuscript/canonical outputs.
- `code/tools/refit_pisco_fish_rules.R`
- `code/tools/report_pisco_fish_rules.py`
- `code/tools/verify_pisco_fish_rules.py`

### `outputs/pisco_painted_cave_20261102/`

- `code/tools/audit_painted_cave_sizes.py`

### `outputs/pisco_pass_crane_20261101/`

- `code/tools/audit_pisco_pass_crane.py` — Exact legacy paired support; do not count every mapped site as current inclusion.
- `code/tools/compare_pisco_pass_policy.py` — Add one isolated pass-completeness contrast on identical years, and one explicit coverage contrast for CRANE.
- `code/tools/refit_pisco_pass_policy.R`
- `code/tools/report_pisco_pass_policy.py`
- `code/tools/verify_pisco_pass_policy.py`

### `outputs/pisco_size_zone_resolution_20261102/`

- `code/tools/report_pisco_size_zones.py`
- `code/tools/verify_pisco_size_zones.R`

### `outputs/pisco_source_reconciliation_20261101/`

- `code/tools/check_reconciled_pisco.py`
- `code/tools/reconcile_pisco_sources.py`
- `code/tools/refit_reconciled_pisco.R`
- `code/tools/report_reconciled_pisco.py`

### `outputs/pisco_timing_donors_20261102/`

- `code/tools/report_pisco_timing_donors.py`
- `code/tools/verify_pisco_timing_donors.R`

### `outputs/point_dume_coverage_20261101/`

- `code/tools/audit_point_dume_coverage.py` — Keep raw rows only at their canonical source; export bounded derived summaries.
- `code/tools/check_point_dume_coverage.py`
- `code/tools/pool_coverage_rules.R`
- `code/tools/report_point_dume_coverage.py`

### `outputs/point_vicente_seasons_20261101/`

- `code/tools/audit_point_vicente_seasons.py` — A present but unmeasured record invalidates its response/event; absent taxa are zeros.
- `code/tools/check_point_vicente_date_pairs.py`
- `code/tools/fit_point_vicente_date_rules.py`
- `code/tools/refit_point_vicente_seasons.R`
- `code/tools/report_point_vicente_seasons.py`
- `code/tools/verify_point_vicente_seasons.py`

### `outputs/pooled_common_estimand_20261101/`

- `code/tools/bridge_common_estimand.R` — Evaluate only the pool function, not the analysis script's top-level code.
- `code/tools/refit_common_estimand.R`
- `code/tools/summarize_common_estimand.py` — Recommended structure

### `outputs/required_scope_area_20261102/`

- `code/tools/audit_area_missingness.py` — Reproduce old zero-then-remove rule and compare explicit missing support without imputing zero.

### `outputs/researcher_scope_20261101/`

- `code/tools/build_researcher_scope_review.py` — Apply user-directed researcher/ecological exclusions before review-only fitting.
- `code/tools/compare_diver_exclusions.R`
- `code/tools/summarize_researcher_scope.py`

### `outputs/satellite_baseline_comparison_20261101/`

- `code/tools/check_satellite_baseline_comparison.py` — Independent NumPy regression and calendar-HAC reconstruction; no primary writes.
- `code/tools/compare_satellite_baselines.R` — Run from repository root; isolated outputs only. See output PLAN.md.
- `code/tools/prepare_satellite_baseline_comparison.py` — Refresh isolated 16-pair inputs; never chooses sites by fitted response.

### `outputs/satellite_primary_20261101/`

- `code/tools/summarize_satellite_primary.py` — Summarize the frozen expanded satellite refit without writing canonical results.

### `outputs/seasonal_zone_audit_20261101/`

- `code/tools/audit_seasonal_zones.py`
- `code/tools/check_seasonal_zones.py`
- `code/tools/check_short_history_influence.R` — Diagnostic only: remove the one previously identified three-year history.
- `code/tools/refit_other_program_units.py` — Survey-unit reconstruction and conservative fixed-unit patch; no canonical writes.
- `code/tools/refit_seasonal_zones.R`
- `code/tools/report_seasonal_zones.py` — Seasonal loss counted only in current candidate histories' actual paired years.

### `outputs/site_assignment_resolution_20261101/`

- `code/tools/resolve_unmatched_site_assignments.py` — Resolve inventory-only unmatched assignments through explicit source-map hierarchy.

### `outputs/site_limitations_current/`

- `code/tools/review_current_site_limitations.py` — Summarize current support and outstanding work without changing selection or estimates.

### `outputs/station_coverage_refits_20261101/`

- `code/tools/report_station_coverage_refits.py`
- `code/tools/verify_station_coverage_refits.R`

### `outputs/upstream_sites_20261101/`

- `code/tools/build_upstream_cohort.py` — Refresh the formerly excluded PISCO urchin branches with the restored size pool.
- `code/tools/check_upstream_raw_fields.py`
- `code/tools/plot_upstream_histories.R`
- `code/tools/reconstruct_upstream_sites.py` — Only remove named whole-MPA omissions. Other expert reference exclusions remain.
- `code/tools/refresh_upstream_pisco_sizes.R`
- `code/tools/summarize_upstream_sites.py`

### `outputs/vrg_source_search_20261102/`

- `code/tools/inspect_vrg_workspace.R`
- `code/tools/reconcile_vrg_archive.py`

## Shared helpers and other scripts

Imported by other scripts or referenced only from code.

- `code/audit_biomass_singleton_refit.py` — Compare a completed singleton-biomass refit with its preserved checkpoint.
- `code/audit_fixed_kelp_refit.R` — Run after the SBC127.5 refit; baseline path is an explicit frozen audit input.
- `code/audit_kelp_calibration.py` — Verify the current SBC 127.5 calibration and the historical release chain; never write analysis outputs.
- `code/audit_kfm_site_continuity.py` — Audit actual KFM predator count-site membership; no primary outputs modified.
- `code/audit_kfm_size_feasibility.py` — Audit KFM predator sizes against retained density histories; never promote biomass.
- `code/audit_landsat_reference_refit.py` — Compare the corrected run with the dated pre-correction checkpoint.
- `code/audit_lter_count_refit.py` — Read-only audit of the authorized LTER sheephead COUNT correction.
- `code/audit_monthly_kelp_refit.py` — Audit the monthly refit against its preserved pre-change input/model snapshot.
- `code/audit_naples_inclusion.py` — Audit one authorised membership addition against the pre-inclusion checkpoint.
- `code/audit_primary_uncertainty.py` — Independently check uncertainty exports and preservation of primary results.
- `code/audit_reference_assignments.py` — Build an assignment ledger without changing sites or fitting models.
- `code/audit_sheephead_allometry_refit.py` — Compare a sheephead-only correction with its preserved input/model checkpoint.
- `code/promote_fixed_kelp.py` — Promote only kelp biomass from a staged full harmonization; preserve other rows.
- `code/promote_lobster_missing_lengths.py` — Promote only LTER lobster biomass from a staged full harmonization; preserve other rows.
- `code/promote_pisco_baseline.py` — Rebuild the authorized six fixed-station PISCO inputs after module 59.
- `code/tools/apply_kfm_reference_review.py` — Expose completed reference tests without adopting their omissions.
- `code/tools/audit_calendar_transition.py` — Compare the surveyed-calendar primary policy with its frozen preceding fit.
- `code/tools/audit_selection_coverage.R` — Reconstruct survey-site continuity without changing abundance or model inputs.
- `code/tools/audit_selection_criteria.py` — Audit observed coverage and delivered pixel overlap; never alter primary data.
- `code/tools/audit_survey_sources.py` — Audit source labels and co-located survey dates without changing inputs.
- `code/tools/audit_tom_delivery_2026_10_20.R` — Diagnostic reconstruction of Tom Bell's delivered arrays; no primary writes.
- `code/tools/build_coauthor_decision_review.py` — Build KFM site-map data for the current analysis decisions page.
- `code/tools/build_coauthor_map_layers.py` — Extract review-map layers from the original satellite and MPA source files.
- `code/tools/build_methods_results_reference.py` — Generate site and response tables for the current analysis page.
- `code/tools/check_joint_biomass_diagnostics.py` — Identification, influence and conditional simulation checks for joint review fits.
- `code/tools/check_kfm_sheephead_biomass.py`
- `code/tools/check_pisco_primary_policy.py` — Independent scope, annual-value, period and primary-route checks for the PISCO rule.
- `code/tools/check_pisco_release.py` — Compare public release 9 with local inputs without replacing analysis data.
- `code/tools/check_processing_labels.py` — Export retained processing routes without relabelling or refitting primary models.
- `code/tools/check_readme_recovery.py` — Reproduce the two README recoveries and audit cleanup evidence, without model runs.
- `code/tools/check_scott_reconciliation.py` — Independent keyed checks of review scope, mass corrections and current guide.
- `code/tools/compare_naples_report.R` — A report-style calendar trend on current inputs, not a report reproduction.
- `code/tools/consolidate_current_guide.py` — Final documentation overlay; no eligibility flags or numerical results change.
- `code/tools/document_workflow_tables.py` — Current documentation tables: explicit workflow scope, no numerical adoption.
- `code/tools/finalize_joint_biomass_review.py` — Verify review outputs, summarize sensitivity, then refresh report and guide.
- `code/tools/inventory_viewer.py` — Focused current-decision viewer; historical audits stay in their original files.
- `code/tools/reconcile_satellite_delivery.R` — Diagnostic reconstruction only; never replaces the delivered primary series.
- `code/tools/repair_upstream_size_assignments.R` — Add source-mapped Anacapa SMCA size memberships absent from the original metadata membership join.
- `code/tools/report_biomass_next_steps.py`
- `code/tools/report_joint_borrowed_biomass.py` — Build a review report without promoting manuscript or primary estimates.
- `code/tools/report_nonlinear_biomass_calibration.py` — Report paired nonlinear calibration separately from primary inference.
- `code/tools/report_whole_dataset_biomass_update.py` — Review-only downstream summaries from the integrated source-to-model build.
- `code/tools/summarize_integrated_review.py` — Plain-language review notes derived from the completed integrated tables.
- `code/tools/test_biomass_observation_model.py` — Independent algebra, estimability and support checks for conditional diagnostics.
- `code/tools/test_borrowing_support.py` — Independent regression tests for portable rebuild and generated CSVs.
- `code/tools/test_joint_borrowed_biomass.py` — Read-only independent checks. Run: python3 -B test_joint.py [--repo PATH].
- `code/tools/test_nonlinear_biomass_calibration.py` — Independent algebra, generator-moment and scope checks for nonlinear review.
- `code/tools/update_biomass_next_steps_guide.py` — Documentation overlays for verified next-step review outputs, no primary promotion.
- `code/tools/update_integrated_guide.py` — Final inventory overlay. Numerical review results never edit frozen paper data.
- `code/tools/update_joint_biomass_guide.py` — Final current-guide overlay for the verified joint working-model review.
- `code/tools/update_nonlinear_biomass_guide.py` — Current guide overlay for nonlinear calibration and verified KFM delivery.
- `code/tools/update_scott_reconciliation_guide.py` — Final overlay: verified source corrections supersede receipt-only holds.
- `code/tools/verify_direct_biomass_process.R` — Independent conditional GLS/REML checks; does not run or modify primary analyses.
- `code/tools/verify_landsat_v2.R` — Reconstruct the corrected masks independently from the cached NetCDF.
- `code/tools/audit_source_metadata.py` — Recover original email workbooks and compare them with the delivered Landsat CSV. (documented in `data/external/source_metadata_2026-10-16/README.md`)
- `code/tools/check_satellite_pixels.R` — Read-only check of the recovered 2022 pixel workbook and linked NetCDF. (documented in `data/external/source_metadata_2026-10-16/README.md`)
- `code/tools/audit_pisco_baseline_policy.R` — Decision audit only: no primary input, classification or manuscript is rewritten.
- `code/tools/check_fish_lter_reconciliation.py` — Verify source hashes, unchanged primary outputs, and comparison invariants.
- `code/tools/check_integrated_mass_functions.py` — Independent R/Python conversion agreement at representative sizes.
- `code/tools/check_parallel_resolution.py` — Coordinator checks of load-bearing parallel-audit claims; no scientific writes.
- `code/tools/check_satellite_baseline_policy.py` — Verify the primary cohort implements the satellite baseline policy.
- `code/tools/run_archived_size_hybrids.py` — Two matched-cohort joint refits separating size-pool point and sampling covariance.
- `code/tools/run_archived_size_matched_refit.py` — Review-only same-cohort model refit, current versus Emily archived size pool.

## Archived (unreferenced)

Historical one-off scripts for the diver-exclusion, common-estimand and KFM-holdings audits (Nov 1–2). Not maintained.


- `code/tools/archive/build_diver_exclusion_audit.py` — Harris omitted ratios escaped the current raw-data protocol/calibration update.
- `code/tools/archive/inventory_kfm_holdings.py`
- `code/tools/archive/plot_common_estimand.R`
- `code/tools/archive/rebuild_common_kfm_kelp.py`
- `code/tools/archive/summarize_diver_exclusions.py`
