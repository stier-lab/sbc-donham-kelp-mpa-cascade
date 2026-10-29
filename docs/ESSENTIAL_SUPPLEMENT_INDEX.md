# Essential supplement and repository results

The coauthor supplement has four sections, six tables and two figures. Full
analysis outputs remain here; removing a display from the supplement does not
remove its code, inputs, diagnostics or interpretation limits.

| Main claim / essential support | Complete repository evidence |
| --- | --- |
| S1 / Table S1: measurements, coverage and exclusions | `tables/table_s_mpa_data_availability.csv`; `tables/table_s_design_reference_sites.csv`; `code/R/01_utils.R`; `data/harmonized/` |
| S1.1 and S2.3: sheephead age definitions and small-fish sensitivity | `tables/table_s_sheephead_yoy_*`; `code/R/44_sheephead_yoy_check.R`; `code/tests/test_sheephead_yoy_check.R`; `docs/SHEEPHEAD_YOY_2026-10-19.md` |
| S2.1 / Table S2: mixed level/change contrasts and model routing | `tables/table_model_selection.csv`; `tables/table_s_processing_routes.csv`; `docs/PRIMARY_PROCESSING_LABELS.md`; `tables/table_s_pbacips_shape_*`; `plots/fig_s_pbacips_selected_histories.pdf`; `code/R/08_effect_sizes.R`, `32_pbacips_shape_display.R`; `code/tools/check_processing_labels.py` |
| S2.2 / Figure S1: individual estimates behind main Table 2 | `plots/fig_s02_forest_plot.png`; `data/sumstats_final.csv`; `tables/table_02_meta_analysis.csv`; `code/R/09_meta_analysis.R` |
| S2.2: sensitivity to near-threshold kelp contrast selection | `tables/table_s_fixed_kelp_selection*.csv`; `code/R/51_kelp_selection_sensitivity.R` |
| S2.2 / Table S3: exact-cohort annual dependence | `tables/table_s_primary_hac_*`; `code/R/36_primary_uncertainty.R` |
| S2.2 / Table S4: same-observation temporal AR1 | `tables/table_s_temporal_*`; `code/R/10_temporal_analysis.R`, `37_temporal_uncertainty.R` |
| S2.3: programme differences, predator sizes and reference choice | `tables/table_s_programme_*`; `table_s_sheephead_programme_*`; `table_s_sheephead_size_*`; `table_s_lobster_size_*`; `table_s_size_refit_*`; `table_s_coauthor_gull_*`; modules 30, 38, 39 and 41 |
| S2.3: change versus post-only differences and shared references | `tables/table_s_design_*`; `code/R/27_design_geography_sensitivity.R` |
| S3.1 / Figure S2: main Table 3 relationships | `plots/fig_s_primary_pathway.png`; `tables/table_s_primary_pathway_membership.csv`; `tables/table_03_cross_taxa_meta_regression.csv`; `tables/table_s_prereview_crosstaxa_permutation.csv` |
| S3.2: programme-separated pathways on identical MPAs and years | `tables/table_s_programme_pathway_*`; `code/R/42_programme_pathway_refits.R`; `code/R/programme_pathway_utils.R`; `code/tests/test_programme_pathway_refits.R`; `docs/PROGRAMME_PATHWAY_REFITS_2026-10-17.md` |
| S3.2: size-filtered predator-urchin relationships | `tables/table_s_size_refit_pathway*`; `code/R/41_predator_size_refits.R`; `docs/PREDATOR_SIZE_REFITS_2026-10-17.md` |
| S3.2: exploratory links, summed predators and disease | `tables/table_s_extended_pathway*`; `table_s_coauthor_combined_predator_*`; `table_compound_pycnopodia_*`; `table_compound_keystone_moderation.csv`; modules 20, 30 and 38 |
| S4.1 / Table S5: Figure 5 paired kelp responses | `tables/table_resistance_recovery_pairs.csv`; `table_resistance_recovery_sensitivity.csv`; `table_figure5_*`; `code/R/22_resistance_recovery.R` |
| S4.2 / Table S6: historical event comparison | `tables/table_s_historical_heatwave_*`; `code/R/28_historical_heatwave_bacips.R` |
| All observed and selected-fit time series | `plots/fig_s_pair_history_atlas.pdf`; `tables/table_s_pair_history_atlas_index.csv`; `plots/fig_s_pbacips_selected_histories.pdf`; modules 29 and 32 |
| Pooling, outliers, zero offsets, distributions and fit diagnostics | `tables/table_s_prereview_*`; `data/model_diagnostics.csv`; `outputs/lmer_*`; `plots/fig_s11_dharma_diagnostics.*`, `fig_s13_lmer_residuals.*`; modules 13 and 25 |

Table-name prefixes without directories refer to `tables/`. Wildcards identify
retained output families, not additional supplementary display numbers. The
manuscript repository retains the longer technical source and a complete
file-by-file `supplement/readability/reproducibility_index.csv`.

## Additional Non-target Comparators

Module `45_nonfished_comparators.R` and `nonfished_comparator_utils.R` implement
the fixed-legacy-reef KFM comparability diagnostic. Plan, limitations and current
results are in `NONFISHED_COMPARISON_PLAN.md` and `NONFISHED_COMPARISON_RESULTS.md`;
outputs are `tables/table_s_nonfished_*` and `plots/fig_s_nonfished_*`.
Tests: `code/tests/test_nonfished_comparators.R`. This analysis was completed
after the V18 circulation and was incorporated into V26. At Adrian's request
on 21 October 2026, it was removed from the V27 manuscript and essential
supplement pending discussion with Bob. The analysis remains in the repository;
it is not currently part of the paper or its reproducibility index.
It is not an additional numbered supplementary display, a validated causal
negative control or a correction to primary estimates.

## Size and legacy-protection checks

Run `Rscript code/R/38_coauthor_programme_checks.R` followed by
`Rscript code/R/39_predator_size_validation.R`; both are wired into `run_all.R`.
Tests are `code/tests/test_coauthor_programme_checks.R` and
`code/tests/test_predator_size_validation.R`.

Module 39 leaves primary data unchanged. Its Catalina check omits Blue Cavern
Onshore, Farnsworth Onshore and Long Point rather than assigning unverified
historical protection ages. All nine directions and adjusted classifications
persist. This does not show that earlier protection is irrelevant or that a
modern polygon equals a former reserve footprint.

Module 39's temporal size summaries are descriptive. Module 41 now fits actual
matched all-size/large-size Figure-4-style models and predator-purple-urchin
regressions, with retained-membership and expanded-reconstruction scopes kept
separate. Run `Rscript code/R/41_predator_size_refits.R` after modules 38/39;
it is wired into `run_all.R`. It records unavailable inputs rather than claiming
a fresh refit when current size inputs cannot be verified. Its tests are
`code/tests/test_predator_size_refits.R`; methods and qualifications are in
`docs/PREDATOR_SIZE_REFITS_2026-10-17.md`. New displays remain repository figures,
not replacements for main Figure 4 or Table 3. They are not consumption tests
or designation-effect estimates. Conditional lobster resampling
uses 1,000 multinomial draws of measured individuals per MPA-status-year pool.
The intervals quantify size-pool sampling alone, not transect counts, survey
clustering, allometric uncertainty or detectability. Single-individual pools
can have zero resampling spread; this is not evidence of known true size.

Expanded lobster biomass overlaps 51 existing PISCO status-year observations.
Following the singleton correction and full-precision allometry reconciliation,
the median absolute difference is 0.02319 g m^-2^; positive-primary ratios range
approximately 0.850-1.568. Matching full-year size pools narrows those ratios to
0.969-1.014. Expanded measurements remain secondary, not replacements for the
primary harmonised release. Module 43 reproduces the primary VRG hierarchy;
complete comparisons are in `tables/table_s_lobster_reconciliation_*` and
`docs/BIOMASS_SINGLETON_REFIT_2026-10-18.md`.

## Remaining source decisions

The subsequent source audit is
`docs/NO_TAKE_SATELLITE_SOURCE_AUDIT_2026-10-16.md`. It recovers native canopy
biomass units and pixel identities, verifies Naples/Campus Point satellite
reference overlap, and identifies six verified empty reference slots encoded
as zero. These six rows have now been removed before reference averaging;
observed zeros in occupied reference slots remain. The corrected satellite
handoff contains 718 paired MPA-years rather than 719. The release marker is
`data/harmonized/landsat_reference_correction.json`; the before/after audit is
`docs/LANDSAT_REFERENCE_REFIT_2026-10-16.md`. The manuscript's V10 reading pair
incorporates this correction; V9 is preserved as a historical export.

The inherited sheephead-only list is not a legal protection class. Current
regulations are documented separately from programme-specific inclusion;
no additional MPA histories were reclassified. Some inclusion rationales and
historical satellite footprints still need contributor confirmation. Corrected
pixel masks and coverage rules now reproduce all 24 MPA and 34 nonempty
reference series exactly; see `LANDSAT_V2_VERIFICATION_2026-10-21.md`.

Gull Island KFM remains primary at Adrian's explicit instruction on 16 October
2026; omission sensitivity remains. Contributor-confirmed de facto no-take
classifications and satellite extraction metadata must not be inferred.

The public 2024 MLPA release is `10.25494/P6/MLPA_kelpforest.9`. Its fish and
swath tables have the same row counts as the local contributor files, but are
not byte-identical or identical CSV records. Local fish data add `LOOKEDFOR`
and change some sites and one non-focal fish count/length record. Swath data
also have record differences; matching row counts do not establish equivalence.
`data/external/pisco_release9/` records metadata, hashes and difference counts.
Reproduce the read-only public-release comparison with
`python3 code/tools/check_pisco_release.py`; it does not replace analysis inputs.

Reference assignment evidence and remaining matching questions: `docs/REFERENCE_ASSIGNMENT_AUDIT.md`; `tables/table_s_reference_assignment_audit.csv`; rebuild with `code/audit_reference_assignments.py`.
