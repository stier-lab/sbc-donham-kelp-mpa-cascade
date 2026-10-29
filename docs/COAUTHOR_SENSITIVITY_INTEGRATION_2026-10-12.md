# Coauthor sensitivity integration

## Scope

Module `code/R/34_coauthor_sensitivities.R` runs in an isolated R process before the results summary in `run_all.R`. It reads tracked annual responses, primary effect estimates and the reference registry. It does not edit raw observations, exclusions, reference assignments or primary effect estimates. Contributor questions remain unresolved; these checks do not validate ecological habitat matching.

The full pipeline completed on 12 October 2026 in 12.31 minutes. Log: `logs/pipeline_20261012_110210.log`. All 144 primary histories and the primary pooled, pathway and modern kelp-resilience tables are identical to the pre-integration copies. Existing diagnostic limitations remain in their exports.

## Gull Island KFM

The four retained KFM histories use declared reference `a-k-10`. The refits omit either all four histories or only lobster density; the distinct PISCO pairing is retained. Each refit uses the same joint REML models, MPA/source random intercepts, t inference and a nine-response BH family as the primary synthesis. The retained-baseline coefficients reproduce the primary export within 1e-7. All nine coefficient signs and adjusted detection classifications are unchanged under both omission rules.

The kelp pooled log contrast changes from 0.807 to 0.771 when all Gull KFM histories are omitted; lobster density changes from 0.951 to 0.994 and remains inconclusive. Complete results: `tables/table_s_coauthor_gull_omission.csv`; membership: `tables/table_s_coauthor_gull_membership.csv`. Scott's assessment of reference suitability is still needed before changing the primary cohort.

## Combined predators

The exploratory check sums standing lobster and sheephead responses before forming lnRRs. It requires primary support for both predators and purple urchins in the same programme/MPA/response and finite inside/reference observations in at least four common after-designation years. Species-specific lnRRs are not summed. Annual ratios use a common half-minimum nonzero offset within each cohort. Temporal means are taken within source and then averaged equally across contributing sources within each MPA.

Biomass has only two MPAs and is not tested. Density has ten MPAs: Spearman rho = -0.091, permutation P = 0.810 from 19,999 seeded permutations plus the observed ordering. This is a descriptive association of post-designation levels, not a designation effect or a measure of consumptive capacity. Shared references, estimated predictors and different source periods limit inference. Temporal complementarity between predator species is not measured.

Exports `table_s_coauthor_combined_predator_*` retain annual inputs, source-level coverage, MPA-level contrasts and tests. Input checksums are in `table_s_coauthor_sensitivity_inputs.csv`.

## Kelp density and taxonomy

The harmonized kelp density route represents adult stipes, not comparable juvenile/adult plant densities across programmes. `table_s_coauthor_kelp_density_coverage.csv` records available processed histories; it is not evidence that original programmes collected no stage-specific data. A juvenile-versus-adult plant analysis would require new raw-data harmonization and contributor confirmation. Absence of satellite density is not the reason to omit a diver-only test.

Sheephead display names are now Bodianus pulcher and B. pulcher, checked against WoRMS AphiaID 1702292. `data/taxon_display_names.csv` records the synonym. Original Semicossyphus/S. pulcher identifiers remain in data, joins, model terms and history keys. `code/R/taxon_display.R` changes display labels only.

## Manuscript connection and checks

The paired core supplement reports these checks in S5.4, Tables S23-S24. The companion manuscript points to that section after its disease discussion. The reference-site lookup covers all 144 retained histories; the two atlases retain all primary histories and all 36 selected pBACIPS fits.

Run `Rscript code/R/34_coauthor_sensitivities.R` and `Rscript code/tests/test_coauthor_sensitivities.R`. Fixtures cover scaling invariance, all-zero responses, common-year coverage, duplicate annual keys, exact small-sample permutations and stable display aliases. The module asserts unchanged input checksums and reproduces the primary baseline before exporting omissions.
