# Orientation — California MPA kelp analysis

Updated 20 October 2026. Read `docs/SCIENTIFIC_APPROACH.md` for the current questions, evidence hierarchy and limitations. Older dated audits are historical records, not current result summaries.

## Project and boundaries

This is the public analysis workspace, `sbc-donham-kelp-mpa-cascade`. The sibling `sbc-donham-kelp-mpa-data` owns raw-data processing; the private `sbc-donham-kelp-mpa-manuscript` owns the live Word draft, essential supplement, private correspondence and dated exports. Former Donham-Stier directory names are compatibility symlinks on this Mac. Do not put private manuscript files or publisher PDFs in this public repository.

The manuscript title is **Kelp resilience and trophic responses in a southern California marine protected-area network**. Target: Journal of Applied Ecology Research Article. Adrian is lead/corresponding author and Emily second author. Dan Reed is acknowledged rather than listed as an author, following his preference. Final contributor roles and affiliations still require confirmation.

## Current scientific interpretation

The primary synthesis retains **145 source-specific histories from 21 MPAs**. It combines mean level differences, before/after step changes and fitted model-time changes. Most diver fits use observation indices, not common elapsed-year intervals. Exponentiating the pooled estimate does not produce a common biomass ratio or percentage protection effect.

Three questions govern interpretation:

1. Population patterns: sheephead biomass, purple-urchin density and kelp have pooled directions consistent with some cascade predictions. Red-urchin biomass is positive and may reflect direct protection of a harvested grazer. Lobster responses and the remaining primary coefficients are inconclusive.
2. Linked responses: the lobster–purple-urchin biomass regression uses four MPAs; its nominal P = 0.040 becomes inconclusive under exact permutation (P = 0.125), with shared-reference dependence. Cross-taxa fits do not establish a complete causal cascade or its absence.
3. Disturbance performance and attribution: diver kelp has a baseline-relative advantage inside MPAs during/after the 2014–2016 heatwave. The historical satellite comparison does not resolve whether designation improved that advantage. Higher post-event biomass is not necessarily regrowth after a decline.

The 17-MPA annual-mean before/after kelp bridge is positive, but the six-MPA diver-only interval and descriptive AR1 GLS interval include zero. The latter re-estimates period differences rather than only changing uncertainty. Non-target KFM comparators cover three kelp-matched/four urchin-matched MPAs and leave reef comparability unresolved. See modules 27, 28 and 45 and `docs/SCIENTIFIC_APPROACH.md`.

## Inputs and authoritative outputs

Tracked harmonized inputs and caches support the core workflow without newly distributing restricted raw records. Optional raw-data analyses require their documented local inputs; do not interpret a skipped optional stage as a completed verification.

- Primary contrasts: `data/sumstats_final.csv`.
- Population synthesis: `tables/table_02_meta_analysis.csv`.
- Cross-taxa fits: `tables/table_03_cross_taxa_meta_regression.csv`; exact-test outputs in `tables/table_s_prereview_crosstaxa_permutation.csv`.
- Figure 4 uncertainty: `tables/table_s_temporal_ar1_slopes.csv`.
- Modern kelp disturbance comparison: `tables/table_resistance_recovery_*.csv`, Figure 5 from module 22.
- Matched design bridge: `tables/table_s_design_matched_*.csv`, module 27.
- Historical disturbance comparison: `tables/table_s_historical_heatwave_*.csv`, module 28.
- Full evidence routing: `docs/ESSENTIAL_SUPPLEMENT_INDEX.md`.
- Current generated results summary: `docs/RESULTS_SUMMARY.md`; legacy percentage columns retain their algebraic meaning, not a universal ecological fold change.

## Reproduction

Core entry point: `Rscript code/R/run_all.R`. Cached figure helper: `code/R/run_figures_only.R`; it is not a complete paper verification. The resilience suite is `run_resilience.R`, with module membership in `resilience_modules.R`.

Modules 35–44 cover source corrections, primary/temporal uncertainty, sizes, programme comparisons, biomass reconciliation and sheephead size proxies. Module 45 is an optional descriptive non-target comparability diagnostic retained in the repository but removed from V27 and its supplement at Adrian's request on 21 October 2026 pending discussion with Bob. Source conventions, input fingerprints and named environment overrides are documented alongside each analysis. Restore the recorded R dependencies before a clean run; use fresh test/run logs, not the obsolete September run status, to claim successful reproduction.

## Open scientific decisions

- Corrected satellite masks and the recovered workspace are verified: all 24 MPA and 34 nonempty reference arrays reproduce exactly, including missingness, and all 60 input rows match. See `docs/LANDSAT_V2_VERIFICATION_2026-10-21.md`. Some splits were recovered against the delivered arrays; a unique historical clustering solution is not established.
- Current harmonization averages available reference arrays. Preferred-reference selection is a different, unperformed sensitivity and needs verified assignments and a source-specific rule; fixed pre-2012 selection is not uniformly pre-designation.
- Judith Rock eligibility, San Miguel satellite identity and historical protection footprints still require contributor confirmation.
- KFM fish age-category provenance and shared biomass size pools limit interpretation; independent counts do not ensure independent biomass estimates.
- Gull Island KFM remains primary at Adrian's explicit instruction, with an omission sensitivity; that decision does not validate the reference habitat.
- Some older readiness logs and exported drafts predate the corrections. Do not infer final coauthor approval, source closure or submission readiness from a completed prose pass.

The private manuscript repository contains the current revision and confirmation logs. No new circulation, push or archive publication is implied by local analysis or writing updates.
