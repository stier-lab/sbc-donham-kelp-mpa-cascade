# Current decision and dataset guide

Updated November 2, 2026. This is the current entry point. These tables describe review analyses; they do not freeze or replace the manuscript analysis.

- [Current rules](current_rules.csv): scope, implementation, limits and superseded choices.
- [Comparison decisions](comparison_decisions.csv): all 206 candidates, with 175 retained review histories and 31 researcher/ecological exclusions. Site inclusion and model eligibility are separate. The paper-membership column describes the existing paper, not a new selection.
- [Measurement support](measurement_support_current.csv): all 94 biomass histories, including corrected KFM donor identities and the latest six-comparison borrowing review.
- [Visit-level support](measurement_visit_decisions.csv): 1,086 PISCO visit/zone/species uses; measured sizes, donor identifiers, matching gaps, calibration support and unsupported observations. These are not independent samples or comparison-years.
- [Decision history](decision_history.csv): all 51 dated entries preserved, with earlier stages explicitly identified.
- [Open questions](unresolved_issues.csv) and [freeze readiness](freeze_readiness.csv): work still required before the primary analysis can be frozen.
- [Interactive tables](index.html) and [full master table](audit_master.csv).

For the six PISCO animal comparisons, the current review retains 63 post-protection density comparison-years. Conditional borrowing supports 44 biomass comparison-years. The earlier 25 local-only years used different temporal-support rules, so 44 versus 25 is not a pure borrowing effect. Of 1,086 visit uses, 624 have same-visit measurements, 131 borrow from another site, 8 from another visit, 285 have zero counts and 38 remain unsupported. Those 38 uses represent 24 site visits; 26 uses fall after protection and block 19 comparison-years. Borrowing is reference-only; protected-site prediction has not been validated. Missing biomass is never converted to zero, and unsupported reference stations are not silently dropped.

Six other biomass histories (four KFM and two PISCO histories using KFM sizes) remain provisional because intentionally excluded KFM donors entered under alternate labels. Density remains usable. Historical and recent original sizes are now reconciled in the corrected six-history review. Conditional uncertainty checks do not yet represent total uncertainty.

## Reproduce and verify

From the analysis repository, run `python3 code/tools/report_dataset_inventory.py` to rebuild these documentation tables from existing audit outputs, then `python3 code/tools/verify_consolidated_inventory.py` and `python3 code/tools/check_current_decisions.py`. The renderer calls `consolidate_current_decisions.py` last so older summaries cannot overwrite current decisions. `current_decision_inputs.json` records SHA-256 hashes of its immediate evidence inputs.

The six-history scientific review has its own entry point, `python3 code/tools/run_pisco_biomass_review.py`; see [its methods and limits](../pisco_biomass_borrowing_20261102/README.md). Earlier all-program station, source and uncertainty audits retain their scripts and dated outputs. This consolidation does not rerun their analyses.

## Parallel follow-up audit

[Four-area findings and recommendations](../parallel_resolution_20261102/README.md) distinguish verified corrections from proposals that still need review. The six-comparison density review uses the biomass roster: all 63 post-protection years remain, but six earlier density comparison-years are omitted. The next-step rebuild restores those density records; full-zone versus shared-zone baseline availability is now explicit. Before–after model eligibility remains separate. Temporal borrowing validation contains eight predictions at two stations in one year, for red urchins and lobster only; no purple-urchin temporal validation is available. Recommendations do not change eligibility or models.

## Archived overviews and retained evidence

Stale inventory overview documents are preserved in `docs/archive/20261102_pre_consolidation/`. Dated scientific audit directories remain at their original paths so source links and reproducibility are preserved. Their findings describe those stages; use this guide and the current rules for present decisions. No manuscript, canonical result or raw-data files changed in this consolidation.

## Implemented follow-up

[Next-step implementation and results](../biomass_next_steps_20261102/README.md): restored density baseline; lossless donor attempts; uniform timing tolerances; exact KFM year holds; conditional size/count models and a separate direct-measurement process diagnostic. Strict borrowing remains44/63;7/14day sensitivity45/63. No new primary policy or manuscript changes.

## Joint biomass recovery review

[Same-year results, figure and model limits](../joint_borrowed_biomass_20261102/index.html). Measured sizes, shared size effects and recovery through time are now fitted in one working likelihood. Six histories retain their existing support; 45 comparison-years in the seven-day sensitivity, 44 with the strict rule. Slopes remain similar and intervals widen; all six current and joint intervals include zero. This is a first-order conditional model, not total uncertainty or a primary-analysis change. Interval calibration remains provisional; use the linked diagnostic results before interpreting interval coverage. KFM original-size holds remain.

## Nonlinear calibration and source receipt

[Current results and figure](../nonlinear_biomass_calibration_20261102/index.html) extend the joint-model review with exact biomass aggregation and sample/count resampling. Same six PISCO histories and45post-protection years; no paper or primary update. Calibration remains conditional, not a complete nonlinear likelihood. Scott’s original urchin-size workbook has arrived and been verified. It covers the requested donor years, and all 54 donor identities are now verified in the subsequent source reconciliation; other support and uncertainty limits remain. Use the current KFM year-hold and receipt tables rather than older waiting-for-file notes.

## Completed KFM source reconciliation

[Corrected results and figure](../kfm_scott_reconciliation_20261102/index.html) resolve all54recent donor records, preserving physical-site exclusions. Together with45historical corrections and one verified zero, all100affected post-protection records now have resolved source identity. Six histories were refitted on unchanged years. Density, primary and paper remain unchanged. The KFM year table separates cleared alias issues from remaining support limits: cross-program sizes, sparse samples and incomplete total uncertainty.

## Source safeguards and older records

[Latest source audit and pipeline safeguards](../kfm_pipeline_guard_20261102/README.md) extend identity checks to904older groups and fix explicit Anacapa memberships. The initial audit identified two older discrepancies and four retained 2023 groups. See the current source-resolution outcome below; earlier six-history review corrections remain valid and the paper is unchanged.

## Remaining source checks: current outcome

[Source resolution and effect checks](../kfm_remaining_sources_20261102/index.html): Hare Rock 1987 is corrected using 102 legacy animals corroborated by all 20 published bins. Pelican Bay 1995 remains unresolved; either source version preserves the positive Scorpion primary before–after result (1.253 vs. 1.270). Four 2023 groups have zero uses in current review/primary biomass histories. Post-protection trends are unchanged. Five raw groups remain unverified; paper unchanged.

## Full comparison with original measurements

[Independent raw and corrected input comparison](../kfm_source_comparison_20261102/index.html) checks all 1,482 red/purple KFM groups: 903 of 904 historical and all 558 Scott-covered groups match after corrections. Pelican 1995 remains unresolved; 20 raw 2023 groups lack source coverage, of which four survive exclusions and none is currently used. The 288 source-only groups are listed separately; availability does not override eligibility. Other taxa are explicitly outside this measurement check.

## Source-first production route

[Raw-to-biomass implementation and current holds](../raw_source_pipeline_20261102/index.html): 23 immutable delivery snapshots; original KFM and UCSB swath size adapters; source selection and biomass donor ledgers. LTER is now integrated; full harmonization remains blocked on VRG source availability and unresolved KFM groups. This replaces opaque combined-file fallback, not the primary results. Missingness, sampling units, area normalization and complete uncertainty remain pre-promotion checks.

## LTER originals located on the NAS

[Verified source comparison](../lter_original_sizes_20261102/README.md): all 2,706 focal LTER rows match counts and sizes in SBC 19.33; four Arroyo Quemado2022 rows differ only in day. Both NAS releases19.32/19.33 pass official checksums. Missing LTER source availability is resolved; the adapter now preserves group-average size semantics and all 897 missing-size rows. Cached UCSB PISCO measurements are reconciled; the subsequent original-source adapter corrects the site alias and preserves the Painted Cave exclusion.

## Current LTER integration and PISCO verification

[Source comparison, changes and remaining holds](../lter_integration_pisco_check_20261102/index.html). LTER original-source intake is implemented: 2,706 positive-count rows, 1,809 supported group means and 897 missing sizes; four source dates corrected. Annual mass proxies agree to numerical precision. Group-average diameters do not recover individual size variation. All 3,970 older UCSB PISCO roving measurements match the cached source, but 99 Painted Cave East rows have a South Point aggregation alias in the old combined file. The subsequent join reconstruction resolves this assignment: all 99 were unmatched in the old workflow and belong to excluded Painted Cave. Historical UCSB roving sizes are now integrated. The 2,915 additional source rows are not adopted; VRG originals remain missing. Paper unchanged.

## Painted Cave correction and archived biomass recipe

[Correction, retained-input comparison and Emily’s archived workflow](../pisco_painted_cave_20261102/index.html): 99 size records representing139urchins receive the correct Painted Cave alias. They were unmatched in the old three-key join, not assigned to South Point. Correcting them changes no retained donor memberships because Painted Cave remains excluded. All3,970historical UCSB roving measurements are independently verified;332Painted Cave rows remain in the exclusion ledger and3,638other rows enter candidate inputs subject to production reference exclusions. VRG originals remain missing. Archived code/notes confirm cross-program size borrowing by comparison, role, year and species; matching suitability and total uncertainty remain separate checks.

## VRG archive reconciliation

[Search and exact workspace comparison](../vrg_source_search_20261102/index.html): all18,284saved VRG size rows match the combined table. The466remaining rows are the three explicitly excluded reference sites (Graveyard Canyon97; Sutil283; Ship Rock86). All18,750focal VRG aliases match curated physical sites. Both Drive combined copies are byte-identical to our input. Independent original VRG delivery and pre-combination assembly steps remain unrecovered; workspace objects are derived evidence, not new raw sources. No numerical changes or primary promotion.

## Current VRG decision: accepted contributor delivery

Adrian approved Emily’s combined VRG size file on November2. The exact-file/VRG-only acceptance is recorded in [the decision](../vrg_source_search_20261102/acceptance_decision.json). The pipeline retains18,750source rows, excludes466researcher reference rows and accepts18,284candidate donors. Missing original delivery/assembly details remain documented but no longer block VRG use. Separate KFM source and biomass-support/uncertainty checks remain. No paper change.

## Current source scope, Pelican and area checks

[Required-year and aggregation review](../required_scope_area_20261102/index.html): four unused2023KFM holds no longer block current matched histories; source support remains unavailable for extensions. Pelican1995is the sole required source hold. Original140versus combined135differs by five sizes; recommendation is original140with sensitivity, pending source choice. All23,236LTERurchin source areas are1m²; explicit normalization changes cached values by zero. The72missing biomass rows stay missing and lead to the same30unsupported site-year-species combinations. PISCO missing counts/stipe sizes require targeted survey-grid tracing next; no paper update.

## Current source gate: clear; Pelican/PISCO repairs

[Accepted Pelican source and PISCO recovery](../pelican_pisco_recovery_20261102/index.html): Adrian approved original140Pelican measurements;135-animal sensitivity retained. Required-source gate now passes, with four unused2023KFM groups still deferred. PISCO1.11recovers219red urchins at CavernPointWest2013; station density1.319→1.928/m², six-station protected mean3.468→3.569/m² under unchanged aggregation. Kelp stipes remain unknown at LittleDume2011andUnderwaterArch2023; biomass was already protected, stipe density now protected too. Source-presence marker separates explicit unknowns from roster absences. Full review harmonization/model comparison still needed before paper updates.
