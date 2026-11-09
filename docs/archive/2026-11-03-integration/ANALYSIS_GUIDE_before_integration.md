# Current data and analysis guide

Updated November 3, 2026. This is the authority for navigating the current decisions. The linked rule and eligibility tables hold the detailed decisions; dated audit reports preserve the evidence behind them.

## Where the work stands

The required-source check passes. Original Pelican measurements are accepted, the recoverable PISCO count is repaired, and researcher exclusions apply to size donors as well as counts. Four unused 2023 KFM size groups remain unavailable. Passing the source check does **not** mean biomass support, uncertainty or final model selection is settled.

The source-processing changes are implemented in the data repository. Coverage and borrowing alternatives have been tested in separate review analyses. A complete harmonization and coordinated model rerun incorporating all corrections is still needed. The manuscript, supplement and frozen results have not been updated with these audits.

## Which workflow is being described?

| Workflow | Status | What it can establish |
|---|---|---|
| Delivered files and corrected source layer | Implemented | Trace observations, physical identities, units, corrections and exclusions. |
| Inherited biomass production | Existing calculation; integrated corrected rebuild pending | Reproduce the existing pooled-size calculation. Its broad donor pools are not automatically validated by source repairs. |
| Six-history PISCO borrowing and uncertainty analysis | Tested for review | Compare stricter visit/site matching and shared-donor uncertainty within its stated scope. |
| Four-history KFM sheephead biomass analysis | Tested for review: 52 comparison-years | Compare post-protection density and biomass using supported visits and conditional joint uncertainty. |
| KFM 55-year matched-visit option | Coverage and point estimates only; not fitted | Show what could be retained by dropping unmatched count visits. It is not a recovered annual dataset. |
| Manuscript and frozen outputs | Unchanged by these audits | Preserve the existing submitted-to-coauthors analysis until a coordinated replacement is reviewed. |

The [workflow table](../outputs/dataset_inventory_20261101/workflow_status.csv) gives each entry command, working directory, input record, prerequisites, output and adoption limit. Run paths marked `.` from this analysis repository; data-repository commands run from the sibling data repository. Private inputs and the NAS are not bundled with the public code. “Implemented” refers to code behavior; “review” means a separate calculation, not an adopted primary analysis. “Unsupported biomass” does not mean that the site or its density data are excluded.

## Follow the data

| Stage | What happens | Record to inspect |
|---|---|---|
| Contributor delivery | Preserve the received bytes, release, units and access conditions. Provider-integrated tables and contributor-derived satellite products are labeled as such. | Data repository `config/source_manifest.json`; immutable `data/raw/` snapshots |
| Source reconciliation | Recover physical identities, retain original row IDs, convert units and record keyed corrections. Never overwrite the delivery. | Data repository source selection, alias, missing-value and count-recovery ledgers |
| Site and response selection | Assign protected/reference membership separately from physical site identity. Apply researcher advice and ecological scope to each response. | [Comparison decisions](../outputs/dataset_inventory_20261101/comparison_decisions.csv) |
| Survey coverage | Check dates, depths, zones, stations, methods and effort before combining surveys. | [Site coverage](../outputs/dataset_inventory_20261101/site_response_coverage.csv), [annual survey records](../outputs/dataset_inventory_20261101/raw_station_year_audit.csv) |
| Density and biomass | Preserve measured counts; distinguish direct sizes, group means, borrowed sizes, imputation and unsupported biomass. Normalize by the surveyed area and effort. | [Measurement support](../outputs/dataset_inventory_20261101/measurement_support_current.csv); data repository biomass traces and donor pools |
| Temporal analysis | Separate before–after changes from post-protection trajectories. Assess model eligibility independently of site inclusion. | [Current rules](../outputs/dataset_inventory_20261101/current_rules.csv), [freeze checklist](../outputs/dataset_inventory_20261101/freeze_readiness.csv) |
| Release | Refit all downstream analyses, verify agreement and show changes before updating the paper. | [Remaining questions](../outputs/dataset_inventory_20261101/unresolved_issues.csv) |

The [technical build guide](../../sbc-donham-kelp-mpa-data/docs/RAW_TO_BIOMASS.md) gives commands, files, formulas and safeguards. The [searchable inventory](../outputs/dataset_inventory_20261101/index.html) brings the linked tables together.

## Apply inclusion decisions in this order

- **Ecological inclusion:** preserve researcher-advised and ecological-scope exclusions, including their program and response scope. Gull Island PISCO is distinct from excluded Gull Island KFM; Anacapa SMR is distinct from Anacapa SMCA. Do not exclude a comparison because its baseline slope, effect direction or significance is inconvenient.
- **Reference membership:** use the documented physical stations and roles. An alias is not a new station. Established assignments with incomplete habitat documentation remain provisional and receive influence checks. Missing documentation does not establish a poor match, and an established assignment does not establish equivalent habitat.
- **Coverage:** use shared comparable zones. Constant-zone and constant-station alternatives have been checked, but no universal final station policy is adopted. Actual depths must not be inferred from identical zone labels across programs.
- **Measurement support:** a valid count can support density even when sizes cannot support biomass. A missing size removes support for that response and sampling unit; it is not automatically a reason to exclude the whole site or other taxa.
- **Model eligibility:** before–after estimation needs usable before data and defensible assumptions. A retained site is not automatically eligible for pBACIPS. Baseline-stability tolerances remain diagnostic; no final tolerance or model roster is frozen. Post-protection trends describe relative trajectories, not uniquely attributable protection effects.

The [inclusion-gate table](../outputs/dataset_inventory_20261101/inclusion_gates.csv) records the unit, rule and consequence at each stage. Keep a physical-site exclusion separate from a missing visit, an unsupported response and an ineligible model. Apply ecological exclusions to observations, donors and uncertainty-calibration records within their documented response scope.

The 206 rows in the comparison table are candidate response histories, not 206 independent sites or a frozen primary cohort. Its `exclude` and `reason` fields preserve 31 researcher/scope exclusions. `in_current_manuscript` describes the existing paper; it is not approval of a revised primary analysis. Rule IDs, measurement status and model decision must also be read before using a history.

## Timing, depth and survey effort

For the six reviewed PISCO animal comparisons at Abalone Cove and Point Vicente, the review uses all valid survey months, shared INNER/MID zones and equally weighted repeated visits. There is no adopted protected/reference date-gap cutoff. These rules must not be described as implemented for every dataset.

For audited PISCO fish histories, retain CRANE surveys, average repeat surveys and omit incomplete passes from the combined fish response. Missing fish levels are not zero counts. For references comprising several stations, compare available-station and constant-station results before finalizing aggregation; otherwise station turnover can resemble biological change.

Use actual calendar years and observed follow-up. Do not silently impose a common nine-year window or fill unsurveyed years. Network dates, earlier protection, designation-year observations and current take restrictions are separate metadata fields; none substitutes for the others.

## What a biomass estimate means

| Size information | Treatment | What remains uncertain |
|---|---|---|
| Individual measurements or size-frequency counts | Convert each measured size to mass; retain frequency weights and original row IDs. Divide totals by the appropriate effort and area. | Sampling, counts, conversion coefficients and dependence |
| LTER survey-group mean diameter | Retain the recorded group mean and count; label the nonlinear mass conversion as a proxy. | A mean diameter cannot recover the mean of individual masses or within-group size variation |
| Sizes from a different visit/site/program | Record the recipient and every eligible donor, matching rule and failed alternative. | Representativeness and shared donor/prediction errors |
| Five audited LTER lobster records lacking length | Use same-site/date mean individual mass as the documented imputation; preserve missing original lengths and donor identities. | The imputed mass is not a measured length |
| Positive count without supported size | Keep the count and mark biomass unavailable. | Do not manufacture zero biomass or average around missing stations without an explicit policy |
| Explicit observed zero count | Retain a zero point estimate for the surveyed unit. | Zero observations do not prove absence or zero uncertainty. An unsurveyed unit or explicit unknown remains missing. The KFM joint review still needs size support at zero-count visits. |

### Calculation and aggregation order

For measured sizes, convert each size to mass first. The mean mass is `sum(frequency × individual_mass) / sum(frequency)`. Applying a nonlinear equation to the mean length generally gives a different answer. LTER’s recorded group-mean diameters are explicitly labeled proxies because individual diameters are unavailable.

| Taxon/measurement | Conversion | Units |
|---|---|---|
| Lobster carapace length | `0.001352821 × length^2.913963113` | mm to g per animal |
| Red urchin diameter | `0.00059 × diameter^2.917` | mm to g per animal |
| Purple urchin diameter | `0.00059 × diameter^2.870` | mm to g per animal |
| Sheephead total length | `0.0144 × length^3.04` | cm to g per animal |
| Adult kelp stipes | Survey-month SBC 127.5 Model1 slope × stipes × 1000 | kg/stipe to g dry foliar mass |

The [technical guide](../../sbc-donham-kelp-mpa-data/docs/RAW_TO_BIOMASS.md) identifies the implementation and source-to-mass traces. These are fixed conversion equations, not estimates of their coefficient uncertainty.

For separate counts and size samples, calculate `biomass_g_m2 = density_per_m2 × mean_mass_g` **within each matched visit**. Then combine visits, zones and stations using the declared workflow weights. Do not multiply annual average density by annual average mass: that loses the association between abundance and size. If the same size survey supplies multiple count visits, it is one measurement pool used several times; report unique measurements separately from uses and share its uncertainty across recipients.

Area-normalized biomass is calculated once. Do not divide density × mean mass by area again. For display, g/m² × 0.1 gives kg/100 m²; individuals/m² × 100 gives individuals/100 m². Preserve dry versus wet mass and adult versus juvenile definitions.

The [calculation table](../outputs/dataset_inventory_20261101/biomass_calculations.csv) distinguishes measured sizes, grouped means, imputation, borrowing and missing values. Each trace should retain the source file/hash and row, physical site, date, response, count state, size unit, area/effort, donor identity, aggregation weight and exclusion/support reason. Missing fields stay unknown rather than being inferred from labels.

KFM per-plant kelp mass requires plant-density information. Diver dry foliar mass and satellite wet canopy mass are different responses. Area units and aggregation order must be retained with each estimate.

## Borrowing: implemented workflow versus tested alternative

**The inherited production pool** selects sizes by comparison, protection role, year and species, with a 25 mm minimum for PISCO urchin donors. It can cross physical stations, dates, zones and programs. The new source layer and donor ledgers make that behavior traceable; they do not by themselves validate it or replace it with stricter matching.

**The six-history PISCO review** uses same-visit sizes first. Other-visit and other-site attempts are recorded separately. Borrowing stays within species, MPA, role, zone and year; cross-site donors and calibration exclude the target site. Actual accepted borrowing is reference-only; protected-site borrowing has not been validated. The lossless attempt table retains rejected alternatives as well as selected donors.

The strict review requires a date gap no larger than its calibration maximum and at least two calibration site-year groups. Two groups is a computational minimum, not proof of adequate biological validation. Small upper-gap tolerances have been tested, with no final primary cutoff selected. Under the strict scenario, 44 of 63 post-protection comparison-years have complete biomass support; the seven-day sensitivity supports 45. These figures describe six histories, not all biomass data.

[Visit decisions](../outputs/dataset_inventory_20261101/measurement_visit_support_detailed.csv) distinguish same-visit measurements, zero counts, borrowed visits, borrowed sites and unsupported values. [Donor attempts](../outputs/dataset_inventory_20261101/donor_attempts_current.csv) retain date gaps, sample sizes, calibration groups, donor keys, selection and rejection reasons. Require both protection roles and the retained reference stations; do not silently drop unsupported stations to obtain an estimate.

Conditional resampling of a fixed size pool does not capture total uncertainty. The joint and nonlinear reviews address some shared-donor and prediction error, but a final integrated treatment of donor sampling, counts, allometry, temporal dependence and shared references remains unfinished. Do not describe their intervals as a complete uncertainty analysis.

## KFM sheephead: exact scope of the current review

The [source reconciliation](../outputs/fish_lter_reconciliation_20261102/README.md) found sizes for 282 of 423 count site-years. Those totals are availability checks, not an eligible biomass cohort. Scott’s fish export is retained in the data repository’s private `data/source_archive/fish_lter_reconciliation_20261102/` review cache and pinned by the review input hashes; it has not been silently added to the production source registry. Sizes begin in 2007 and cannot support a before–after biomass analysis of the 2003 designation. Scott confirms that KFM does not measure lobster carapace lengths.

The [restricted biomass review](../outputs/kfm_sheephead_biomass_20261103/README.md) uses the existing four fish-comparison rosters and follows these rules:

- Use same-date, same-physical-site sizes first. Where unavailable, allow the nearest same-site size survey within seven days. Do not borrow across sites or years. The seven-day rule matches counts to sizes; it is not a protected-versus-reference survey-date tolerance.
- Require size support for **every retained count visit**, including zero-count visits, to propagate biomass uncertainty. Keep observed count zeros as point zeros; missing size support does not justify zero uncertainty.
- Multiply density by mean mass within visits. Preserve the existing count-record weights within a site and equal weights among sites in each protection role. A repeated size survey is not a new independent sample each time it is used.
- Keep both roles and all stations in that year’s existing roster. “Complete” means complete for that annual roster, not identical station membership across all years. Constant-station results remain a separate sensitivity.
- Restrict size-variance and count-dispersion calibration to the retained fish cohort. Keep response-specific exclusions: a benthic exclusion does not automatically exclude fish.
- Fit changes in protected minus reference means over actual calendar years. These are post-protection trajectories, not effects uniquely attributable to protection.

This supports **52 comparison-years**: Anacapa 13, Harris Point 13, Scorpion 12 and South Point 14. The [seven unmatched visits](../outputs/kfm_missing_size_recovery_20261103/missing_visit_evidence.csv) affect five additional years. No lengths were recovered from the checked sources. A [55-year matched-visit option](../outputs/kfm_missing_size_recovery_20261103/five_year_options.csv) would omit some count visits while retaining annual site membership; it has not been fitted or adopted. Harris Point 2010 and Scorpion 2023 would still require unsupported borrowing or site removal. Valid density observations are retained throughout.

### What the uncertainty calculation includes

The KFM review propagates count-rate uncertainty and frequency-weighted size-sample uncertainty through the nonlinear biomass product, sharing draws for reused size surveys. Sparse size pools receive an explicit regularization assumption. A working temporal model then estimates additional annual variation, with calendar-correlation, count-prior, size-pooling, count–size dependence and station/visit sensitivities. The [reproduction driver](../code/tools/run_kfm_sheephead_biomass_review.py) rebuilds the variants, figures, tables and checks in order.

These intervals remain conditional on the observation models and estimated sampling covariance. They do not measure allometric uncertainty, detection error or whether roving sizes represent the fish counted on visual transects. All four main biomass intervals include zero, and Anacapa is interval-method sensitive. Do not promote a favorable interval or call this total uncertainty. Scott’s new reference-site concerns also remain unresolved.

## Current source decisions and exceptions

| Dataset or issue | Current decision | Evidence |
|---|---|---|
| KFM Pelican Bay red urchins, 1995 | Use all 140 original measurements. Retain the 135-animal combined-file sensitivity. | [Pelican and PISCO repairs](../outputs/pelican_pisco_recovery_20261102/README.md) |
| Four unused KFM 2023 groups | Deferred, not accepted or imputed. Current required histories do not use them; an extension requires source verification. | [Required scope](../outputs/required_scope_area_20261102/README.md) |
| KFM excluded size donors | Enforce physical-site exclusions through aliases; recent and historical reconciliations completed. Final biomass support and full refit remain separate. | [KFM reconciliation](../outputs/kfm_scott_reconciliation_20261102/README.md) |
| LTER urchin sizes | Original group means integrated: 1,809 supported and 897 missing among 2,706 positive-count rows in the selected scope. No fabricated individual sizes. | [LTER integration](../outputs/lter_integration_pisco_check_20261102/README.md) |
| PISCO Painted Cave | Correct the 99 erroneous aggregation aliases; retain the ecological exclusion. Those records were unmatched in the legacy join, not South Point donors. | [Alias and membership audit](../outputs/pisco_painted_cave_20261102/README.md) |
| VRG sizes | Exact combined-file subset accepted by Adrian. Retain 466 researcher-excluded rows in the ledger and 18,284 candidate donors. Independent original delivery/assembly remains unavailable. | [Acceptance and reconciliation](../outputs/vrg_source_search_20261102/README.md) |
| Cavern Point West red urchins, 2013 | Recover 219 animals from the matching provider record only; preserve both source identities. | [Keyed recovery](../outputs/pelican_pisco_recovery_20261102/pisco_count_recovery.csv) |
| Missing kelp stipes | Little Dume West 2011 and Underwater Arch 2023 remain unknown. Preserve plant counts; do not infer stipes or biomass. | [Missing-value safeguards](../outputs/pelican_pisco_recovery_20261102/README.md) |
| South Point sheephead, 2005 | Three fish have no recovered lengths; retain the missing-size flag and distinguish density from biomass support. | [Original source check](../outputs/pisco_source_reconciliation_20261101/README.md) |
| KFM predators | Full-history density remains primary. The 52-year sheephead biomass calculation is review-only; no KFM lobster lengths are available. | [Restricted biomass](../outputs/kfm_sheephead_biomass_20261103/README.md) |

## Which table answers which question?

| Table | One row represents | Use it for |
|---|---|---|
| [Workflow status](../outputs/dataset_inventory_20261101/workflow_status.csv) | Implemented, reviewed or candidate workflow | Commands, scope, outputs and adoption boundaries |
| [Inclusion gates](../outputs/dataset_inventory_20261101/inclusion_gates.csv) and [biomass calculations](../outputs/dataset_inventory_20261101/biomass_calculations.csv) | A decision stage or measurement case | Rules, formulas, missingness and limitations |
| [Dataset inventory](../outputs/dataset_inventory_20261101/dataset_inventory.csv) | An audited delivery/dataset | Release, units, methods, missing-value conventions and access; this older delivery scan is not the expanding raw-source registry |
| [MPA metadata](../outputs/dataset_inventory_20261101/mpa_metadata.csv) | An MPA | Protection date and its meaning, earlier protection, take restrictions, geographic metadata and evidence |
| [Site coverage](../outputs/dataset_inventory_20261101/site_response_coverage.csv) | Program/site/response-method coverage entry, including explicitly labeled group entries | Years, gaps, before/after coverage, months, depth, zones, method, effort, assignments and missing data |
| [Comparison decisions](../outputs/dataset_inventory_20261101/comparison_decisions.csv) | Candidate program/MPA/taxon/response history | Included/excluded scope, reasons, reference membership, applicable rules and remaining issues |
| [Measurement support](../outputs/dataset_inventory_20261101/measurement_support_current.csv) | Audited biomass history | Sparse/shared pools and support limitations; some summary metrics describe earlier pools, as labeled |
| [Detailed visits and donor attempts](../outputs/dataset_inventory_20261101/measurement_visit_support_detailed.csv) | Recipient visit use (attempts are in their separate linked table) | How a particular reviewed biomass value was measured, borrowed or left missing |
| [Current rules](../outputs/dataset_inventory_20261101/current_rules.csv) | A decision with scope and implementation status | What is approved, tested or pending |
| [Questions](../outputs/dataset_inventory_20261101/unresolved_issues.csv) and [freeze checklist](../outputs/dataset_inventory_20261101/freeze_readiness.csv) | An unresolved issue or release requirement | What still prevents adoption or needs disclosure |

A table row is not automatically an independent sample. A group label is not a physical station; a dataset hash verifies bytes, not ecological suitability. Unknown habitat, depth, coordinates or exact dates stay unknown.

## Next release and historical records

The [notes and correspondence audit](../outputs/notes_email_audit_20261102/README.md) and its [finding ledger](../outputs/notes_email_audit_20261102/findings.csv) preserve contributor explanations. They support decisions; old status statements do not override the current rules above.

SBC50.21 is also a pinned review input, not a replacement production source. The newer LTER adult kelp source contains the same observations through 2023 as our registered input. SBC50.21 combines adult and juvenile kelp; its reported density must not be treated as adult stipe density. Its adult biomass values imply a fixed coefficient, whereas our pipeline uses the observed survey month. Keep our adult foliar response and actual-month conversion separate from that provider aggregate. The linked provider protocol and implied coefficient need clarification before substituting the aggregate. New 2024–2026 observations do not enter the current cohort automatically. These are source and formula comparisons, not new fitted protection effects.

Emily clarified that her historical minimum was three before observations at the same physical stations, assessed separately for fish and benthic surveys, and that she counted the implementation year as before. A stricter calendar-year convention is a deliberate change, not a recovered original rule. Keep the matched station-year counts and model eligibility visible when finalizing that choice. Her old `Sites_Years.xlsx` conveys coverage by cell colors; blank cell values must not be interpreted as absent surveys.

Run the complete review build, inspect donor and missingness ledgers, finalize coverage/model/uncertainty rules, and compare updated density, biomass and downstream results. Resolve permissions before redistribution. Show those changes before altering the paper.

[Archived documents](archive/2026-11-02-consolidation/README.md) preserve superseded instructions. The [document catalog](DOCUMENT_CATALOG.csv) distinguishes archived documents from technical references and hash-bound historical evidence. Old reports remain valid accounts of their own runs; their cohort sizes, holds and recommendations are not current policy. The [audit register](AUDIT_FINDINGS_REGISTER.md) preserves the sequence of decisions without serving as operating instructions.
