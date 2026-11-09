# Current data and analysis guide

Updated November 4, 2026. This is the authority for navigating the current decisions. The linked rule and eligibility tables hold the detailed decisions; dated audit reports preserve the evidence behind them.

The [current analysis page](../outputs/dataset_inventory_20261101/index.html) combines inclusion decisions, biomass support, time series, pooled results and the site explorer. For KFM areas, the site explorer maps response-specific survey stations, MPA outlines, 2021 satellite kelp detections and Tom Bell's protected/reference extraction pixels where available. The [map-layer builder](../code/tools/build_coauthor_map_layers.py) reads the original 2021 kelp NetCDF, Tom's MPA shapefile and his pixel-mask workbook; [the extracted layers](../outputs/dataset_inventory_20261101/map_layers.json) record the sources and hashes. The kelp marks are detections on a 30 m grid, not photographic imagery or evidence of local habitat equivalence. Filled station points contributed to supported years; hollow points are excluded or candidate assignments. [Land context](../data/external/map_context/README.md) is generalized. Reviewers can download their site comments from the page; it does not send them or change an analysis.

## Where the work stands

The required-source check passes. Original Pelican measurements are accepted, the recoverable PISCO count is repaired, and researcher exclusions apply to size donors as well as counts. Four unused 2023 KFM size groups remain unavailable. The integrated review has completed the biomass-support assessment and fitted the agreed models and uncertainty sensitivities. Missing measurements and unestimable uncertainty remain explicitly identified; adopting these results in the paper is a separate step.

The corrected source build now runs in an isolated folder. One coordinated review pipeline connects source repairs, response-specific exclusions, survey units, density, supported biomass, model eligibility and joint uncertainty. The November 5 full rebuild passed all 18 integration checks and reproduced the prior numerical tables. Its [model results](../outputs/integrated_review_20261103/index.html) and [whole-dataset review](../outputs/whole_dataset_biomass_update_20261105/index.html) include matched-year pathway and conditional heatwave refits. The manuscript, supplement and frozen paper results remain unchanged. Paper promotion still requires review of coverage, uncertainty limits, diagnostics and the downstream estimands, followed by final figure and prose concordance.

## Which workflow is being described?

| Workflow | Status | What it can establish |
|---|---|---|
| Delivered files and corrected source layer | Implemented | Trace observations, physical identities, units, corrections and exclusions. |
| Inherited biomass production | Corrected source reconstruction completed in isolation | Reproduce the existing pooled-size calculation. Its broad donor pools are not automatically validated by source repairs. |
| Integrated source-to-model review | One reproducible command; review only | Rebuild raw survey units, coverage, biomass and joint temporal models. |
| Six-history PISCO borrowing and uncertainty analysis | Incorporated within its reviewed target/calibration scope | Compare stricter visit/site matching and shared-donor uncertainty within its stated scope. |
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
| Density and biomass | Preserve measured counts; distinguish direct sizes, group means, borrowed sizes, imputation and unsupported biomass. Normalize by the surveyed area and effort. | [Integrated visit support](../outputs/dataset_inventory_20261101/integrated_visit_biomass_support.csv); data repository biomass traces and donor pools |
| Temporal analysis | Separate before–after changes from post-protection trajectories. Assess model eligibility independently of site inclusion. | [Current rules](../outputs/dataset_inventory_20261101/current_rules.csv), [freeze checklist](../outputs/dataset_inventory_20261101/freeze_readiness.csv) |
| Release | Refit all downstream analyses, verify agreement and show changes before updating the paper. | [Remaining questions](../outputs/dataset_inventory_20261101/unresolved_issues.csv) |

Use the [current analysis page](../outputs/dataset_inventory_20261101/index.html) when writing methods or results. Start by choosing a protected area or searching for a site name or code. The short area summary shows assigned protected and reference stations, kept and excluded comparisons, and recorded questions. The site card gives source-year coverage; the timeline gives usable years for each species and measurement. The decision flow and detailed tables are available below the site views. The [site table](../outputs/methods_results_reference_20261104/site_decisions.csv) and [response table](../outputs/methods_results_reference_20261104/response_histories.csv) hold the same decisions as downloadable records. Use these to verify a site-specific sentence. The [technical build guide](../../sbc-donham-kelp-mpa-data/docs/RAW_TO_BIOMASS.md) gives commands, files, formulas and safeguards. Blank timeline years mean no row in the current roster, not proof that no survey occurred.

## Apply inclusion decisions in this order

- **Ecological inclusion:** preserve researcher-advised and ecological-scope exclusions, including their program and response scope. Gull Island PISCO is distinct from excluded Gull Island KFM; Anacapa SMR is distinct from Anacapa SMCA. Do not exclude a comparison because its baseline slope, effect direction or significance is inconvenient.
- **Reference membership:** use the documented physical stations and roles. An alias is not a new station. Established assignments with incomplete habitat documentation remain provisional and receive influence checks. Missing documentation does not establish a poor match, and an established assignment does not establish equivalent habitat.
- **Coverage:** use shared comparable zones. Use documented stations available each year; constant-zone and constant-station fits are sensitivities. Actual depths must not be inferred from identical zone labels across programs.
- **Measurement support:** a valid count can support density even when sizes cannot support biomass. A missing size removes support for that response and sampling unit; it is not automatically a reason to exclude the whole site or other taxa.
- **Model eligibility:** use actual calendar years, with the designation year after protection. Three matched before years permit the flat-baseline level/recovery model; four is a sensitivity. Allowing a pre-protection trend is a second sensitivity. Other comparisons contribute post-protection trends only, with at least three observed post years for a slope and residual uncertainty. Use all supported observed years within the registered source scope; there is no nine-year cap or requirement for consecutive surveys. Three post-protection observations do not imply a long recovery history: first/last years, gaps and the observed span remain visible. These rules do not certify pBACIPS assumptions or causal attribution.

The [inclusion-gate table](../outputs/dataset_inventory_20261101/inclusion_gates.csv) records the unit, rule and consequence at each stage. Keep a physical-site exclusion separate from a missing visit, an unsupported response and an ineligible model. Apply ecological exclusions to observations, donors and uncertainty-calibration records within their documented response scope.

The 206 rows in the comparison table are candidate program × protected–reference comparison × species × response histories across 25 protected areas, not independent sites or fitted effects. Thirty-one fail researcher/ecological scope rules. Of the 175 that pass scope, 62 qualify for a before–after model, 100 for a post-protection trend only, and 13 lack enough supported post years for a fit. The `exclude` and `reason` fields preserve scope decisions; `in_current_manuscript` describes the existing paper, not approval of a revised primary analysis. Rule IDs, measurement support and model eligibility must also be read before using a history.

## Timing, depth and survey effort

The integrated review retains valid survey dates across months, averages dated visits equally within zones, zones equally within stations, and available documented stations equally within protection roles. PISCO retains CRANE and requires BOT plus MID or CNMD for the combined fish response; incomplete passes are omitted, not zero-filled. LTER quadrats are averaged within transects before equal-transect averaging. KFM areas come from each source record, preserving protocol changes.

Shared program-zone labels define the PISCO comparison footprint; physical depth and habitat equivalence are not inferred across programs. There is no protected/reference survey-date cutoff. Dates, gaps, available zones and lost coverage remain visible. KFM urchin and kelp sizes may match the nearest survey within seven days at the same physical site and year; exact-date matching is a sensitivity. Larger gaps remain unsupported. The KFM sheephead secondary review retains its own documented seven-day rule and complete-visit coverage.

The configuration files [survey and model policy](../config/integrated_review_policy.json) and [candidate assignments](../config/integrated_comparison_scope.csv) are fixed before fitting. The [annual roster](../outputs/integrated_review_20261103/comparison_year_roster.csv) and [model eligibility](../outputs/integrated_review_20261103/model_eligibility.csv) record their consequences.

## What a biomass estimate means

| Size information | Treatment | What remains uncertain |
|---|---|---|
| Individual measurements or size-frequency counts | Convert each measured size to mass; retain frequency weights and original row IDs. Divide totals by the appropriate effort and area. | Sampling, counts, conversion coefficients and dependence |
| LTER survey-group mean diameter | Retain the recorded group mean and count; label the nonlinear mass conversion as a proxy. | A mean diameter cannot recover the mean of individual masses or within-group size variation |
| Sizes from a different visit/site/program | Record the recipient and every eligible donor, matching rule and failed alternative. | Representativeness and shared donor/prediction errors |
| Five audited LTER lobster records lacking length | Use same-site/date mean individual mass as the documented imputation; preserve missing original lengths and donor identities. | The imputed mass is not a measured length |
| Positive count without supported size | Keep the count and mark biomass unavailable. | Do not manufacture zero biomass or average around missing stations without an explicit policy |
| Explicit observed zero count | Retain a zero point estimate for the surveyed unit. | Zero observations do not prove absence or zero uncertainty. An unsurveyed unit or explicit unknown remains missing. The separate restricted KFM sheephead review requires size support for its modeled zero-count uncertainty; the integrated observed-zero biomass point estimate itself does not require an invented length. |

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

## Current size matching and borrowing

**The historical production pool** selects sizes by comparison, protection role, year and species, with a 25 mm minimum for PISCO urchin donors. It can cross physical stations, dates, zones and programs. The isolated R reconstruction preserves that behavior for comparison. The integrated review uses supported visit-level measurements and declared donor rules instead; source repair alone does not validate historical pooling.

**The integrated PISCO calculation** incorporates the six-history review within its approved target and calibration scope and uses same-visit sizes first. Other-visit and other-site attempts are recorded separately. Borrowing stays within species, MPA, role, zone and year; cross-site donors and calibration exclude the target site. Actual accepted borrowing is reference-only; protected-site borrowing has not been validated. The lossless attempt table retains rejected alternatives as well as selected donors.

The integrated PISCO borrowing rule requires a donor–count date gap no larger than the largest gap observed in that comparison’s calibration data, with zero extra days of tolerance and at least two calibration site-year groups. The maximum varies by comparison and reaches 311 days among eligible donor attempts; it is not the separate KFM seven-day matching rule. Two groups is a computational minimum, not proof of adequate biological validation. The integrated review uses the strict calibrated-gap rule within the previously reviewed PISCO targets and calibration samples. The earlier six-history audit—not the current integrated roster—found 44 of 63 post-protection years under the strict rule and 45 with a seven-day margin; those historical counts do not replace the current integrated roster.

[Current visit decisions](../outputs/dataset_inventory_20261101/integrated_visit_biomass_support.csv) distinguish same-visit measurements, zero counts, borrowed visits, borrowed sites and unsupported values. [Current donor attempts](../outputs/dataset_inventory_20261101/integrated_donor_attempts.csv) record tested alternatives, eligibility and rejection reasons. The visit ledger links selected donor pools; the [size-pool inventory](../outputs/integrated_review_20261103/size_pool_inventory.csv) gives their source, date and sample support. Require both protection roles and the retained reference stations; do not silently drop unsupported stations to obtain an estimate.

The integrated model propagates shared count units, size records, imputed components and calibrated donor predictions through the nonlinear biomass calculation. Joint models estimate calendar-time dependence and shared-reference process variation; verified satellite pixel overlap informs an explicitly stated working kernel. Pooled summaries keep estimands separate, include source fixed effects and MPA heterogeneity, and average represented program means equally. Constant coverage, four-before-year, sparse-size and exact-date sensitivities remain separate. These are conditional intervals: group-mean size uncertainty, detection, coefficient covariance and incompletely documented spatial dependence are not invented. See the [uncertainty limits](../outputs/integrated_review_20261103/uncertainty_limits.csv) and [diagnostics](../outputs/integrated_review_20261103/joint_model_diagnostics.csv).

### Current biomass decision and evidence

Use a size sample from the count visit when one is available. For the six reviewed PISCO histories, allow the documented, calibrated reference-only borrowing rule when that direct sample is missing. Apply the same rule before looking at the estimated effect. Keep a valid density history when biomass lacks adequate sizes, and do not treat a missing biomass value as a zero or a site exclusion.

The [site and response table](../outputs/biomass_policy_review_20261103/comparison_support.csv) records 1,449 supported biomass comparison-years among 1,990 candidates. Twenty supported PISCO years depend on borrowed sizes. The [direct-size-only refit](../outputs/biomass_policy_review_20261103/README.md) removes those 20 years; every retained annual estimate and measurement draw is unchanged. Pooled red-urchin biomass changes from −13.19% to −18.38% per year as the fitted comparisons fall from 11 to nine. That difference is due to coverage and composition; it does not justify selecting a borrowing rule by whether an interval crosses zero. Show both cohorts when presenting a biomass result that depends on borrowing.

The [matched comparison with Emily's archived size pool](../outputs/emily_biomass_bridge_20261103/README.md) holds 30 supported PISCO urchin comparison-years fixed. In four joint-model fits, the pooled purple-urchin biomass slope moves from −6.92% to −11.84% per year when the size pool changes and current measurement covariance is held fixed. Changing only the size-driven covariance gives −6.50%. This separates the main source of the shift within that selected cohort; it does not reproduce Emily's original 1,000 draws or extend the comparison to unsupported years.

LTER records supply urchin group mean diameters, not individual sizes. The [20% and 30% within-group variation scenarios](../outputs/biomass_policy_review_20261103/lter_group_mean_scenarios.csv) barely change pooled trends when the same variation is assumed on protected and reference sides. Those percentages are illustrative assumptions, not measured LTER variation or confidence bounds. The [conversion-coefficient ledger](../outputs/biomass_policy_review_20261103/conversion_uncertainty_status.csv) records that the current files lack the coefficient covariance needed for probabilistic allometry uncertainty. Report fitted intervals as conditional on the recorded size data, conversion equations and stated covariance model. These review results have not been adopted in the manuscript.

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
| KFM excluded size donors | Enforce physical-site exclusions through aliases; recent and historical reconciliations completed. Integrated visit-level support and recovery refits are complete; downstream paper analyses remain on hold. | [KFM reconciliation](../outputs/kfm_scott_reconciliation_20261102/README.md) |
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
| [Integrated visit support](../outputs/dataset_inventory_20261101/integrated_visit_biomass_support.csv) | Retained biomass visit | Direct, local, short-gap, calibrated or unsupported sizes |
| [Integrated model eligibility](../outputs/dataset_inventory_20261101/integrated_model_eligibility.csv) | Retained comparison | Supported calendar years and eligible models |
| [Integrated donor attempts](../outputs/dataset_inventory_20261101/integrated_donor_attempts.csv) | Tested donor alternative | Eligibility, rejection reasons and calibration support; link to selected pools through the integrated visit table |
| [Current rules](../outputs/dataset_inventory_20261101/current_rules.csv) | A decision with scope and implementation status | What is approved, tested or pending |
| [Questions](../outputs/dataset_inventory_20261101/unresolved_issues.csv) and [freeze checklist](../outputs/dataset_inventory_20261101/freeze_readiness.csv) | An unresolved issue or release requirement | What still prevents adoption or needs disclosure |

A table row is not automatically an independent sample. A group label is not a physical station; a dataset hash verifies bytes, not ecological suitability. Unknown habitat, depth, coordinates or exact dates stay unknown.

Measurement errors from independent survey years have zero cross-year covariance in the conditional sampling model; shared sizes, calibration errors and LTER permanent-transect resampling retain their documented links. Biological dependence through time is estimated separately using calendar gaps. Variance optimization checks multiple starting values and the temporal-correlation profile; zero-variance and nonidentifiable components are flagged.

## KFM reference tests completed; ecological decision pending

The [reference-site sensitivity report](../outputs/kfm_reference_sensitivity_20261103/index.html) tests Scott’s flagged stations one at a time and together. The [coverage table](../outputs/dataset_inventory_20261101/kfm_reference_coverage_changes.csv), [local sheephead trends](../outputs/dataset_inventory_20261101/kfm_reference_sheephead_trends.csv) and [pooled decomposition](../outputs/dataset_inventory_20261101/kfm_reference_pooled_decomposition.csv) separate lost coverage from changes in reference definition. These are sensitivity results; the current roster remains unchanged.

Wyckoff is the sole reference for six Harris Point responses. For sheephead, comparing the same 15 post years gives −16.5% per year with the current references and −13.8% with Miracle Mile alone. Removing both references removes all Harris Point histories. Removing Webster’s Arch affects SBI sheephead density only: +3.4% becomes −0.4% per year on the same 20 years, with both intervals spanning zero. Rodes has no current measurement uses. Renormalizing the remaining station weights is explicit; retained within-station observations and ratio offsets stay fixed.

Combined omission removes seven fitted histories and 217 supported comparison-years. Pooled kelp and purple-urchin post trends change little. The flat-baseline before–after lobster estimate grows from +229% (95% interval −40 to +1695) to +424% (+77 to +1451), entirely through removing Harris Point; the pre-trend sensitivity remains very uncertain. Combined SBI sheephead uncertainty is dependence-sensitive and reaches the AR upper bound. Do not select references from these results. Scott’s ecological recommendation and Webster’s fish scope remain open. Adrian sent Scott the site-code key and these questions on 2026-11-05; no reply as of 2026-11-09. Email drafts are kept outside this public repository.

## Earlier audits and remaining release work

The [current site-limitations review](../outputs/site_limitations_current/README.md) separates work possible with existing files, contributor questions and limitations to disclose. Its [site-by-program table](../outputs/site_limitations_current/site_program_actions.csv) covers all 25 protected areas and 46 area–program groups.

Earlier [visit decisions](../outputs/dataset_inventory_20261101/measurement_visit_support_detailed.csv) and [donor attempts](../outputs/dataset_inventory_20261101/donor_attempts_current.csv) are retained as historical evidence for the six-history audit. Their counts and status statements do not describe the integrated run. Use the integrated tables above for current support and eligibility.

The [notes and correspondence audit](../outputs/notes_email_audit_20261102/README.md) and its [finding ledger](../outputs/notes_email_audit_20261102/findings.csv) preserve contributor explanations. They support decisions; old status statements do not override the current rules above.

SBC50.21 is also a pinned review input, not a replacement production source. The newer LTER adult kelp source contains the same observations through 2023 as our registered input. SBC50.21 combines adult and juvenile kelp; its reported density must not be treated as adult stipe density. Its adult biomass values imply a fixed coefficient, whereas our pipeline uses the observed survey month. Keep our adult foliar response and actual-month conversion separate from that provider aggregate. The linked provider protocol and implied coefficient need clarification before substituting the aggregate. New 2024–2026 observations do not enter the current cohort automatically. These are source and formula comparisons, not new fitted protection effects.

Emily clarified that her historical minimum was three before observations at the same physical stations, assessed separately for fish and benthic surveys, and that she counted the implementation year as before. A stricter calendar-year convention is a deliberate change, not a recovered original rule. The integrated review now freezes that revised convention; its matched station-year counts and eligibility remain visible. Her old `Sites_Years.xlsx` conveys coverage by cell colors; blank cell values must not be interpreted as absent surveys.

The coordinated source-to-recovery review has run with frozen coverage and model rules. Review its donor, missingness, model diagnostics and sensitivity tables before adopting it. [The November 5 downstream review](../outputs/whole_dataset_biomass_update_20261105/README.md) pairs pathway responses on the same years, stations and zones and refits heatwave periods with updated annual ratios. Those conditional models do not automatically replace the manuscript's Table 3 or full resilience suite. Final paper figures and manuscript concordance remain to be rebuilt after deciding which estimands to report. Ratio-offset sensitivity, satellite quarterly timing and unresolved habitat/permission limits remain explicit in the inventory. Resolve permissions before redistribution; show results before altering the paper.

