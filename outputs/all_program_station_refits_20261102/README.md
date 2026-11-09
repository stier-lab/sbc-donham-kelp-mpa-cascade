# Station coverage and reference checks across programs

November 2, 2026. Review only. Researcher exclusions and the approved 175-history review cohort are preserved; no manuscript, supplement, canonical results or source data were changed.

## Scope and reconstruction

Rebuilt the existing station-level PISCO/KFM/LTER responses in an isolated copy of the processing pipeline, retaining cached biomass calculations and calibrated coefficients. Preserved the approved complete-pass/repeat/CRANE sheephead reconstruction and the separately reconstructed Point Dume membership. Recovered satellite annual values at the delivered polygon level after the documented empty-reference removals. All175 histories reproduce the approved annual ratios to <1e-7 with a single nonnegative additive correction per history. Corrections were recovered algebraically from the existing ratios, held fixed and checked over every year; none was tuned to a desired sensitivity result. This validates the aggregation/reconstruction, not every historical survey convention or size donor.

There are35 KFM,18 LTER,17 satellite and105 PISCO histories. Histories are response/program/MPA combinations, not independent sites. Six PISCO histories lack three observed post-protection years;169 recovery slopes are reproducible. Current counts/size support and inherited seasonal/protocol filters are retained so the comparisons isolate reference handling, rather than mixing new measurement policies.

## What changed

Reference membership changes during the post-protection record in58 PISCO,7 KFM,1 satellite and0 LTER histories. Across the entire before-and-after record the corresponding counts are62,7,3 and2. Requiring a reference to be present before protection as well as after can therefore discard a perfectly continuous recovery-period reference. Both full-record and post-only persistent-reference alternatives are supplied, separately labeled.

PISCO Scorpion kelp is the largest newly identified persistent-reference sensitivity: the annual protected/reference slope changes0.2357→0.0044 on the same20 post-protection years. The intervals are[0.1093,0.3620] and[−0.0790,0.0877]. However, adding station intercepts on the separate log-station scale changes0.1829→0.1824. The narrow reference subset represents different reference locations; its result alone does not establish that changing station membership caused the original trend. The separate KFM Scorpion kelp slope remains0.3884 under these reference checks.

At Campus Point, LTER sheephead reference membership is stable during recovery. Requiring full-record continuity removes Goleta Beach because of baseline coverage and changes biomass slope−0.4624→−0.2441. Requiring only recovery-period continuity leaves−0.4624 unchanged. This is a reference-population choice, not evidence of turnover during recovery.

At Santa Barbara Island, KFM lobster density changes0.1284→0.2049 when restricted to the full-record persistent reference set. Omitting Arch Point in a separate same-year comparison changes0.1392→0.0018 across18 post-years. Those previous slopes differ because the omission comparison is restricted to its own common years. Do not conflate year changes and reference effects.

## Reference identity, even with stable coverage

Removing Mohawk from Campus Point LTER reverses the red-urchin biomass point estimate from−0.1557 to+0.1302 on the same11 years; density changes−0.1235→+0.0117. At Santa Barbara Island KFM, removing Arch Point changes purple-urchin biomass−0.0681→+0.1053 on17 years and density−0.0472→+0.1259 on18 years. These signs describe exploratory point estimates, not new declarations of significance.

Satellite Naples is sensitive to which of its two reference polygons represents outside canopy. Omitting reference polygon1 changes the slope+0.0092→+0.4442 on the same9 post-years. Swamis also differs substantially between references. Polygon totals describe different spatial footprints; these are existing delivered references, not interchangeable diver stations or independent replicates. Empty masks remain unavailable, not zero canopy. No polygon boundary or annual canopy reconstruction is changed.

The source/response table gives every check, including stable cases. No influential reference has been removed and no effect direction determines eligibility.

## Interpretation and recommendation

Keep the researcher-supported assignments. Carry station/polygon identity through analysis and distinguish changing membership from differences in the locations represented. A permanent-reference subset is a useful sensitivity, not automatically a better primary comparison. Existing independent programs can disagree because they sample different footprints; do not merge their station identities or force agreement.

For the main analysis, these results support explicit reference dependence and transparent influence checks. They do not yet select a new primary station model. Log-station models average transformed station values and answer a different question from log ratios of arithmetic means; only their adjusted-minus-unadjusted contrast isolates the station-intercept change. Constant intercepts cannot remove changing habitat, season or station-specific biological trends.

## Biomass, missingness and other limits

Reproducing biomass does not certify sparse or shared size samples. The expanded count-to-size linkage/uncertainty inventory is still required; no new biomass was fabricated for unsupported KFM predators. This check retains existing measurement processing, including unresolved swath visit weighting and source missingness questions, rather than silently fixing them along with reference selection. Different sample sizes for biomass and density remain explicit.

52 histories have only one reference unit, so a reference-omission comparison is unavailable.24 lack any reference present throughout their entire approved record; these remain in the available-reference analysis, not excluded. The post-only persistent subset has broader coverage. No depth/timing rule, before–after estimator or primary cohort was changed. These fits do not establish whether protection caused the observed differences.

## Validation

All169 existing recovery slopes and calendar-gap HAC3 SEs reproduce to <5e-14. R independently reproduces338 log-station coefficients to <1e-12 and447 reference-omission coefficients and interval endpoints to <1e-12. There are950 diagnostic fits plus447 paired omission checks. All487 protected manuscript/result files remain unchanged. The log-station models use equal role/year weights, station/year intercepts and year-summed HAC3 scores; their small-sample intervals remain approximate. No new regional pooled, pathway, resilience or before–after results are claimed.

## Reproduce

Run stage_all_program_station_refits.py; run Rscript run.R inside its ignored stage; run prepare_all_program_stations.py, refit_all_program_stations.py, verify_all_program_stations.R, verify_all_program_omissions.R and report_all_program_stations.py from the analysis root. Source hashes, annual-reconstruction gates, exact model rows, role values, coverage and omissions are saved beside this report. Staging contains source symlinks plus private copies of caches; original source files are not edited. No emails, commit or push.
