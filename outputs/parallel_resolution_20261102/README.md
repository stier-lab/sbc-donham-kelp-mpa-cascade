# Four-area audit: findings and recommended next steps

November 2, 2026. Four agents independently checked source files, code and dated audit outputs. The coordinator reproduced the major counts and source checks below. Recommendations are not adopted inclusion rules. The primary analysis and manuscript remain unchanged.

## 1. KFM donor identities

The hold is justified. Four KFM and two PISCO biomass histories use the affected KFM sizes. The local audit found no additional original post-2011 measurements that safely resolve the aliases. Historical checks reproduce 65 size-distribution identities, four 2015 grouped checks and 14 physical count-site comparisons.

The six histories contain 100 post-protection comparison-years: 45 already have historical diagnostic replacements, 54 still require exact original sizes, and one is a verified zero-reference year needing no size replacement. These are overlapping comparison-years, not 54 independent missing samples. Binned 2015 summaries confirm identities but cannot reconstruct nonlinear biomass from individual sizes.

**Recommendation:** keep density, keep these biomass histories provisional, and use the already-sent Scott request for original site-specific sizes and dates. Prioritize permitted sites 1 and 14, with sites 21/37 or construction code to resolve the combined labels. No duplicate email is needed. This audit checked cached correspondence, not live inbox freshness.

[Detailed audit](kfm/report.md) · [Exact year table](kfm/history_years.csv)

## 2. Unsupported PISCO biomass

Current coverage is correctly reported as **44 of 63 post-protection comparison-years**. The 38 unsupported species-zone-date uses span 24 site-date combinations; 26 post-protection uses block 19 comparison-years. Across all years, 20 uses lack eligible donors, 17 fail calibration requirements and one fails the date-gap guard. No missed eligible same-visit measurements or explanatory bad join was found in the examined delivery.

One case deserves a less brittle rule: Point Fermin red urchins in 2019 have local sizes 294 days after counts, four days beyond the observed calibration maximum of 290 days. A bounded diagnostic allowing this excess restores one paired year, giving **45/63**, but does not validate the uncertainty. The calibration uses only two site-year groups. Four other local candidates are only one or two days apart; all are before protection and restore no post years.

Point Vicente protected-site borrowing fails partly because nested validation leaves too few protected stations, not because every donor is absent. Five post-protection comparison-years are blocked solely by that calibration problem. Reference-site residuals should not automatically stand in for protected-site validation.

**Recommendation:** retain density; distinguish missing measurements from unavailable calibration; review a small, consistently applied date tolerance with exact gaps exposed. Keep the 45/63 result diagnostic until the rule and uncertainty treatment are reviewed. Preserve both failed borrowing attempts in the ledger so “no donors” does not conceal available local sizes that failed calibration.

[Detailed audit](support/report.md) · [Cell-by-cell causes](support/unsupported_use_trace.csv) · [Blocked years](support/blocked_comparison_years.csv)

## 3. Biomass uncertainty

Arithmetic checks passed. They do not establish total uncertainty. Eight available temporal validation predictions represent two stations in one year, for red urchins and lobster; none validates purple-urchin temporal borrowing. Of 624 measured visit uses, 84 contain one animal; the empirical resampling cannot estimate population size variation from those singletons.

Shared-donor, persistent station and calibration-distribution uncertainty remain incomplete. The earlier joint count/size audit and latest borrowing audit use different observation policies and must not be combined by simply adding their variances. Existing code does not make that mechanical error. KFM kelp visit variation is also largely unidentifiable: 216 of 220 station-years with plant measurements have just one survey date.

**Recommendation:** first build one observation-identity registry on the current survey roster, then couple reused count units, sizes and references in one uncertainty calculation. Separate sample noise from size-transfer error and combine observation and temporal-process uncertainty once. Keep temporal borrowing exploratory while its calibration is this sparse. Use coefficient uncertainty only when source covariance/calibration data support it; otherwise state the conditional limit. Do not invent singleton variance or call fixed-offset sensitivity a measured random uncertainty.

[Detailed audit and staged implementation](uncertainty/report.md) · [Calibration replication](uncertainty/calibration_groups.csv) · [Validation replication](uncertainty/validation_groups.csv)

## 4. Depth, timing and reference consistency

The approved all-month/shared INNER–MID treatment remains scoped to six PISCO animal comparisons and their density/biomass responses. Comparable zone labels do not establish identical physical depth across programs. KFM/LTER depth gaps are metadata limitations, not reasons by themselves to discard valid density records. Satellite comparisons require spatial-footprint checks rather than diver depth rules.

For the 105 retained PISCO response histories in the station audit, requiring a persistent reference through the whole post-protection record leaves 17 histories without a common reference, including two already lacking post observations, and removes 187 of 1,152 post-protection comparison-years. Requiring persistence through before and after years leaves 23 histories without an intersection and removes 269 comparison-years. These counts include density/biomass counterparts and are not independent sites. They describe the station audit's cohort, before the later six-comparison biomass restrictions.

**Additional baseline correction:** the six-comparison density rerun inherited the biomass station/year roster. It omits 33 predesignation station-year uses, representing six density comparison-years across four histories: Abalone lobster 2011; Point Vicente red and purple 2008; Point Vicente lobster 2008–2010. No post-protection density years are lost. Rebuild from response-specific density records before any before–after analysis; the shared-zone claim currently applies to the restricted roster.

**Recommendation:** keep researcher-approved reference assignments and use fixed-reference and station-adjusted comparisons to reveal sensitivity. Decide any restriction from survey design and coverage, not its effect on the slope. Use shared comparable zones within each program; do not impose a universal numeric-depth or date-gap cutoff. Keep all-month expansion outside the six comparisons a separate decision.

[Detailed coverage audit](coverage/report.md) · [Reference coverage by program](coverage/reference_by_program.csv)

## Documentation corrected now; scientific choices still pending

Verified factual corrections are included in the current guide generator: R04 now links to all twelve affected response histories, R05 to the 34 retained PISCO sheephead histories, and diver depth rule R03 no longer attaches to satellite rows. Count units for unsupported biomass are clarified. No inclusion, missing-data treatment, uncertainty model or scientific estimate changed.

The next table update should add separate fields for donor identity, exact unresolved years, local versus cross-site failure reasons, calibration stations/years, uncertainty components and recommended action. A proposed rule must remain distinct from an implemented rule. The linked `recommended_actions.csv` records those dispositions.

## Reproduction and limits

Each area includes an independent check script, report and evidence tables. Scripts identify local input paths; some retain their original scratch output directory, while others write beside the script. Inspect the output path before rerunning. `agent_evidence_manifest.json` hashes the collected evidence; `audit_manifest.json` records the initial guide hashes. `coordinator_checks.json` records independent coordinator checks. No raw email bodies or bulk source deliveries were copied.

Use `python3 code/tools/collect_parallel_resolution.py` to collect reports from the isolated scratch directories; it does not adopt recommendations. Existing primary/manuscript files remain protected by `outputs/researcher_scope_20261101/protected_hashes.json`. The consolidated guide and its renderer continue to be the current policy entry point.
