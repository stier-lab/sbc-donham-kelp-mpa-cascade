# Unsupported PISCO biomass support audit — 2026-11-02

Read-only audit under parent session `393ab715dc104db49d989bacdd9f2f3e`. All scripts and diagnostics are confined to `/tmp/donham_swarm_20261102/support`. No repository, guide, primary result, source data, manuscript, status, or continuity edits; no messages, pushes, extra agents, or full model runs.

## Verified result and count units

**44/63 is accurate for the current strict borrowing implementation.** The 38 unsupported positive entries are comparison/species/zone/site/date uses, spanning **24 distinct physical site-date combinations**, 29 site-zone-date combinations, and 34 species-site-date combinations. These are descriptive keys, not verified independent campaigns or replicates. **26 post-protection uses block 19 comparison-years**; the other 12 uses block seven non-post comparison-years and do not reduce the 63-year post denominator. Across all years, 26 uses are in the protected role and 12 in the reference role (a different partition from post/non-post).

The independent raw reconstruction reproduced 682 size sample-use inventory rows, their full keys, frequency-weighted masses, all 38 unsupported positive count densities, all annual completeness flags, and equality of the guide's 1,086 visit rows with the scientific ledger. The 20 missing-data uses, 17 calibration-guard uses, and one date-gap-guard use sum to 38. Among post uses these counts are 19, six, and one. Counts of affected years must be deduplicated: one purple2016 comparison-year has both missing-data and calibration failures.

No valid same-visit recovery was found in the delivered combined size table or direct swath size field. No bad join was found that explains these 38 uses. Rebuilding the inventory from raw records finds no omitted eligible size sample, and the `site_new` alternative key offers no rescue. Adult urchin sizes below25mm remain ineligible: LONG_POINT_E purple MID2022 has two such raw rows, not an adult size measurement usable under the retained protocol.

## Point Fermin: four days beyond an empirical boundary

Abalone red, POINT_FERMIN reference, INNER, count date2019-01-04 (row80), has four valid local animals sized2019-10-25. The local gap is294days; calibration support ends at290days. Local temporal calibration comprises just **two site-year groups**, represented by four reciprocal prediction rows: 3_PALMS_EAST2015 at283days and BUNKER_POINT2015 at290days. Both calibration groups are MID2015; the prediction donor itself stays at the same site/species/zone/year as required. Current calibration residuals pool zones and years within series/role; donor restrictions and calibration pooling are distinct and should be disclosed.

The alternative cross-site candidate BUNKER_POINT is224days away, versus its cross-site calibration maximum150days; that is a different, much larger support extrapolation. The actionable near-boundary case is the local294vs290day match.

**Recommendation, not an adopted rule:** honor the user's direction that a few days should not cause exclusion by replacing the exact observed-maximum cliff with an explicitly reviewed small tolerance or coarse rounded boundary sensitivity. A tolerance of at least four days admits this case; a seven-day margin is an illustrative diagnostic, not an established biological threshold. Keep the actual gap, excess, sample count, and calibration group count visible. The four-day excess is not a strong reason by itself to discard the local source, but the roughly ten-month transfer and two-group calibration remain weak evidence for predictive uncertainty. A tolerant boundary does not manufacture adequate temporal validation, and four reciprocal rows are not four independent groups.

A bounded scratch-only diagnostic fills row80 from its existing local donor and current median calibration correction, leaving all other sites/years/exclusions and biomass values fixed. Every required cell in Abalone red2019 then becomes finite: **Abalone red7/11→8/11; total44/63→45/63**. This is exactly one restored paired comparison-year, with 18 post years still blocked. It neither certifies intervals nor estimates a new primary effect. See `point_fermin_diagnostic.json` and `point_fermin_calibration.csv`.

## Other recoverable local evidence and calibration limitations

Four other unsupported uses already have valid local, same-zone/year measurements on nearby dates:

- LONG_POINT_E red MID and purple INNER: count2010-07-21, sizes2010-07-20 (one day).
- ROCKY_POINT_N red and purple INNER: count2009-08-19, sizes2009-08-21 (two days).

These all fail temporal calibration because zero qualifying groups remain after the outer-site exclusion. None contributes to the post-protection63year denominator. They are legitimate local candidates for a campaign-date review or transparent diagnostic sensitivity; they are not verified same-visit matches, and their recovery restores zero post years. `local_candidate_source_rows.csv` preserves exact raw row positions and original frequencies. No records were deduplicated by size value.

At Point Vicente, 15 protected uses have other protected-site donors but zero nested cross-site calibration groups. There are only two retained protected stations: after outer target exclusion one remains, and that station has no other protected station from which to obtain its own inner prediction. This is structural lack of calibration, not lack of all donor data. Six of the15 uses are post. Five post comparison-years are blocked only by this guard (red2015/2019; purple2017/2022/2023). Purple2016 also has two no-data cells, so solving its calibration cell alone does not restore the year. Do not claim protected borrowing has been validated or borrow reference residuals automatically.

## Exact post-year impact

| Comparison | Strict coverage | Blocked post years |
|---|---:|---|
| Abalone Cove SMCA / M. franciscanus / PISCO / Bio | 7/11 | 2014;2019;2022;2023 |
| Abalone Cove SMCA / P. interruptus / PISCO / Bio | 10/11 | 2014 |
| Abalone Cove SMCA / S. purpuratus / PISCO / Bio | 7/11 | 2014;2018;2021;2022 |
| Point Vicente SMCA / M. franciscanus / PISCO / Bio | 5/10 | 2015;2016;2019;2021;2022 |
| Point Vicente SMCA / P. interruptus / PISCO / Bio | 9/10 | 2022 |
| Point Vicente SMCA / S. purpuratus / PISCO / Bio | 6/10 | 2016;2017;2022;2023 |

Full per-cell site/species/year/zone/role/date coverage and both failed attempts are in `unsupported_use_trace.csv`. `blocked_comparison_years.csv` links each distinct comparison-year to row IDs, sites, roles, and reasons. Do not sum separately tabulated reason-specific year counts without removing overlaps.

## Minimal next actions and precise guide changes needed

1. Keep current44/63 in the guide as the current implementation's coverage; add45/63 only as the explicit Point Fermin local-date diagnostic. Retain the warning that25local-only versus44borrow-supported uses different temporal-support rules.
2. Update `code/tools/propagate_pisco_borrowing.py` in a future authorized change to retain both attempted paths. Lines75–79 overwrite the same-site failure when testing other-site borrowing; unsupported rows lose candidate and calibration details. In particular ROCKY_POINT_N red/purple2009 currently say no donors although local sizes exist two days later. This is diagnostic information loss, not an incorrect same-visit join or evidence to impute zero.
3. In `outputs/dataset_inventory_20261101/measurement_visit_decisions.csv`, add separate local/other-site candidate availability and failure fields, dates, donor animals/sites, calibration groups/max gap, target gap/excess, post indicator, and blocked comparison-year identifier. Use physical/species/date keys in addition to the current build-specific row IDs. Update the generating consolidation code so these persist on rebuild.
4. Amend the guide `README.md`, `unresolved_issues.csv` U08, and support-summary prose generated for `measurement_support_current.csv` to say38uses across24site-dates;26post uses block19years. Split no-data versus calibration versus date-gap failures rather than calling every case missing measurements or protected-site gaps. The guide README already correctly calls the1,086 entries uses; U08 and the borrowing README still call the38 entries visits.
5. `current_rules.csv` R11 should distinguish the implemented exact empirical-maximum guard from the proposed tolerant boundary and explicitly disclose sparse calibration. R10 should distinguish same-zone/year donors from calibration residuals pooled across zones/years. No new rule is adopted by this audit.
6. For remaining missing-data cases, use the row-key ledger for a narrowly scoped original-size-record request if later authorized. Do not relax researcher exclusions, cross species/MPA/role/zone/year for donors, or treat unavailable data as zero. No external request is needed to complete this audit.

## Unresolved evidence and bounds

The combined delivered size file is the operational raw input; it does not prove whether absent measurements were never collected or omitted from a delivery. This audit did not obtain original field sheets or search remote holdings. The four one-/two-day matches need campaign identity evidence before being called same-visit. Two-group temporal calibration is too sparse to establish predictive accuracy, and count, allometric, seasonal and persistent dependence uncertainties remain outside this bounded check. No new slope fit, uncertainty claim, primary promotion, or ecological exclusion is proposed from these results.

## Reproduce and files examined

Run `python3 /tmp/donham_swarm_20261102/support/check_support.py` (pandas/numpy). It writes only this scratch directory. Then run `write_report.py` to regenerate this report and `findings.csv`. All assertions passed in the final bounded run. A one-off exploratory unfiltered raw-date parse initially encountered an invalid date outside the target check; the reproducible audit uses coercion plus exact target keys, and all38 target count checks pass. No target dates were repaired or invented.

`examined_data_files.json` lists exact data paths, resolved symlink targets, and SHA256 hashes. Substantive data examined:
- `/Users/adrianstier/sbc-donham-kelp-mpa-cascade/outputs/pisco_biomass_borrowing_20261102/visit_biomass_support.csv`
- `/Users/adrianstier/sbc-donham-kelp-mpa-cascade/outputs/pisco_timing_donors_20261102/size_sample_inventory.csv`
- `/Users/adrianstier/sbc-donham-kelp-mpa-data/data/ALL_sizefreq_2024.csv`
- `/Users/adrianstier/sbc-donham-kelp-mpa-cascade/outputs/pisco_biomass_borrowing_20261102/outer_site_calibration.csv`
- `/Users/adrianstier/sbc-donham-kelp-mpa-cascade/outputs/pisco_biomass_borrowing_20261102/biomass_annual.csv`
- `/Users/adrianstier/sbc-donham-kelp-mpa-data/data/PISCO/MLPA_kelpforest_swath_2024.csv`
- `/Users/adrianstier/sbc-donham-kelp-mpa-cascade/outputs/dataset_inventory_20261101/measurement_visit_decisions.csv`

Additional code/documentation examined (paths below relative to analysis root unless absolute):

- `AGENTS.md`, `CLAUDE.md`, `.ai-continuity/CURRENT.md`, recent `.ai-continuity/ACTIVITY.jsonl`.
- `outputs/dataset_inventory_20261101/README.md`, `current_rules.csv`, `unresolved_issues.csv`.
- `outputs/pisco_biomass_borrowing_20261102/README.md`, `input_hashes.json`, `biomass_coverage.csv`.
- `code/tools/run_pisco_biomass_review.py`, `propagate_pisco_borrowing.py`, `check_pisco_timing_donors.py`, `refit_pisco_all_months.py`, `audit_pisco_size_zones.py`; relevant report-generator matches in `report_pisco_borrowing.py`.
- `/Users/adrianstier/sbc-donham-kelp-mpa-data/DATA_SOURCES.md` source inventory and relevant import/filter sections of `code/R/03_data_import.R`; directory inventory of `data/PISCO/`. No AGENTS.md/CLAUDE.md existed at the data-companion root.
- `/Users/adrianstier/sbc-donham-kelp-mpa-manuscript/AGENTS.md`, `CLAUDE.md`, `.ai-continuity/CURRENT.md` for applicable instructions only; no manuscript content reviewed.
- `/Users/adrianstier/AI-Workspace/WORKFLOW.md`, `ACTIVITY-WORKFLOW.md`, matching project record in `projects.json`, and current matching handoff `/Users/adrianstier/AI-Workspace/handoffs/Donham-Stier-CA-MPA-2026.md`.
- `/Users/adrianstier/.codex/skills/data-provenance-audit/SKILL.md`.

## Cell-level appendix

| Row ID | Comparison | Site | Role | Year | Zone | Date | Post | Root reason |
|---:|---|---|---|---:|---|---|---|---|
| 43 | Abalone Cove SMCA / M. franciscanus / PISCO / Bio | ABALONE_COVE_KELP_W | mpa | 2022 | MID | 2022-07-06 | True | no eligible local or same-role donor data |
| 68 | Abalone Cove SMCA / M. franciscanus / PISCO / Bio | BUNKER_POINT | reference | 2023 | INNER | 2023-09-21 | True | no eligible local or same-role donor data |
| 74 | Abalone Cove SMCA / M. franciscanus / PISCO / Bio | POINT_FERMIN | reference | 2014 | INNER | 2014-07-18 | True | no eligible local or same-role donor data |
| 80 | Abalone Cove SMCA / M. franciscanus / PISCO / Bio | POINT_FERMIN | reference | 2019 | INNER | 2019-01-04 | True | date-gap guard |
| 94 | Abalone Cove SMCA / M. franciscanus / PISCO / Bio | WHITES_POINT | reference | 2014 | INNER | 2014-10-02 | True | no eligible local or same-role donor data |
| 194 | Abalone Cove SMCA / P. interruptus / PISCO / Bio | WHITES_POINT | reference | 2014 | INNER | 2014-10-02 | True | no eligible local or same-role donor data |
| 216 | Abalone Cove SMCA / S. purpuratus / PISCO / Bio | 3_PALMS_EAST | reference | 2018 | INNER | 2018-11-16 | True | no eligible local or same-role donor data |
| 217 | Abalone Cove SMCA / S. purpuratus / PISCO / Bio | 3_PALMS_EAST | reference | 2018 | MID | 2018-11-16 | True | no eligible local or same-role donor data |
| 249 | Abalone Cove SMCA / S. purpuratus / PISCO / Bio | ABALONE_COVE_KELP_W | mpa | 2021 | MID | 2021-06-09 | True | no eligible local or same-role donor data |
| 251 | Abalone Cove SMCA / S. purpuratus / PISCO / Bio | ABALONE_COVE_KELP_W | mpa | 2022 | MID | 2022-07-06 | True | no eligible local or same-role donor data |
| 257 | Abalone Cove SMCA / S. purpuratus / PISCO / Bio | BUNKER_POINT | reference | 2012 | MID | 2012-07-10 | False | no eligible local or same-role donor data |
| 282 | Abalone Cove SMCA / S. purpuratus / PISCO / Bio | POINT_FERMIN | reference | 2014 | INNER | 2014-07-18 | True | no eligible local or same-role donor data |
| 302 | Abalone Cove SMCA / S. purpuratus / PISCO / Bio | WHITES_POINT | reference | 2014 | INNER | 2014-10-02 | True | no eligible local or same-role donor data |
| 358 | Point Vicente SMCA / M. franciscanus / PISCO / Bio | LONG_POINT_E | mpa | 2010 | INNER | 2010-07-21 | False | calibration guard |
| 359 | Point Vicente SMCA / M. franciscanus / PISCO / Bio | LONG_POINT_E | mpa | 2010 | MID | 2010-07-21 | False | calibration guard |
| 368 | Point Vicente SMCA / M. franciscanus / PISCO / Bio | LONG_POINT_E | mpa | 2016 | INNER | 2016-08-12 | True | no eligible local or same-role donor data |
| 378 | Point Vicente SMCA / M. franciscanus / PISCO / Bio | LONG_POINT_E | mpa | 2021 | INNER | 2021-07-30 | True | no eligible local or same-role donor data |
| 380 | Point Vicente SMCA / M. franciscanus / PISCO / Bio | LONG_POINT_E | mpa | 2022 | INNER | 2022-09-07 | True | no eligible local or same-role donor data |
| 409 | Point Vicente SMCA / M. franciscanus / PISCO / Bio | POINT_VICENTE_W | mpa | 2009 | INNER | 2009-11-03 | False | calibration guard |
| 412 | Point Vicente SMCA / M. franciscanus / PISCO / Bio | POINT_VICENTE_W | mpa | 2011 | INNER | 2011-10-12 | False | calibration guard |
| 413 | Point Vicente SMCA / M. franciscanus / PISCO / Bio | POINT_VICENTE_W | mpa | 2011 | MID | 2011-10-12 | False | calibration guard |
| 414 | Point Vicente SMCA / M. franciscanus / PISCO / Bio | POINT_VICENTE_W | mpa | 2012 | INNER | 2012-08-10 | False | calibration guard |
| 418 | Point Vicente SMCA / M. franciscanus / PISCO / Bio | POINT_VICENTE_W | mpa | 2015 | INNER | 2015-09-23 | True | calibration guard |
| 420 | Point Vicente SMCA / M. franciscanus / PISCO / Bio | POINT_VICENTE_W | mpa | 2016 | INNER | 2016-06-22 | True | no eligible local or same-role donor data |
| 426 | Point Vicente SMCA / M. franciscanus / PISCO / Bio | POINT_VICENTE_W | mpa | 2019 | INNER | 2019-06-14 | True | calibration guard |
| 430 | Point Vicente SMCA / M. franciscanus / PISCO / Bio | POINT_VICENTE_W | mpa | 2021 | INNER | 2021-07-28 | True | no eligible local or same-role donor data |
| 507 | Point Vicente SMCA / M. franciscanus / PISCO / Bio | ROCKY_POINT_N | reference | 2009 | INNER | 2009-08-19 | False | calibration guard |
| 685 | Point Vicente SMCA / P. interruptus / PISCO / Bio | POINT_VICENTE_W | mpa | 2022 | MID | 2022-07-28 | True | no eligible local or same-role donor data |
| 863 | Point Vicente SMCA / S. purpuratus / PISCO / Bio | LONG_POINT_E | mpa | 2009 | MID | 2009-05-26 | False | calibration guard |
| 864 | Point Vicente SMCA / S. purpuratus / PISCO / Bio | LONG_POINT_E | mpa | 2010 | INNER | 2010-07-21 | False | calibration guard |
| 874 | Point Vicente SMCA / S. purpuratus / PISCO / Bio | LONG_POINT_E | mpa | 2016 | INNER | 2016-08-12 | True | no eligible local or same-role donor data |
| 876 | Point Vicente SMCA / S. purpuratus / PISCO / Bio | LONG_POINT_E | mpa | 2017 | INNER | 2017-09-08 | True | calibration guard |
| 887 | Point Vicente SMCA / S. purpuratus / PISCO / Bio | LONG_POINT_E | mpa | 2022 | MID | 2022-09-07 | True | calibration guard |
| 915 | Point Vicente SMCA / S. purpuratus / PISCO / Bio | POINT_VICENTE_W | mpa | 2009 | INNER | 2009-11-03 | False | calibration guard |
| 926 | Point Vicente SMCA / S. purpuratus / PISCO / Bio | POINT_VICENTE_W | mpa | 2016 | INNER | 2016-06-22 | True | no eligible local or same-role donor data |
| 927 | Point Vicente SMCA / S. purpuratus / PISCO / Bio | POINT_VICENTE_W | mpa | 2016 | MID | 2016-06-22 | True | calibration guard |
| 940 | Point Vicente SMCA / S. purpuratus / PISCO / Bio | POINT_VICENTE_W | mpa | 2023 | INNER | 2023-08-18 | True | calibration guard |
| 1013 | Point Vicente SMCA / S. purpuratus / PISCO / Bio | ROCKY_POINT_N | reference | 2009 | INNER | 2009-08-19 | False | calibration guard |
