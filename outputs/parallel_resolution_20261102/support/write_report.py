from pathlib import Path
import pandas as pd,json
O=Path('/tmp/donham_swarm_20261102/support');R=Path('/Users/adrianstier/sbc-donham-kelp-mpa-cascade')
x=pd.read_csv(O/'unsupported_use_trace.csv');b=pd.read_csv(O/'blocked_comparison_years.csv');s=json.loads((O/'summary.json').read_text())
f=[
('SUP01','Six PISCO animal histories; all years','summary.json; blocked_comparison_years.csv','38 unsupported species/zone/comparison uses represent 24 distinct site-date combinations (29 site-zone-dates; 34 species-site-dates), not 38 independent visits. 26 post uses block 19 comparison-years; remaining 12 uses affect seven non-post comparison-years.','Use explicit count units and link uses to blocked comparison-years.','verified'),
('SUP02','Six histories; post-protection years','check_support.py; summary.json','Independent completeness reconstruction confirms 44/63 biomass comparison-years and exactly the saved supported flags; no unsupported reference station is silently omitted.','Retain 44/63 as current strict-review coverage, with diagnostic exceptions separate.','verified'),
('SUP03','38 unsupported uses','unsupported_use_trace.csv; raw_count_checks.csv; raw_swath_size_checks.csv','20 uses lack both eligible local-year sizes and eligible same-role other-site donors; 17 have candidate data but fail calibration-group guards; one fails date-gap guards. All 38 positive densities reproduce from raw counts; none has eligible same-visit sizes in either size inventory or swath.','Separate missing source data from computational guards; do not describe all 38 as missing measurements.','verified'),
('SUP04','POINT_FERMIN; Abalone red; INNER; 2019-01-04','point_fermin_diagnostic.json; point_fermin_calibration.csv; local_candidate_source_rows.csv','Four local animals sized 2019-10-25 give gap294 vs empirical max290. Diagnostic allowance restores Abalone red2019:7/11 to8/11 and total44/63 to45/63. Only two calibration site-years/four reciprocal rows exist, both MID2015.','Recommend a small explicit rounding/tolerance sensitivity consistent with user instruction that a few days should not exclude; retain sparse temporal-calibration warning. Do not claim validation or adopt a threshold here.','verified diagnostic; recommendation not adopted'),
('SUP05','Point Vicente protected uses','unsupported_use_trace.csv; available_donor_evidence.csv','15 protected uses have same-role cross-site donors but zero nested outer-site calibration groups; two protected stations leave only one after holding target out, so cross-site calibration is structurally unavailable. Six such uses are post, and only five post years are blocked solely by this guard; purple2016 also has data gaps.','Retain donor availability independently of calibration eligibility. Prefer original missing local samples; any alternative uncertainty model requires explicit review, not removal of the guard by counting rows as replication.','verified'),
('SUP06','Five local other-date candidates','local_candidate_source_rows.csv; unsupported_use_trace.csv','Besides Point Fermin, LONG_POINT_E red MID and purple INNER2010 have one-day local matches; ROCKY_POINT_N red/purple INNER2009 have two-day local matches. All four fail temporal calibration (zero groups), and all are non-post.','Flag these short-date matches for campaign metadata review; they restore no post years in the 63-year denominator. Do not call them same-visit measurements without evidence.','verified'),
('SUP07','Unsupported ledger reasons and consolidated guide','code/tools/propagate_pisco_borrowing.py lines75-79; measurement_visit_decisions.csv','Final reason overwrites same-site attempt with cross-site failure. ROCKY_POINT_N red/purple2009 are labeled no donors despite local samples two days later. Unsupported rows also lose available donors and calibration diagnostics.','Preserve per-attempt reason, candidate dates/animals/sites, target gap, calibration maximum and group count; derive combined reason and post-year impact. Update generator before regenerating guide.','verified'),
('SUP08','Eligible size inventory and joins','examined_data_files.json; check_support.py','All682 inventory sample-use rows reconstructed from combined raw sizes and fixed retained mapping; keys and weighted means match. No alternate site_new rescue found for38 gaps. LONG_POINT_E purple MID2022 contains only sub25mm local sizes, ineligible for adult counts.','No join repair or cutoff relaxation justified; retain frequency multiplicity and adult cutoff. Scope is delivered local holdings, not unavailable original field records.','verified negative')]
pd.DataFrame(f,columns=['finding_id','scope','evidence','finding','recommended_action','status']).to_csv(O/'findings.csv',index=False)
text='''# Unsupported PISCO biomass support audit — 2026-11-02

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

'''
text+='| Comparison | Strict coverage | Blocked post years |\n|---|---:|---|\n'
cv=pd.read_csv(R/'outputs/pisco_biomass_borrowing_20261102/biomass_coverage.csv')
for t in cv.itertuples():text+=f'| {t.series_id.replace(chr(124), " / ")} | {t.borrowing_post_years}/{t.original_post_years} | {t.unsupported_years} |\n'
text+='''
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
'''
for p in json.loads((O/'examined_data_files.json').read_text()):text+='- `'+p['path']+'`\n'
text+='''
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
'''
for t in x.itertuples():text+=f'| {t.row_id} | {t.series_id.replace("|"," / ")} | {t.site} | {t.status} | {t.year} | {t.zone} | {str(t.date)[:10]} | {t.post} | {t.root_reason} |\n'
(O/'report.md').write_text(text)
print('Wrote report.md and findings.csv',len(text),len(f))
