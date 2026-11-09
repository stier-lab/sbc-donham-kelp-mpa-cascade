from pathlib import Path
import pandas as pd,json,csv,hashlib
O=Path(__file__).resolve().parent;R=Path('/Users/adrianstier/sbc-donham-kelp-mpa-cascade');M=Path('/Users/adrianstier/sbc-donham-kelp-mpa-manuscript')
def md(df):
 return '| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n'+'\n'.join('| '+' | '.join(str(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))
checks=json.loads((O/'checks.json').read_text());ref=pd.read_csv(O/'reference_by_program.csv');resp=pd.read_csv(O/'reference_by_response.csv');dep=pd.read_csv(O/'pisco_depth_by_response.csv');ds=pd.read_csv(O/'density_scope_differences.csv').fillna('')
report='''# Coverage, depth, timing and references: final read-only audit

2026-11-02. Parent session: `393ab715dc104db49d989bacdd9f2f3e`. All agent writes are in `/tmp/donham_swarm_20261102/coverage`. No source, primary output, guide, manuscript or shared-status edits; no messages, pushes or additional agents. Recommendations below are not adopted scientific rules.

## Actionable result

The reference and timing arithmetic checked cleanly. Preserve the current scoped rules and researcher exclusions; no new global depth, date-gap or duration cutoff is supported by this audit. The remaining concrete guide change is to identify the roster behind the six-comparison density/shared-zone review: it uses the biomass station/year roster, preserving recovery coverage but not the complete original density baseline. Add a program/response coverage table rather than implying all depth/reference questions are resolved.

Three rule-link defects were verified in the initial decision table: missing R04 for six density histories, missing R05 for 34 retained PISCO fish histories, and R03 incorrectly attached to 18 satellite candidates (17 retained, one excluded). Parent reports these fixed. I independently read the corrected renderer logic; the generated `comparison_decisions.csv` still had the old links at my final spot check. These are therefore **renderer fixes verified; regenerated-table verification pending**, not requests for additional scientific changes. Parent owns rendering and final consolidated verification.

## Completed independent checks

- Recomputed all 175 station/reference coverage summaries from 12,840 station-response-year rows; reference union, full-record intersection, post-only intersection and changing-membership flags all agree.
- Independently recomputed 4,082 available/persistent annual log ratios from arithmetic role means and frozen offsets; maximum absolute error 4.45e-16. This certifies aggregation, not donor validity or causality.
- Rebuilt all 417 PISCO fish station-response-year zone intersections from complete physical survey units; exact agreement.
- Read the original PISCO swath and fish deliveries and reconstructed method-preserving depth-event availability across all 105 retained PISCO response histories. Inspected headers of the other eight source deliveries, spanning all ten inventory datasets; no claim of a fresh full-content audit of those eight.
- Independently enumerated date assignments for 8,292 station-pair/zone/year/tolerance strata, checking maximum cardinality then minimum total gap. All 4,206 saved pair rows have correct gaps and unique within-pair use; all 36 timing coverage rows agree. No repository matching/model code was executed.
- Checked all-month status applies to exactly twelve response histories (six MPA/taxon comparisons, biomass and density). Current README/current_rules appropriately restrict adoption to these comparisons and keep primary promotion on hold.

Reproduce with `python3 /tmp/donham_swarm_20261102/coverage/check_coverage.py`. `checks.json`, `check_run.log`, `input_hashes.json` and detailed CSVs record results and exact inputs. The first run stopped at a scratch-only string/integer transect-ID join mismatch; normalized IDs and the two subsequent full runs completed successfully. No fitting jobs remain.

## Available versus constant references

Counts below are **response/program/MPA histories** and **post-protection comparison-years**, not independent sites or samples. Full-record fixed references must occur in every approved year, including unused baseline/designation years; post-only fixed references need occur in every observed post-protection year. The protected mean, years and offset remain unchanged wherever a reference intersection exists.

'''
q=ref[['Source','histories','changing_full','changing_post','single_reference','no_full_fixed','no_post_fixed','post_years','post_years_full_fixed','post_years_post_fixed']]
report+=md(q)+'''

Available references retain 2,185 post comparison-years. Requiring full-record persistence retains 1,897 (288 lost); post-only persistence retains 1,979 (206 lost). Full persistence leaves 24 histories with no reference intersection (23 PISCO, one KFM). The 18 rows with an empty post-only set include **two Point Dume biomass histories with zero post years**, so only 16 histories with observed recovery years lose all post-persistent references. Distinguish structural absence of post observations from loss under a reference rule.

Post-only versus full-record persistence recovers 82 comparison-years in eight PISCO histories: Abalone sheephead biomass/density (22 years total) and Point Vicente kelp, red/purple biomass/density and lobster density (60 total). LTER Campus Point sheephead has stable recovery references; two full-record turnover flags reflect earlier coverage, not recovery-period turnover. Stable membership still does not validate habitat matching. All 52 single-reference histories lack a leave-one-reference-out diagnostic; they are not automatically ineligible.

Program/response detail:

'''+md(resp[['Source','Taxa','Resp','histories','changing_post','no_full_fixed','no_post_fixed','post_years','post_years_full_fixed','post_years_post_fixed']])+'''

These counts describe the dated all-program station reconstruction, not the later corrected/borrowed biomass cohort. In particular, original Point Vicente kelp 2023 and donor-alias biomass limitations are not cured by successful ratio reconstruction. Do not relabel these estimates as current verified biomass support.

## Depth and zone losses

PISCO raw availability checks find at least one missing positive numeric depth in 40 of 105 response histories, totaling 733 of 53,688 **event uses across histories**. These are not unique dives: the same events recur across taxa, responses and sometimes MPA assignments. Missing depths are metadata flags, not demonstrated habitat mismatches. Fish events combine level/observer records, so multiple recorded depths within one event are not automatically conflicts. No numeric-depth trimming was adopted or outcome-refitted.

'''+md(dep)+'''

The raw swath column above is an **all-month availability diagnostic outside the approved six comparisons**, not a global authorization to change seasonal filters. Fish uses approved complete units. Under unchanged full station/year membership, no constant zone exists for six PISCO response histories: Abalone sheephead biomass/density, and Point Vicente kelp biomass plus red/purple/lobster density. The latter four differ from the six-comparison result because their full density/kelp records include additional earlier station-years.

The fish audit independently confirms Abalone lacks an all-station shared zone in 2010 and 2014; deleting those years would be a new coverage decision. The constant fish zones narrow to OUTER at Anacapa, Point Vicente and South Point and INNER at Santa Barbara Island. These are protocol strata, not interchangeable numerical depths; a blanket INNER/MID policy would remove valid fish coverage. Do not conflate CNMD fish survey level with the MID habitat zone.

The older measured-shallow/90-day/five-year proposal retains MPA-year unions of Abalone 6, Gull 16, Point Dume 6, Vicente 11, South Point 14. These are years with at least one qualifying pair, **not complete-reference annual years**, and therefore cannot be compared directly to the newer timing grid. Its 15/16.25/20-m alternatives have identical saved pair/year membership in this scoped proposal. Allowing unknown depths retains one additional Gull year (2004), and one more South Point reference with no additional MPA-years. These are historical sensitivity costs only; their numeric cutoffs were not adopted. Program-wide station/zone/depth-standardized outcome refits remain unfinished.

## Six-comparison timing and baseline scope

All 63 original post-protection density comparison-years survive the all-month review. At date tolerances 0/14/30/60/90/366 days, complete station-pair plus INNER/MID coverage retains Abalone 0/0/1/1/2/11 years and Point Vicente 0/0/0/1/6/10 years **for each of three taxa**. The 90-day version therefore retains 24 of 63 taxon comparison-years, loses 39, and leaves Abalone too short for a fitted trend. No cutoff is recommended. Maximum-cardinality pairing is correctly implemented; date matching does not recover verified campaigns or remove seasonal confounding.

The density rerun deliberately restricts its baseline to biomass roster support (`refit_pisco_all_months.py`: `base` filtered to `v` years and station rows built from biomass IDs). Here are the exact losses relative to the original density roster:

'''+md(ds)+'''

The 33 omitted density station-year uses are all prior to designation, across four density histories; eight density comparison-years are omitted in total (Abalone lobster 2011; Vicente red/purple 2008 each; Vicente lobster 2008–2010). **Correction: this list sums to six density comparison-years, not eight**: 1+1+1+3=6. No recovery-year or designation-year station uses are lost. The raw shared INNER/MID assertion is verified on the biomass-derived roster, not the full original density baseline. Minimal guide fix: add `roster_basis`, `baseline_years_omitted` and `post_years_preserved` or equivalent explicit prose to the twelve relevant status rows. Restore/use density-specific baseline records in any future before–after integration; do not silently generalize a recovery review to baseline eligibility.

The historical depth code treats `survey_year` as a campaign. Raw retained-site swath data contain 55 site/calendar-year/survey-year/zone groups with multiple dates, spanning up to 293 days (e.g., South Point West/Cluster Point North, January–November 2014). A campaign label may span dates, but this does not verify a single comparable field visit. The newer timing report correctly says campaign IDs are unverified. Also, old depth scenario `d0` means unrestricted gap, whereas new timing `tolerance_days=0` means exact-date matching. Rename/annotate these when presenting historical results; do not pool those labels as the same scenario.

## Program-specific recommendations (not new adopted rules)

- **PISCO swath:** keep all valid months and constant INNER/MID only in the six approved comparisons, with equal visit/zone/station weighting and full reference support. Preserve their 63 density recovery-years. Label the restricted baseline roster. For other swath histories, show constant-zone/station losses and measured-depth availability before extending the rule; retain unknown depth pending metadata. Separate count support from size support.
- **PISCO fish:** retain the approved CRANE, repeat-visit and complete-pass rules. Use fish-specific zone definitions and actual available station-years. Report Abalone's two zone-gap years and location-specific constant-zone subsets. Do not impose benthic depth targets or automatic INNER/MID-only selection on fish.
- **KFM:** retain permanent physical site/transect identities where supported by field methods, separate visual-fish from benthic protocols, and preserve researcher exclusions. Compare available, full-record and post-only sets; one Scorpion sheephead density history loses all fixed references. Count-file schemas contain no survey-depth field; later sensor-site depths are not survey-year quadrat depths. Obtain compatible metadata if available; do not exclude 35 retained histories because that field is absent. Keep unresolved size-donor identities separate.
- **SBC-LTER:** all 18 histories have post-stable reference station membership, but unit/transect stability and habitat equality are separate checks. Retain the documented expert assignments provisionally. Full-record restrictions can change Campus Point sheephead's represented reference population without diagnosing recovery turnover. Preserve fish 80-m², lobster 300-m² and kelp/benthic sampling identities; no global seasonal cutoff or borrowed PISCO depth labels. Missing count-file depth metadata applies to all 18 histories, not proof that 18 comparisons are invalid.
- **Landsat:** diver depth matching is not applicable. Retain verified masks and nonempty polygons, distinguish changing polygon availability from unequal reference footprints, and retain available-versus-fixed and omission diagnostics. Three histories change polygon membership over the full record, one during recovery. Preserve annual maximum/quarterly support interpretation; do not impose diver-date windows or equate polygon totals to station densities.
- **Shared size and KFM plant-size datasets:** these are donor measurements, not extra reference-station replication. Keep exact-date and annual support distinct; numeric size-depth fields do not supply missing count-survey depths. Conditional borrowing and six KFM alias histories remain the other agents' substantive scope.

Cached incoming correspondence confirms Scott Hamilton's habitat/geography rationale and Emily's distinction between same-site before–after and after-only comparisons. Her statement is recollection, not implementation verification. The cached November 1 evidence record likewise does not authorize 90 days, 16.25 m or five years. Do not infer site deletion from missing rationale or environmental concerns.

## Exact guide/table changes needed

1. Parent's R04/R05/R03 renderer repairs: regenerate and verify `comparison_decisions.csv`, the viewer and any manifests; no scientific rerun.
2. Add the density-roster/baseline limitation to `audit_status_by_comparison.csv`, current guide README and R03/R04 limitation text, with the exact six-comparison scope and `density_scope_differences.csv` evidence. Avoid suggesting the full before–after baseline has passed the shared-zone rule.
3. Add response-specific reference coverage counts to a linked table, including full/post period distinction and zero-post-year N/A status. Current R06 is appropriately provisional; no change to selection policy is required.
4. Keep U02 (KFM/LTER depth), U03/U04 (reference rationale), U06 (broader coverage policy) and U11 (satellite design) open. Clearly distinguish full-history data audit from final physical habitat validation. Older depth proposal needs a historical/not-adopted banner and `survey_year`/`d0` definitions if it remains linked as decision evidence.

## Unfinished and limits

No new outcome model, total-uncertainty analysis, final program-wide depth standardization, all-program seasonal adjustment, physical-depth equivalence validation, missing campaign metadata recovery or before–after baseline rebuild was attempted. No need to run them to adopt the documentation corrections. Source-specific methods were read through the cached project evidence and code; no fresh external methods retrieval is claimed. No fresh full contents check of all eight non-PISCO deliveries or source-donor correction is claimed. Parent's regenerated documentation remains to be verified. Available source evidence cannot resolve unrecorded field metadata; no automatic exclusion is recommended.

`examined_files.csv` lists exact read inputs with inspection scope; `input_hashes.json` records computational input bytes. Pre-existing repo modifications were observed and left untouched. The agent did not open an overlapping continuity session or update parent status.
'''
report=report.replace('eight density comparison-years are omitted in total','six density comparison-years are omitted in total').replace(' **Correction: this list sums to six density comparison-years, not eight**: 1+1+1+3=6.','')
(O/'report.md').write_text(report)
findings=[
('C01','Six PISCO density histories','consolidate_current_decisions.py; rule_link_discrepancies.csv','Initial table omitted R04 for six density counterparts; parent renderer now corrected','Parent regenerate tables and verify all twelve links','fixed_in_renderer_pending_render'),
('C02','34 retained PISCO fish histories','consolidate_current_decisions.py; rule_link_discrepancies.csv','Initial table omitted R05; parent renderer now corrected','Parent regenerate and verify 34 retained fish links','fixed_in_renderer_pending_render'),
('C03','18 satellite candidates / 17 retained','consolidate_current_decisions.py; rule_link_discrepancies.csv','Diver-zone R03 attached to satellite rows; parent renderer now corrected','Parent regenerate and verify satellite R03 removal','fixed_in_renderer_pending_render'),
('C04','Six-comparison density baseline','refit_pisco_all_months.py; density_scope_differences.csv','33 predesignation density station-year uses and six density comparison-years omitted because review uses biomass roster; all 63 post comparison-years preserved','Label roster scope now; preserve density-specific baseline in future BA integration','verified_action_needed'),
('C05','175 station histories','reference_checks.csv; checks.json','All reference summaries and 4082 arithmetic ratio cells reproduce; max error 4.44e-16','Keep aggregation claims; do not equate them with donor/habitat certification','verified_clean_negative'),
('C06','All-program reference policy','reference_by_program.csv; reference_by_response.csv','Full fixed reference sets lose 288/2185 post comparison-years; post-only sets lose 206; 24 histories lack full intersection','Publish period-specific coverage; keep available-reference and constant-set sensitivity pending selection','verified_recommendation_not_adopted'),
('C07','PISCO depth/zone support','pisco_depth_by_response.csv; fish_zone_year_checks.csv','40/105 histories contain missing positive numeric depths; full-history constant zones absent in six response histories','Retain metadata flags; use response-specific zones; do not impose global cutoff','verified_recommendation_not_adopted'),
('C08','Six PISCO timing comparisons','timing_checks.csv; checks.json','4206 date-pair rows and 36 coverage rows validated; 90-day complete coverage retains 24/63 comparison-years','Keep timing sensitivities without adopting a cutoff','verified_clean_negative'),
('C09','Historical depth proposal','audit_depth_reference_policy.py; survey_year_multi_date_groups.csv','d0 means unrestricted not zero-day; 55 raw site/year/survey-year/zone groups contain multiple dates up to 293-day spread','Label survey_year as unverified campaign grouping and historical thresholds not adopted','verified_action_needed'),
('C10','KFM/LTER and satellite','source_schema_checks.csv; cached correspondence','No count-survey depth fields in inspected KFM/LTER schemas; satellite diver depth not applicable; reference rationale incomplete','Retain provisional researcher assignments and program-specific methods; seek metadata without automatic exclusions','verified_metadata_gap'),
('C11','Current guide adoption scope','current_rules.csv; README.md; checks.json','Current rule scope correctly limits all months to six comparisons and does not promote review into paper','Preserve scope and review hold','verified_clean_negative')]
with (O/'findings.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['finding_id','scope','evidence','finding','recommended_action','status']);w.writerows(findings)
# Computational files plus the exact contextual code/correspondence used.
manifest=json.loads((O/'input_hashes.json').read_text());entries=[]
for p,h in manifest.items():
 scope='header/schema only' if p.startswith('/Users/adrianstier/sbc-donham-kelp-mpa-data') and 'PISCO/' not in p else 'full CSV read; analysis subset documented in check_coverage.py'
 entries.append(dict(path=p,inspection=scope,sha256_at_check=h))
extra=['AGENTS.md','CLAUDE.md','code/tools/consolidate_current_decisions.py','code/tools/refit_all_program_stations.py','code/tools/prepare_all_program_stations.py','code/tools/audit_station_coverage_refits.py','code/tools/refit_pisco_all_months.py','code/tools/check_pisco_timing_donors.py','code/tools/audit_depth_reference_policy.py','code/tools/audit_other_program_footprints.py','code/tools/build_dataset_inventory.py','code/tools/inventory_latest_audits.py','code/tools/report_dataset_inventory.py']
extra+=['outputs/'+p+'/README.md' for p in ['dataset_inventory_20261101','all_program_station_refits_20261102','station_coverage_refits_20261101','pisco_all_months_20261102','pisco_timing_donors_20261102','depth_reference_policy_20261101','seasonal_zone_audit_20261101','anacapa_comparability_20261101']]
extra+=['outputs/dataset_inventory_20261101/'+p for p in ['current_rules.csv','unresolved_issues.csv','user_decisions.csv','freeze_readiness.csv','dataset_inventory.csv','dataset_audit.csv','measurement_protocols.csv','raw_station_response_scan.csv','survey_calendar_campaign_audit.csv']]
context=[R/p for p in extra]+[M/'AGENTS.md',M/'CLAUDE.md',M/'correspondence/coauthor_reviews/2026-11-01_depth_reference/EVIDENCE.md',M/'correspondence/coauthor_reviews/2026-10-29_reply_integration/incoming.json',M/'correspondence/coauthor_reviews/2026-10-29_reply_integration/replies.json']
for p in context:
 if str(p) not in manifest and p.exists():entries.append(dict(path=str(p),inspection='context/code inspection; relevant sections or parsed fields, not blanket full-source certification',sha256_at_check=hashlib.sha256(p.read_bytes()).hexdigest()))
pd.DataFrame(entries).to_csv(O/'examined_files.csv',index=False)
print('Final report and findings written:',len(findings),'findings;',len(entries),'source paths')
