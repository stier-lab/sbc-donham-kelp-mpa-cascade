from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,html
O=Path('outputs/dataset_inventory_20261101');A=Path.cwd()
def read(n):return pd.read_csv(O/(n+'.csv'),dtype={'site_id':str},low_memory=False)
def save(d,n):d.to_csv(O/(n+'.csv'),index=False)
sites=read('site_response_coverage');rawyear=read('raw_station_year_audit');pairs=read('candidate_pair_coverage');comp=read('comparison_eligibility');ds=read('dataset_inventory');mpa=read('mpa_metadata')
from inventory_site_resolution import apply as apply_site_resolution
sites=apply_site_resolution(sites)
from inventory_measurement_review import apply as apply_measurement_review
ds,sites=apply_measurement_review(ds,sites)
fish_review=Path('outputs/pisco_fish_effort_20261101/comparison_summary.csv')
if fish_review.exists():
 fr=pd.read_csv(fish_review).set_index('series_id')
 for field in ['repeat_changes_fitted_recovery','complete_level_changes_fitted_recovery']:
  comp[field]=comp.series_id.map(fr[field])
 comp['fish_effort_audit']=comp.series_id.map(lambda x:'../pisco_fish_effort_20261101/README.md' if x in fr.index else '')

# Supplemental measurements and conversion evidence are indexed without pretending they have been re-extracted here.
aux=[]
for f,role,status in [
('data/external/sbc1275/Macrocystis_Biomass_Relationships_20261016.csv','Kelp monthly conversion coefficients','Li confirmed SBC127.5 Table2 Model1; conversion does not remove ecological seasonality'),
('tables/table_s_kfm_size_inputs.csv','Separate KFM predator-size source manifest','Previous raw matching audit; no new size reconstruction this turn'),
('tables/table_s_primary_size_inputs.csv','Lobster/purple-urchin pool source manifest','15 canonical histories; not all expanded candidates'),
('data/harmonized/landsat_v2_verification.json','Satellite pixel/mask reconstruction manifest','Previously verified58 nonempty arrays; masks are not newly rerun here'),
('data/harmonized/lobster_missing_lengths.json','Five missing LTER lobster lengths correction','Missing flags retained; same-site/date mean mass, not recovered lengths'),
('docs/PRIMARY_SIZE_UNCERTAINTY_2026-10-30.md','Size uncertainty propagation evidence','Existing sensitivity; not a refit for newly proposed survey matching')]:
 p=Path(f);aux.append(dict(source=f,role=role,exists=p.exists(),sha256=hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else '',audit_status=status))
save(pd.DataFrame(aux),'supporting_source_inventory')
ylookup={k:g for k,g in rawyear.groupby(['dataset_id','program','site_id','taxon','method'])}
changes=[]
for _,s in sites.iterrows():
 g=ylookup.get((s.dataset_id,s.program,s.site_id,s.taxon,s.method),pd.DataFrame());d=dict(site_record_id=s.site_record_id)
 if len(g):
  d.update(audited_years=len(g),month_set_changes=g.months.fillna('').nunique()>1,zone_set_changes=g.zones.fillna('').nunique()>1,survey_unit_set_changes=g.unit_ids.fillna('').nunique()>1,years_without_complete_date=int((g.dates.isna()|g.dates.eq('')).sum()),years_with_missing_or_negative_counts=int(((g.missing_counts+g.negative_counts)>0).sum()),years_without_measured_depth=int((g.positive_depth_records==0).sum()))
 changes.append(d)
changes=pd.DataFrame(changes);sites=sites.drop(columns=[c for c in changes if c!='site_record_id'],errors='ignore');sites=sites.merge(changes,on='site_record_id',validate='one_to_one');save(sites,'site_response_coverage')
slookup=sites.set_index('site_record_id')
def dateyears(v):
 out={}
 if pd.isna(v):return out
 for d in v.split(';'):
  t=pd.Timestamp(d);out.setdefault(t.year,[]).append(t)
 return out
match=[]
for _,p in pairs.iterrows():
 a=slookup.loc[p.protected_record_id];b=slookup.loc[p.reference_record_id];aa=dateyears(a.get('survey_dates'));bb=dateyears(b.get('survey_dates'));gaps={y:min(abs((x-z).days) for x in aa[y] for z in bb[y]) for y in set(aa)&set(bb)}
 # Broad all-years source depth range diagnostic; NOT within-zone/visit eligibility.
 depthknown=all(pd.notna(x) for x in [a.depth_min_m,a.depth_max_m,b.depth_min_m,b.depth_max_m])
 match.append(dict(pair_id=p.pair_id,years_with_exact_dates_both=len(gaps),years_nearest_dates_within60=sum(v<=60 for v in gaps.values()),years_nearest_dates_within90=sum(v<=90 for v in gaps.values()),max_nearest_date_gap_days=max(gaps.values(),default=np.nan),broad_depth_ranges_overlap=(max(a.depth_min_m,b.depth_min_m)<=min(a.depth_max_m,b.depth_max_m)) if depthknown else np.nan,date_depth_limit='All raw methods/recorded years within this row; no campaign/stratum/effort validation, not primary eligibility'))
match=pd.DataFrame(match);pairs=pairs.drop(columns=[c for c in pairs if any(c==k or c.startswith(k+'_') for k in match if k!='pair_id')]);pairs=pairs.merge(match,on='pair_id',validate='one_to_one');save(pairs,'candidate_pair_coverage')
# Explicit unresolved issues: distinguish missing evidence from an analysis choice.
issues=[
('U01','KFM','Permission and redistribution','Scott Gabara said Oct29 that park permit and upload conditions need confirmation; no newer answer found in targeted searches.','Contributor confirmation','Scott Gabara / park data manager','Keep raw records private; confirm before submission/public archive','docs/KFM_PREDATOR_SIZE_FEASIBILITY_2026-10-29.md'),
('U02','KFM and LTER','Physical depth/habitat match','Count deliveries lack compatible survey-depth fields; fixed station/transect IDs do not establish depth equivalence.','Missing metadata','Program investigators / field metadata','Recover survey-level depth if available; disclose if unavailable','outputs/seasonal_zone_audit_20261101/README.md'),
('U03','LTER','Individual reference rationale','Emily confirms consultation with Bob/Li; station-by-station ecological reasons remain undocumented.','Contributor confirmation','Emily / relevant LTER investigators','Retain documented map and state limits; do not call matching validated','docs/SELECTION_CRITERIA_AUDIT_2026-10-29.md'),
('U04','PISCO Naples','Community matching evidence','Jenn confirms diver expertise and nMDS; underlying ordination not recovered, and no close match exists.','Missing supporting analysis','Jenn / Spiecker technical-report materials','Cite qualitative basis; avoid claiming independent quantitative validation','docs/SELECTION_CRITERIA_AUDIT_2026-10-29.md'),
('U05','PISCO','Three missing source measurements','Point Vicente kelp fronds, Scorpion red-urchin count and South Point sheephead lengths remain unrecovered.','Missing measurements','Original field records / program data managers','Retain missing flags and complete-observation sensitivities','outputs/pisco_source_reconciliation_20261101/README.md'),
('U06','PISCO','Season/depth matching and reference dependence','90-day, shallow-depth and five-year rules remain proposals. Abalone depends on date tolerance; Dume changes season/depth and has a long gap.','Analysis decision','Adrian / analysis review','Freeze coverage rules before outcome refits; use shared-reference dependence and common prediction support','outputs/depth_reference_policy_20261101/README.md'),
('U07','KFM','Predator biomass support','Separate sheephead sizes cover a limited period with footprint/sparsity limitations; no verified lobster sizes in inspected releases.','Missing measurements / population match','Scott Gabara','Keep density distinct; no fabricated full-history biomass','docs/KFM_PREDATOR_SIZE_FEASIBILITY_2026-10-29.md'),
('U08','PISCO/KFM/LTER','Sparse and shared size pools','Canonical lobster/purple-urchin pools audited; expanded candidate coverage and newly matched depths/years need corresponding size linkage.','Analysis work','Analysis pipeline','Extend size-pool ledger and shared-sample propagation to final roster','docs/PRIMARY_SIZE_POOL_INVENTORY_2026-10-30.md'),
('U09','All','Historical protection and sampled boundaries','Network dates and current restriction snapshots do not establish first protection or historical membership of every transect. Five previously unverified current take classes were resolved against CDFW pages; historical rules still require separate evidence.','Missing historical metadata','Official archived designations / program records','Separate network date, prior protection, legal class and ecological scope','mpa_metadata.csv'),
('U10','PISCO','LECHUZA coordinates','Source map supplies conflicting coordinates; distance must not be chosen by row order.','Conflicting metadata','Source mapping / program data manager','Resolve coordinate identity before distance-based selection','outputs/depth_reference_policy_20261101/coordinate_conflicts.csv'),
('U11','Landsat','Reference design and interpretation','Masks/annual arrays verified; approximate boundary band and historical clustering recovery do not establish causal comparability. Pre2012 selection is not predesignation for island MPAs.','Design limitation','Analysis review / Tom if historical detail needed','Keep saved masks; explain design and avoid silently rerunning clustering','docs/LANDSAT_V2_VERIFICATION_2026-10-21.md'),
('U12','All','Unmapped source records','19 formerly unmatched assignments are resolved group labels linked to 40 existing stations; only Trailer/Willows remain unassigned. Other-region/method source groups remain outside automatic eligibility.','Scope/crosswalk review','Analysis pipeline then contributors if needed','Use explicit group-to-station links; keep Trailer/Willows unassigned unless source evidence is recovered','site_group_resolution.csv'),
('U13','Size surveys','Incomplete exact dates','Many size records lack a valid complete date although year/month fields exist; annual coverage must not be discarded or claimed as exact visit matching.','Metadata limitation','Source size metadata','Keep annual size support and exact-date support separate','raw_dataset_scan.csv'),
('U14','All','Repeated records','Selected-field duplicates omit some taxon/observer/replicate details and may be legitimate equal-sized individuals or repeated samples.','Record identity screen','Source schema / field keys','Do not deduplicate by values or row order; resolve full physical keys before modifying data','raw_dataset_scan.csv')]
issue=pd.DataFrame(issues,columns=['issue_id','scope','issue','finding','type','who_or_where','next_action','evidence']);save(issue,'unresolved_issues')
# Comparison station lists contain observed station identifiers, not parent/map-only entries.
for col in ['mapped_protected_sites','mapped_reference_sites']:
 old='historical_'+col+'_including_map_labels'
 if old not in comp:comp[old]=comp[col]
for idx,r in comp.iterrows():
 st=sites[(sites.MPA==r.MPA)&(sites.program==r.Source)]
 obs=st[st.inventory_record_type=='Observed source/response record']
 comp.loc[idx,'mapped_protected_sites']=';'.join(sorted(obs[obs.role.isin(['mpa','protected'])].site_id.astype(str).unique()))
 comp.loc[idx,'mapped_reference_sites']=';'.join(sorted(obs[obs.role=='reference'].site_id.astype(str).unique()))
 comp.loc[idx,'mapped_group_labels']=';'.join(sorted(st[st.inventory_record_type.str.startswith('Grouped')].site_id.astype(str).unique()))
 comp.loc[idx,'station_list_scope']='Observed candidate station IDs only; exclusions remain explicit and these lists are not fitted membership. Group labels linked separately.'
 comp.loc[idx,'group_resolution_table']='site_group_resolution.csv'
# Link issues to every comparison; broad labels represent unverified scope, not known defects everywhere.
links=[]
for _,r in comp.iterrows():
 ids=['U09']
 if r.Source in ['KFM','LTER']:ids+=['U02']
 if r.Source=='KFM':ids+=['U01']+(['U07'] if r.Taxa in ['S. pulcher','P. interruptus'] else [])
 if r.Source=='LTER':ids+=['U03']
 if r.Source=='PISCO':ids+=['U06']+(['U04'] if r.MPA=='Naples SMCA' else [])+(['U10'] if r.MPA=='Point Dume SMCA' else [])
 if r.Source=='Landsat':ids+=['U11']
 if r.Resp=='Bio' and r.Taxa in ['P. interruptus','S. purpuratus','M. franciscanus']:ids+=['U08','U13']
 if (r.MPA,r.Source,r.Taxa) in [('Point Vicente SMCA','PISCO','M. pyrifera'),('Scorpion SMR','PISCO','M. franciscanus'),('South Point SMR','PISCO','S. pulcher')]:ids+=['U05']
 links.append(';'.join(ids))
comp['issue_ids']=links;save(comp,'comparison_eligibility')
prev=read('supplement_table_preview');prev['issue_ids']=links
if fish_review.exists():
 pk=prev[['MPA','Taxa','Source','Resp']].agg('|'.join,axis=1)
 for field in ['repeat_changes_fitted_recovery','complete_level_changes_fitted_recovery']:
  prev[field]=pk.map(fr[field])
save(prev,'supplement_table_preview')
# Program summary: count comparison histories, not stations or independent samples.
progs=[]
for k,g in comp.groupby('Source'):
 st=sites[sites.program==k];progs.append(dict(program=k,current_paper_histories=int(g.in_current_manuscript.sum()),reconstructed_candidates=len(g),researcher_excluded=int(g.exclude.sum()),retained_review_candidates=int((~g.exclude).sum()),mapped_station_assignments=len(st[['MPA','site_id','role']].drop_duplicates()),observed_station_assignments=len(st[st.inventory_record_type=='Observed source/response record'][['MPA','site_id','role']].drop_duplicates()),retained_with_three_prior_years=int(((~g.exclude)&g.three_strict_prior_years).sum()),retained_passing_proposed_recovery_duration=int(((~g.exclude)&g.proposed_recovery_duration_pass).sum())))
prog=pd.DataFrame(progs);save(prog,'program_summary')
def table(d):
 return '| '+' | '.join(map(str,d.columns))+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+''.join('| '+' | '.join(str(x).replace('|',' / ').replace('\n',' ') for x in r)+' |\n' for r in d.itertuples(index=False,name=None))
counts=json.loads((O/'build_checks.json').read_text())
text=f'''# Dataset, site and comparison audit

November 1, 2026. Review only. No manuscript, supplement, canonical data, models or source records changed. No messages sent. No inclusion policy was adopted.

## Summary

The audit inventories all four programs used in the paper: PISCO, KFM, SBC-LTER and Landsat. It reads nine biological source CSVs plus the annual satellite delivery ({counts['source_rows_scanned']:,} source rows; the large integrated files are narrowed to focal taxa after reading). Supporting coefficient, size-source and satellite reconstruction manifests are indexed separately. These files are the local delivered releases, not a claim that no newer release exists.

The tables retain all 337 mapped candidate station assignments, expanded to {len(sites):,} source/response/method rows, all 206 reconstructed candidate comparisons, and all 143 current manuscript histories. Some stations occur in more than one MPA assignment; assignment counts retain cross-program rows from the inherited master map and are not counts of PISCO survey sites. Rows are not independent replicates. There are {len(pairs):,} candidate raw same-year pair records, before ecological exclusions or depth/date qualification.

'''+table(prog)+'''
The 175 retained candidates are the previous researcher/ecological scope roster, not a newly approved primary cohort. Three prior years and five recovery years are coverage screens only. Neither establishes baseline stability, habitat equivalence, identifiable recovery or sufficient size measurements. Designation-year observations are counted separately; existing model routing is not replaced by this table.

## Read the tables

- [Dataset-by-dataset measurement audit](DATASET_AUDIT.md), [screen results](measurement_design_scan.csv), and [additional checks before freezing](additional_analysis_checks.csv). These distinguish verified source flags from untested numerical consequences.
- [Dataset inventory](dataset_inventory.csv): source/version/hash, sampling design, measured quantity, units, missingness conventions and data-use status.
- [Resolved site groups](site_group_resolution.csv) and [group-to-station links](site_group_links.csv) distinguish group records from sampled stations and preserve exclusions.
- [Site and response coverage](site_response_coverage.csv): every mapped station assignment, source method/response, actual years, internal gaps, months, depth fields, coordinates and changes in sampling coverage. Raw coverage includes observations outside historical seasonal filters; it is not final biomass-valid coverage.
- [Comparison eligibility](comparison_eligibility.csv): all reconstructed candidates, current-paper membership, researcher exclusions, baseline/recovery coverage, reference evidence, size limitations and issue IDs.
- [Compact supplement table preview](supplement_table_preview.csv): one row per comparison, including excluded candidates. This is a review artifact, not inserted into the supplement.
- [Candidate pair coverage](candidate_pair_coverage.csv): separate protected/reference identities with same-year coverage and nearest-date diagnostics. These are raw opportunities, not approved pairs. Broad all-years depth-range overlap is not a matched-stratum test.
- [Full satellite delivery](satellite_full_delivery.csv) includes all 24 supplied MPAs and all 72 polygon slots, with 14 empty masks explicitly marked unavailable.
- [Satellite candidate coverage](satellite_candidate_coverage.csv) covers all delivered polygons, including those without reconstructed comparisons.
- [MPA metadata](mpa_metadata.csv): dates, current-law snapshots, historical uncertainty and take-by-taxon fields. A reference row inherits the comparison MPA date for analysis context, not a claim the reference was protected.
- [Unresolved issues](unresolved_issues.csv): specific missing evidence, analysis choices and next actions.
- [Supporting sources](supporting_source_inventory.csv), [unmapped source groups](unmapped_source_records.csv) and [annual sampling audit](raw_station_year_audit.csv) retain the evidence behind coverage flags.

## Main findings

1. **The programs measure different things.** Diver density, size-derived standing biomass and satellite surface-canopy biomass are kept distinct. Fish, benthic surveys and separate size collections do not inherit one another's coverage.
2. **Source-map namespaces matter.** Of 64 assignments without a direct source match, 43 are KFM/LTER rows embedded in the PISCO master map; they are not evidence of missing PISCO surveys. The subsequent hierarchy audit resolves 19 of the remaining 21 rows as grouped-site labels linked to 40 existing stations; only Trailer and Willows remain unassigned. Group labels carry no additional observations; identity and exclusion are separate. [Station exceptions](station_link_exceptions.csv) and [resolved hierarchy](site_group_resolution.csv); the earlier [legacy-name-only check](alias_evidence.csv) is historical distinguish these cases. Provider campaign years also differ from calendar years in 130 distinct metadata cells; [the date ledger](survey_calendar_campaign_audit.csv) preserves both.
3. **Long time series can conceal changing sampling.** Raw annual station metadata now identifies changes in months, zones and unit labels. Changes are flags, not automatic exclusions: non-fixed transects can have legitimate changes in identifiers. Depth summaries are record-weighted descriptions, not area-weighted habitat estimates.
4. **A year with no record is not a biological zero.** Counts, explicit zeros, negative sentinels and missing measurements remain separate. Incomplete full dates retain a recorded year/month. Repeated selected-field values are never silently removed.
5. **Some reference rationales remain incomplete.** Scott supports specific PISCO habitat exclusions; Jenn supports the Naples reference choice but reports no close match; Emily confirms LTER consultation but not each station's rationale. This does not validate every inherited pairing. The user's newer exclusion of Gull KFM supersedes older documents that retained it; Gull PISCO remains a separate candidate.
6. **Size support is uneven.** Existing detailed lobster/purple-urchin inventories cover 15 canonical histories, with sparse/shared pools explicitly linked. This cannot certify expanded candidate years. KFM counts do not become a full biomass history simply because a limited separate size survey exists.
7. **Satellite mask discrepancies are resolved in the corrected delivery.** Existing verification covers 24 MPA and 34 nonempty reference arrays. This audit reads the current annual delivery and links that verification rather than pretending to repeat pixel reconstruction. Known empty reference slots remain unavailable, not zeros. Surface pixels are not assigned dive-depth fields.
8. **Five current protection labels were resolved.** The audit now distinguishes Point Vicente and Casino Point no-take SMCAs from the limited-take rules at Abalone Cove, Arrow Point–Lion Head Point and Crystal Cove. [Dated source links](protection_metadata_resolutions.csv) accompany each correction. These metadata corrections do not change the cohort.
9. **Historical metadata and proposed eligibility remain distinct.** Current protection labels do not reconstruct historical boundaries. The depth/date/five-year recommendation is still a proposal, not a filter retrospectively selected to preserve an effect.

## Primary-analysis freeze

See [freeze readiness](FREEZE_READINESS.md) for required decisions and checks, conditional comparison-level limits, and items that can remain disclosed limitations. Resolving a site name does not adopt new analysis rules.

## Questions still open

'''+table(issue[['issue_id','scope','issue','type','who_or_where']])+'''
Detailed findings, next actions and evidence are in unresolved_issues.csv. Do not send new questions already answered: Li specified SBC127.5 Table2 Model1; Jenn explained the Naples matching rationale; Scott explained Harris/Cuyler and Yellowbanks; Emily explained intended station continuity. The remaining needs are narrower: missing measurements, unrecovered supporting records, historical boundaries, data permission and final analysis choices.

## Email and source checks

Targeted Gmail searches on November 1 covered recent replies from Emily, Jenn, Scott Hamilton, Scott Gabara, Li and Tom. The latest relevant reference replies found were October 29; no newer response to these specific questions appeared in the searched results. Scott's full permission reply was reread: park permission and redistribution terms remain to be confirmed. Li's full reply confirms the fixed kelp coefficients. Earlier reference replies were read in the preceding depth audit; its private dated evidence record is retained in the manuscript repository. Search coverage is not proof that no other private conversation occurred. No private email bodies are copied into this analysis report.

Historical inclusion labels in the station table are explicitly named historical; they do not override the current researcher-scope decision column. The older 145-history coverage snapshot is not used as the current primary roster; current membership is rebuilt from the live 143-history result keys. Exact raw sampling dates, source program and method are kept separate from the provider campaign year.

## Reproduction and boundaries

Run audit_all_dataset_sources.py, build_dataset_inventory.py, extend_dataset_audit.py and report_dataset_inventory.py from code/tools, then verify_dataset_inventory.py. The scan is read-only and writes only this audit directory. Original source hashes and table input hashes are recorded. Auxiliary files are indexed, not automatically certified by their filenames. No raw individual animal records or emails are copied into the deliverable.

The audit is a systematic inventory and comparability screen across the paper's datasets. It does not certify every historical field convention, reconstruct unprovided measurements, validate every site pairing or refit every biological result. Unresolved items are visible instead of being treated as passed checks. The next analysis should use this inventory to freeze a final cohort and verify response-specific survey and size support before any manuscript update.
'''
text='''# Consolidated audit — November 2, 2026

Use audit_master.csv as the consolidated comparison record (206 candidates;175 retained review histories). audit_status_by_comparison.csv tracks all nine audit domains. The viewer now includes Audit status, Station checks and Decisions tabs. Latest annual coverage is labeled audited_; original inventory coverage remains separately recorded. Station/reference findings are integrated, but depth/timing, expanded size support and final model integration remain open. No manuscript changes.

[Latest station/reference findings and figures](../all_program_station_refits_20261102/index.html). Historical November1 source and email checks below retain their original dates; no new email search is claimed.

'''+text
(O/'README.md').write_text(text)
# Standalone searchable table viewer, no remote dependencies.
views={'Summary':prog,'Dataset audit':read('dataset_audit'),'Additional checks':read('additional_analysis_checks'),'Freeze readiness':pd.read_csv(O/'freeze_readiness.csv'),'Site groups':pd.read_csv(O/'site_group_resolution.csv').fillna(''),'Comparisons':prev,'Datasets':ds[['dataset_id','programs','source_file','release_label','sampling_design','measured_response','units','missing_zero_convention','data_use_status','audit_id','measurement_audit_status','measurement_audit_next_action']],'Sites':sites[['program','MPA','site_id','role','source_namespace','record_link_status','inventory_record_type','group_names','linked_constituent_stations','resolved_group_scope','inherited_pisco_reference_exclusion','inherited_kfm_nonfish_exclusion','dataset_id','taxon','method','first_year','last_year','n_years','n_missing','longest_internal_gap','months','depth_min_m','depth_max_m','review_scope']],'Questions':issue}
from inventory_latest_audits import apply as apply_latest_audits
views=apply_latest_audits(comp,prev,sites,ds,issue,views)
from consolidate_current_guide import apply as consolidate_current_guide
views=consolidate_current_guide(views)
import subprocess,sys
subprocess.run([sys.executable,'code/tools/build_coauthor_map_layers.py'],check=True)
subprocess.run([sys.executable,'code/tools/build_coauthor_decision_review.py'],check=True)
subprocess.run([sys.executable,'code/tools/build_methods_results_reference.py'],check=True)
payload={k:json.loads(v.replace({'':None,np.nan:None}).to_json(orient='records',double_precision=15)) for k,v in views.items()}
from inventory_viewer import build
page=build(payload)
(O/'index.html').write_text(page)
print(prog.to_string(index=False));print('Issues',len(issue),'site rows',len(sites),'pair rows',len(pairs))
