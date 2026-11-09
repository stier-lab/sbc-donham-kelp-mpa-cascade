"""Consolidate dated review evidence without promoting it to the manuscript."""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
O=Path('outputs/dataset_inventory_20261101');N=Path('outputs/all_program_station_refits_20261102')
def apply(comp,prev,sites,ds,issue,views):
 used=[]
 def read(p):
  used.append(p);return pd.read_csv(p,low_memory=False)
 coverage=read(N/'coverage_summary.csv');checks=read(N/'reconstruction_checks.csv');results=read(N/'comparison_summary.csv');om=read(N/'reference_omission.csv');annual=read(Path('outputs/pisco_fish_rules_20261101/approved_annual_input.csv'));zones=read(Path('outputs/station_coverage_refits_20261101/station_zone_summary.csv'))
 fields=['n_reference_stations','n_persistent_references','reference_membership_changes','post_reference_membership_changes','n_post_persistent','recovery_fit_status','size_support_status']
 latest=coverage[['series_id']+fields].rename(columns={x:'audit_'+x for x in fields})
 latest=latest.merge(checks[['series_id','verified','max_ratio_error']].rename(columns={'verified':'annual_reconstruction_verified','max_ratio_error':'annual_reconstruction_error'}),on='series_id',validate='one_to_one')
 rescols=['Approved','Persistent_reference','Persistent_post_reference','Log_station_unadjusted','Log_station_adjusted','station_adjustment','reference_subset_change'];latest=latest.merge(results[['series_id']+rescols].rename(columns={x:'recovery_'+x for x in rescols}),on='series_id',how='left',validate='one_to_one')
 omstats=om.groupby('series_id').agg(reference_omissions_checked=('omitted_reference','size'),largest_absolute_omission_shift=('difference',lambda x:x.abs().max())).reset_index();latest=latest.merge(omstats,on='series_id',how='left',validate='one_to_one');latest['reference_omissions_checked']=latest.reference_omissions_checked.fillna(0).astype(int)
 records=[]
 for sid,g in annual.groupby('series_id'):
  yy=sorted(set(g.year));records.append(dict(series_id=sid,audited_first_year=min(yy),audited_last_year=max(yy),audited_n_years=len(yy),audited_years=';'.join(map(str,yy)),audited_missing_years=';'.join(map(str,sorted(set(range(min(yy),max(yy)+1))-set(yy)))),audited_before_years=int((g.t<0).sum()),audited_post_years=int((g.t>0).sum())))
 latest=latest.merge(pd.DataFrame(records),on='series_id',validate='one_to_one');latest=latest.merge(zones[['series_id','constant_zones_every_station','years_no_shared_zone']].rename(columns={'constant_zones_every_station':'fish_station_shared_zones','years_no_shared_zone':'fish_years_no_all_station_shared_zone'}),on='series_id',how='left',validate='one_to_one')
 latest['latest_audit_date']='2026-11-02';latest['station_reference_evidence']='../all_program_station_refits_20261102/README.md';latest['model_scale_note']='Arithmetic role-mean ratios and log-station models differ; compare station-adjusted with unadjusted log-station model only'
 comp=comp.drop(columns=[x for x in latest if x!='series_id'],errors='ignore').merge(latest,on='series_id',how='left',validate='one_to_one');assert len(comp)==206 and comp.annual_reconstruction_verified.fillna(False).sum()==175
 assert not comp.loc[comp.exclude,'annual_reconstruction_verified'].fillna(False).any()
 biomass_path=Path('outputs/biomass_support_20261102/biomass_history_status.csv')
 if biomass_path.exists():
  bio=read(biomass_path);bio=bio[['series_id','support_status','singleton_pools','sparse_pools_lt10','shared_role_years']].rename(columns={'support_status':'latest_biomass_support','singleton_pools':'latest_singleton_size_pools','sparse_pools_lt10':'latest_small_size_pools','shared_role_years':'latest_shared_size_role_years'});comp=comp.drop(columns=[c for c in bio if c!='series_id'],errors='ignore').merge(bio,on='series_id',how='left',validate='one_to_one')
 comp['audit_scope']=np.where(comp.exclude,'Researcher/ecological exclusion retained; no new refit','Retained review candidate; not promoted to manuscript')
 # Nine domains: completed evidence is distinguished from open design decisions.
 status=[]
 for r in comp.itertuples():
  excluded=bool(r.exclude)
  status.append(dict(series_id=r.series_id,MPA=r.MPA,Source=r.Source,Taxa=r.Taxa,Resp=r.Resp,scope=r.audit_scope,
   identity_inclusion='Crosswalk and exclusions recorded; unassigned Trailer/Willows remain outside cohort',
   protection='Dates and current take class recorded; historical boundary/first-protection questions remain',
   time_coverage='Audited response-specific years updated' if not excluded else 'Historical coverage retained; excluded from current review',
   habitat_depth=('PISCO sheephead all-station zone screen complete; measured-depth equivalence unresolved' if r.Source=='PISCO' and r.Taxa=='S. pulcher' else 'Satellite spatial masks; diver depth not applicable' if r.Source=='Landsat' else 'Depth/habitat comparability remains unresolved'),
   survey_timing='Coverage recorded; no new seasonal cutoff selected',
   methods_effort=('CRANE/repeat/complete-pass rules implemented in review' if r.Source=='PISCO' and r.Taxa=='S. pulcher' else 'Empty-mask correction verified; quarterly support audit pending' if r.Source=='Landsat' else 'Source conventions screened; final visit/effort safeguards still pending'),
   biomass_support=(getattr(r,'latest_biomass_support',None) if r.Resp=='Bio' and pd.notna(getattr(r,'latest_biomass_support',None)) else 'Not a biomass response' if r.Resp!='Bio' else 'Satellite canopy; no animal/plant size conversion' if r.Source=='Landsat' else 'Existing conversion reproduced; expanded size-pool adequacy not certified'),
   station_references='Excluded; not newly fitted' if excluded else 'Station coverage and reference diagnostics complete; single-reference/absent persistent sets explicit',
   inference='No model admission while excluded' if excluded else 'Recovery diagnostics checked; before-after and final pooled/downstream integration pending',
   evidence='../all_program_station_refits_20261102/README.md' if not excluded else 'comparison_eligibility.csv'))
 status=pd.DataFrame(status);status.to_csv(O/'audit_status_by_comparison.csv',index=False);comp.to_csv(O/'comparison_eligibility.csv',index=False);comp.to_csv(O/'audit_master.csv',index=False)
 key=['MPA','Source','Taxa','Resp'];add=['audit_scope','audited_n_years','audited_before_years','audited_post_years','audit_n_reference_stations','audit_n_persistent_references','audit_post_reference_membership_changes','reference_omissions_checked','latest_audit_date','station_reference_evidence']
 prev=prev.drop(columns=add,errors='ignore').merge(comp[key+add],on=key,how='left',validate='one_to_one');prev.to_csv(O/'supplement_table_preview.csv',index=False)
 sites['latest_comparison_audit']='audit_master.csv; join by program/MPA/taxon, preserving response and method distinctions';sites.to_csv(O/'site_response_coverage.csv',index=False)
 ds['station_reference_audit']='../all_program_station_refits_20261102/README.md';ds['integration_status']='Review diagnostics consolidated; measurement/size-specific open items retained';ds.to_csv(O/'dataset_inventory.csv',index=False)
 issue.loc[issue.issue_id.eq('U06'),'finding']='Shared comparable zones approved in principle; station-level fish check finds Abalone2010/2014 have no all-station shared zone. Timing/depth restriction choices still pending. All-program reference influence checks completed.'
 issue.loc[issue.issue_id.eq('U06'),'evidence']='outputs/all_program_station_refits_20261102/README.md; outputs/station_coverage_refits_20261101/station_zone_summary.csv'
 issue.loc[issue.issue_id.eq('U08'),'next_action']='Extend size-donor adequacy and shared uncertainty to expanded roster; station reconstruction is not certification of size support'
 issue.to_csv(O/'unresolved_issues.csv',index=False)
 for name in ['dataset_audit','additional_analysis_checks']:
  d=pd.read_csv(O/(name+'.csv'))
  if name=='dataset_audit':
   d['latest_station_reference_check']='Completed2026-11-02 for retained response histories; see audit_master.csv. Does not resolve listed measurement/size gaps.'
  else:
   d['latest_status']=d.check_id.map({'Q01':'PISCO fish corrected; other visit/effort checks remain','Q02':'All-program station/reference diagnostics complete; primary weighting and station-zone standardization pending','Q03':'Annual biomass reconstruction verified; expanded donor adequacy remains open','Q04':'Detection limitations remain; no detection model fitted','Q05':'Existing offsets reproduced and held fixed; offset sensitivity remains open','Q06':'Polygon reference influence checked; quarterly timing and pixel support remain open'})
  d.to_csv(O/(name+'.csv'),index=False)
 freeze=pd.read_csv(O/'freeze_readiness.csv');freeze['latest_status']=freeze.item_id.map({'F1':'Shared-zone principle approved; timing and constant-station/depth restrictions undecided','F2':'Recovery diagnostics complete; final model admission remains separate','F3':'Current biomass reproduces; expanded size support pending','F4':'Shared-reference influence quantified; final joint uncertainty specification pending','F5':'950 station diagnostics and447 omissions completed; not final primary refit','F6':'Not completed; paper and downstream results unchanged','F7':'PISCO fish rule correction complete; other measurement checks remain'}).fillna('Existing requirement/limitation retained; not newly resolved');freeze.to_csv(O/'freeze_readiness.csv',index=False)
 views.update({'Comparisons':prev,'Audit status':status,'Station checks':comp[['MPA','Source','Taxa','Resp']+list(latest.columns[1:])],'Dataset audit':pd.read_csv(O/'dataset_audit.csv'),'Additional checks':pd.read_csv(O/'additional_analysis_checks.csv'),'Freeze readiness':freeze,'Questions':issue})
 pd.DataFrame([dict(decision='Shared comparable depth zones',status='Approved in principle; constant station/zone losses require review'),dict(decision='Timing',status='Show actual coverage before choosing cutoff'),dict(decision='Biomass',status='Keep supported density; require justified size support for biomass'),dict(decision='Reference composition',status='Compare fixed and available sets before choosing'),dict(decision='Incomplete reference documentation',status='Retain provisionally; test influence; preserve researcher exclusions'),dict(decision='Paper integration',status='On hold pending review of consolidated analyses')]).to_csv(O/'user_decisions.csv',index=False)
 if biomass_path.exists():
  views['Biomass support']=read(biomass_path)
  issue.loc[issue.issue_id.eq('U08'),'finding']='Expanded39 animal histories and38 direct/kelp histories screened:23 singleton pools,571 shared positive role-years; no missing positive animal pool. Donor representativeness and uncertainty remain open.'
  issue.loc[issue.issue_id.eq('U08'),'evidence']='outputs/biomass_support_20261102/README.md'
  issue.to_csv(O/'unresolved_issues.csv',index=False)
 resolution=Path('outputs/biomass_resolution_20261102')
 if (resolution/'verification.json').exists():
  loss=read(resolution/'donor_coverage_losses.csv');loss=loss[loss.scenario.eq('same_count_site')][['series_id','n_original','n_supported','lost_years']].rename(columns={'n_original':'size_check_original_post_years','n_supported':'size_check_count_site_supported_years','lost_years':'size_check_unsupported_years'})
  comp=comp.drop(columns=[c for c in loss if c!='series_id'],errors='ignore').merge(loss,on='series_id',how='left',validate='one_to_one')
  comp['biomass_resolution_evidence']=np.where(comp.Resp.eq('Bio'),'../biomass_resolution_20261102/README.md','')
  comp.to_csv(O/'audit_master.csv',index=False);comp.to_csv(O/'comparison_eligibility.csv',index=False)

  views['Biomass refits']=read(resolution/'separated_changes.csv')
  views['Red urchins']=read(resolution/'red_urchin_decisions.csv')
  issue.loc[issue.issue_id.eq('U08'),'finding']='35 pooled-size recovery histories checked with shared-record draws and singleton stress tests; 46 matched-donor refits verified. Missing-stipe paired year withheld in separate review input. Red-urchin donor representativeness, count/size thresholds and KFM dates remain unresolved.'
  issue.loc[issue.issue_id.eq('U08'),'next_action']='Resolve count-size protocol compatibility and original KFM date crosswalk; complete joint count/size and KFM kelp uncertainty before freezing biomass. Missing-stipe guard is implemented upstream; full production rebuild pending.'
  issue.loc[issue.issue_id.eq('U08'),'evidence']='outputs/biomass_resolution_20261102/README.md'
  issue.to_csv(O/'unresolved_issues.csv',index=False)
  views['Questions']=issue
 protocol=Path('outputs/biomass_protocols_20261102')
 if (protocol/'verification.json').exists():
  views['Biomass assumptions']=read(protocol/'assumptions_and_changes.csv')
  views['KFM size dates']=read(protocol/'retained_history_date_support.csv')
  comp['biomass_site_identity_alert']=''
  sbi=comp.MPA.eq('Santa Barbara Island SMR') & comp.Source.eq('KFM') & comp.Resp.eq('Bio') & comp.Taxa.isin(['M. franciscanus','S. purpuratus'])
  comp.loc[sbi,'biomass_site_identity_alert']='2005–2011 count site14 uses sizes matching site37; prior physical same-site support claim superseded. Later years unresolved.'
  comp.loc[sbi,'latest_biomass_support']='Source-site mismatch confirmed; biomass provisional; partial same-observation refit available'
  hp=comp.MPA.eq('Harris Point SMR') & comp.Resp.eq('Bio') & comp.Taxa.isin(['M. franciscanus','S. purpuratus'])
  comp.loc[hp,'biomass_site_identity_alert']='Combined a-k-01 size label merges legacy1+21; label-only count-site support requires reinterpretation'
  comp.to_csv(O/'audit_master.csv',index=False);comp.to_csv(O/'comparison_eligibility.csv',index=False)
  views['Size identity alerts']=comp.loc[comp.biomass_site_identity_alert.ne(''),['MPA','Source','Taxa','Resp','biomass_site_identity_alert']]
  issue.loc[issue.issue_id.eq('U08'),'finding']='KFM species/units verified;58369 frequency-record dates recovered as overlay. SBI2005–2011 counts match legacy14 but sizes match37. Prior same-site support superseded; isolated partial refit materially changes red-urchin slope.'
  issue.loc[issue.issue_id.eq('U08'),'next_action']='Resolve recent donor site/date identities; recompute physical coverage; complete joint uncertainty. No primary promotion of partial SBI refit.'
  issue.loc[issue.issue_id.eq('U08'),'evidence']='outputs/biomass_protocols_20261102/README.md; docs/BIOMASS_ASSUMPTIONS_AND_CHANGES.csv'
  issue.to_csv(O/'unresolved_issues.csv',index=False);views['Questions']=issue
 ex=Path('outputs/biomass_exclusion_reconciliation_20261102')
 if (ex/'verification.json').exists():
  affected=read(ex/'affected_histories.csv')
  views['Donor exclusion checks']=affected
  views['Donor exclusion refits']=read(ex/'effects.csv')
  comp['donor_exclusion_reconciliation']=''
  comp.loc[comp.series_id.isin(affected.series_id),'donor_exclusion_reconciliation']='Intended excluded donors survive alternate labels; historical and2015source confirmation. Partial<=2011refits only; recent identities unresolved.'
  comp.to_csv(O/'audit_master.csv',index=False);comp.to_csv(O/'comparison_eligibility.csv',index=False)
  issue.loc[issue.issue_id.eq('U08'),'finding']='Six biomass histories affected by excluded-donor alias leakage. Original code excludes sites21and37; legacy and2015source distributions identify them under labels1and14. Partial intended-donor refits verified24fits.'
  issue.loc[issue.issue_id.eq('U08'),'next_action']='Obtain/recover recent original-site sizes and upstream harmonization code; rebuild all affected donor pools before primary freeze.'
  issue.loc[issue.issue_id.eq('U08'),'evidence']='outputs/biomass_exclusion_reconciliation_20261102/README.md; docs/BIOMASS_ASSUMPTIONS_AND_CHANGES.csv'
  issue.to_csv(O/'unresolved_issues.csv',index=False);views['Questions']=issue
 joint=Path('outputs/biomass_joint_checks_20261102')
 if (joint/'verification.json').exists():
  views['Count-size matching']=read(joint/'count_size_history_summary.csv')
  views['Joint biomass draws']=read(joint/'conditional_joint_results.csv')
  views['Biomass draw holds']=read(joint/'resampling_admission.csv')
  admission=read(joint/'resampling_admission.csv').set_index('series_id')
  comp['joint_sampling_admitted']=comp.series_id.map(admission.admitted)
  comp['joint_sampling_scope']=comp.series_id.map(admission.reason)
  comp['joint_sampling_evidence']=comp.series_id.map(lambda x:'../biomass_joint_checks_20261102/README.md' if x in admission.index else '')
  comp.to_csv(O/'audit_master.csv',index=False);comp.to_csv(O/'comparison_eligibility.csv',index=False)
  issue.loc[issue.issue_id.eq('U08'),'finding']='Six alias-affected histories pending original sizes. Conditional joint count-size draws now cover34histories;23PISCO size/count zone flags. Total uncertainty and donor representativeness remain open.'
  issue.loc[issue.issue_id.eq('U08'),'next_action']='Resolve recent original donor identities; assess mismatched zones/borrowed sizes; coordinate full uncertainty before primary freeze.'
  issue.to_csv(O/'unresolved_issues.csv',index=False);views['Questions']=issue
 zr=Path('outputs/pisco_size_zone_resolution_20261102')
 if (zr/'summary.csv').exists():
  zs=read(zr/'summary.csv');views['Size-zone resolution']=zs
  views['Size-zone source dates']=read(zr/'source_zone_details.csv')
  comp['size_zone_resolution']=comp.series_id.map(zs.set_index('series_id').source_resolution)
  comp.to_csv(O/'audit_master.csv',index=False);comp.to_csv(O/'comparison_eligibility.csv',index=False)
  issue.loc[issue.issue_id.eq('U08'),'finding']='Six alias biomass histories unresolved. All23PISCO zone flags explained by counts filtered toMay–October but sizes retained from other months;66matched-year sensitivity fits checked.'
  issue.loc[issue.issue_id.eq('U08'),'next_action']='Resolve original KFM donors; choose survey timing policy before matching counts/sizes consistently; complete total uncertainty.'
  issue.to_csv(O/'unresolved_issues.csv',index=False);views['Questions']=issue
 am=Path('outputs/pisco_all_months_20261102')
 if (am/'summary.csv').exists():
  ams=read(am/'summary.csv');views['All-month rerun']=ams
  views['All-month timing']=read(am/'timing_summary.csv')
  comp['all_month_review_supported_years']=comp.series_id.map(ams.set_index('series_id').supported_post_years)
  comp['all_month_review_status']=comp.series_id.map(ams.set_index('series_id').status)
  comp.to_csv(O/'audit_master.csv',index=False);comp.to_csv(O/'comparison_eligibility.csv',index=False)
  issue.loc[issue.issue_id.eq('U08'),'next_action']='All-month six-history rerun retains63density years but25locally supported biomass years; assess borrowing/long size-date gaps and total uncertainty; KFM originals still pending.'
  issue.to_csv(O/'unresolved_issues.csv',index=False);views['Questions']=issue
 td=Path('outputs/pisco_timing_donors_20261102')
 if (td/'verification.json').exists():
  views['Survey date matching']=read(td/'timing_coverage.csv')
  views['Size donor validation']=read(td/'donor_validation_summary.csv')
  issue.loc[issue.issue_id.eq('U08'),'next_action']='624held-out size samples tested; same-role donor availability120/153at90days is not validated imputation. Address prediction uncertainty and Abalone protected-site validation gap; KFM originals pending.'
  issue.to_csv(O/'unresolved_issues.csv',index=False);views['Questions']=issue
 br=Path('outputs/pisco_biomass_borrowing_20261102')
 if (br/'summary.csv').exists():
  bsum=read(br/'summary.csv');views['Borrowing uncertainty']=bsum
  views['Borrowing validation']=read(br/'interval_validation_summary.csv')
  comp['conditional_borrowing_post_years']=comp.series_id.map(bsum.set_index('series_id').borrowing_post_years)
  comp.to_csv(O/'audit_master.csv',index=False);comp.to_csv(O/'comparison_eligibility.csv',index=False)
  issue.loc[issue.issue_id.eq('U08'),'next_action']='Count-led review keeps63density years; conditional borrowing sensitivity supports44biomass years. Protected-site calibration/38unsupported visits and total uncertainty remain open; KFM originals pending.'
  issue.to_csv(O/'unresolved_issues.csv',index=False);views['Questions']=issue
 views['Decisions']=pd.read_csv(O/'user_decisions.csv')
 (O/'consolidation_input_hashes.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in used},indent=2))
 from consolidate_current_decisions import apply as consolidate
 return consolidate(views)
