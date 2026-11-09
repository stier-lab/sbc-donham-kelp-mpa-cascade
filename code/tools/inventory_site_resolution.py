"""Apply reviewed site hierarchy to inventory only; never change analytic membership."""
from pathlib import Path
import pandas as pd,json,re,hashlib
ROOT=Path(__file__).resolve().parents[2];O=ROOT/'outputs/dataset_inventory_20261101';R=ROOT/'outputs/site_assignment_resolution_20261101'
def apply(sites):
 res=pd.read_csv(R/'assignment_resolution.csv');links=pd.read_csv(R/'constituent_station_evidence.csv')
 utils=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/code/R/01_utils.R');s=utils.read_text()
 def const(n):return set(re.findall(r'"([^"]+)"',re.search(n+r' <- c\((.*?)\)',s,re.S).group(1)))
 exref=const('EXCLUDED_REFERENCE_SITES');exmpa=const('EXCLUDED_MPAS');exk=const('EXCLUDED_KFM_SITES')
 sites=sites.copy();sites['inventory_record_type']=sites.dataset_id.map(lambda x:'Observed source/response record' if x!='unmatched' else 'Unresolved map entry')
 sites['group_names']='';sites['linked_constituent_stations']='';sites['resolved_group_scope']='';sites['resolved_group_scope_reason']=''
 sites['inherited_pisco_reference_exclusion']=(sites.program=='PISCO')&sites.site_id.isin(exref)
 sites['inherited_whole_mpa_exclusion']=sites.MPA.isin(exmpa)
 sites['inherited_kfm_nonfish_exclusion']=(sites.program=='KFM')&sites.site_id.isin(exk)
 sites['exclusion_interpretation']='Inherited rule flags, not a new filter; newer researcher scope and response-specific fish exceptions control analysis'
 cross=(sites.program=='PISCO')&sites.source_namespace.isin(['KFM','LTER'])&(sites.dataset_id=='unmatched');sites.loc[cross,'inventory_record_type']='Cross-program map entry'
 for r in res.itertuples():
  mpa='' if r.MPA=='Unassigned' else r.MPA
  hit=(sites.program=='PISCO')&(sites.site_id==r.assigned_group)&(sites.MPA.fillna('')==mpa);assert hit.sum()==1
  isgroup=r.n_children>0
  sites.loc[hit,'inventory_record_type']='Grouped-site label; not a sampled station' if isgroup else 'Unassigned map label; no records recovered'
  sites.loc[hit,'record_link_status']=r.identity_class;sites.loc[hit,'linked_constituent_stations']=r.raw_children if pd.notna(r.raw_children) else ''
  sites.loc[hit,'resolved_group_scope']=r.analysis_status;sites.loc[hit,'resolved_group_scope_reason']=r.reason
  for child in str(r.raw_children).split(';') if isgroup else []:
   ch=(sites.program=='PISCO')&(sites.site_id==child)&(sites.MPA.fillna('')==mpa)
   sites.loc[ch,'group_names']=r.assigned_group
  # Do not fill a parent's missing years with a child union: that would count a second sampled station.
  if isgroup:assert sites.loc[hit,'n_years'].eq(0).all()
 res.to_csv(O/'site_group_resolution.csv',index=False)
 physical=links[['MPA','group_name','raw_station','campus','explicit_reference_exclusion','whole_mpa_exclusion']].drop_duplicates();physical.to_csv(O/'site_group_links.csv',index=False)
 sites[sites.dataset_id=='unmatched'][['program','MPA','site_id','source_namespace','inventory_record_type','record_link_status','linked_constituent_stations','resolved_group_scope']].to_csv(O/'station_link_exceptions.csv',index=False)
 (O/'site_resolution_input_hashes.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [R/'assignment_resolution.csv',R/'constituent_station_evidence.csv',utils]},indent=2)+'\n')
 return sites
