"""Resolve cross-program mapping labels and cover all24 delivered satellite polygons."""
from pathlib import Path
import pandas as pd,numpy as np,hashlib,json
O=Path('outputs/dataset_inventory_20261101');D=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data');used={}
def rd(p,**kw):
 p=Path(p);used[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();return pd.read_csv(p,**kw)
site=pd.read_csv(O/'site_response_coverage.csv',dtype={'site_id':str});raw=pd.read_csv(O/'raw_station_response_scan.csv',dtype={'site_id':str});master=rd(D/'PISCO/master_site_table_Emilyedit.csv')
ns=master.groupby('site').campus.agg(lambda v:';'.join(sorted(v.dropna().unique()))).to_dict();source_lookup=raw.groupby('site_id').program.agg(lambda v:';'.join(sorted(v.unique()))).to_dict()
site['source_namespace']=site.apply(lambda r:ns.get(r.site_id,'Unknown') if r.program=='PISCO' else r.program,axis=1)
site['other_source_programs_with_same_id']=site.site_id.map(source_lookup).fillna('No exact ID in scanned records')
site['record_link_status']=np.where(site.dataset_id!='unmatched','Exact source/program/station match',np.where(site.source_namespace.isin(['KFM','LTER'])&(site.program=='PISCO'),'Cross-program row in PISCO master map; not a missing PISCO survey','No exact source match; alias or unavailable historical record requires review'))
site.to_csv(O/'site_response_coverage.csv',index=False)
site[site.dataset_id=='unmatched'][['program','MPA','site_id','source_namespace','other_source_programs_with_same_id','record_link_status']].to_csv(O/'station_link_exceptions.csv',index=False)
# Exact legacy-name crosswalk evidence; never fuzzy-match or silently merge split stations.
import re
norm=lambda x:re.sub(r'[^A-Z0-9]','',str(x).upper()) if pd.notna(x) else ''
legacy=[]
for rel in ['PISCO/MLPA_kelpforest_fish_2024.csv','PISCO/MLPA_kelpforest_swath_2024.csv']:
 z=rd(D/rel,usecols=['site','site_name_old'],encoding='latin1').drop_duplicates();z['legacy_key']=z.site_name_old.map(norm);legacy.append(z)
legacy=pd.concat(legacy).drop_duplicates();ae=[]
for r in site[(site.dataset_id=='unmatched')&~site.record_link_status.str.startswith('Cross-program')].itertuples():
 m=master[master.site==r.site_id];keys={norm(r.site_id)}
 oldcol=next((c for c in master if c.strip()=='site_name_old'),None)
 if oldcol:keys|={norm(x) for x in m[oldcol].dropna()}
 keys.discard('');hit=legacy[legacy.legacy_key.isin(keys)];targets=sorted(hit.site.dropna().unique())
 ae.append(dict(MPA=r.MPA,assigned_site=r.site_id,exact_legacy_name_targets=';'.join(targets),n_targets=len(targets),resolution='No exact legacy-name evidence recovered' if not targets else 'Legacy-name evidence only; verify physical station/zone identity before joining'))
pd.DataFrame(ae).to_csv(O/'alias_evidence.csv',index=False)
# All delivered satellite polygons, including those excluded before harmonization.
labels=rd('data/external/tom_bell_2026-10-20/03_Site_List.csv');checks=rd('outputs/landsat_v2_verification/mask_array_checks.csv');known=rd('outputs/landsat_v2_verification/input_checks.csv');sel=rd('outputs/researcher_scope_20261101/selection_decisions.csv');names=known.groupby('site').MPA.first().to_dict();rows=[]
for suffix,grp in [('MPAs','MPA'),('MPAC1','Ref_1'),('MPAC2','Ref_2')]:
 d=rd('data/external/tom_bell_2026-10-21/rebuilt_'+suffix+'.csv')
 for _,r in d.iterrows():
  i=int(r['Site_#']);q=checks[(checks.site==i)&(checks.group==grp)];assert len(q)==1;q=q.iloc[0];label=str(labels.loc[labels['Site_#']==i,'Site'].iloc[0]).strip("'");name=names.get(i,label)
  yy=[int(c) for c in d if c.isdigit() and pd.notna(r[c])] if not q['empty'] else []
  rows.append(dict(site_number=i,MPA=name,source_label=label,group=grp,pixels=int(q.n_pixels),empty_mask=bool(q['empty']),n_years=len(yy),observed_years=';'.join(map(str,yy)),candidate_reconstruction_present=((sel.MPA==name)&(sel.Source=='Landsat')).any(),status='Unavailable: empty polygon; values not biological zeros' if q['empty'] else 'Nonempty saved mask; previous pixel reconstruction verified; cohort selection separate'))
x=pd.DataFrame(rows);assert len(x)==72 and x.empty_mask.sum()==14;x.to_csv(O/'satellite_full_delivery.csv',index=False)
# Provider survey-year differences are metadata facts, not new time coding.
cc=[]
for did,rel in [('pisco_fish','PISCO/MLPA_kelpforest_fish_2024.csv'),('pisco_swath','PISCO/MLPA_kelpforest_swath_2024.csv'),('invertebrate_sizes','ALL_sizefreq_2024.csv')]:
 d=rd(D/rel,usecols=['site','year','survey_year','month','day'],encoding='latin1');d['dataset_id']=did;g=d.groupby(['dataset_id','site','year','survey_year','month','day'],dropna=False).size().reset_index(name='source_rows');g['calendar_campaign_differ']=g.year!=g.survey_year;cc.append(g)
cc=pd.concat(cc);cc.to_csv(O/'survey_calendar_campaign_audit.csv',index=False)
(O/'additional_input_hashes.json').write_text(json.dumps(used,indent=2)+'\n')
print('Cross-program map rows',sum(site.record_link_status.str.startswith('Cross-program')),'full satellite polygons',x.site_number.nunique(),'nonempty arrays',sum(~x.empty_mask),'campaign/calendar differing metadata cells',sum(cc.calendar_campaign_differ))
