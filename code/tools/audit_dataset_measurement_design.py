"""Read-only measurement/design screen; raw scope, never automatic exclusions."""
from pathlib import Path
import pandas as pd, numpy as np, hashlib, json
O=Path('outputs/dataset_inventory_20261101');D=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data')
inv=pd.read_csv(O/'dataset_inventory.csv');metrics=[];protocol=[];visits=[];hashes={}
T=['Macrocystis pyrifera','Panulirus interruptus','Semicossyphus pulcher','Bodianus pulcher','Strongylocentrotus purpuratus','Mesocentrotus franciscanus','Strongylocentrotus franciscanus']
def num(x,c):return pd.to_numeric(x[c],errors='coerce') if c in x else pd.Series(np.nan,index=x.index)
for r in inv.itertuples():
 p=D/r.source_file;hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();frames=[];total=0
 for x in pd.read_csv(p,encoding='latin1',chunksize=200000,low_memory=False):
  total+=len(x)
  if r.dataset_id.startswith('mbon'):x=x[x.taxon_name.isin(T)]
  elif r.dataset_id=='lter_fish':x=x[x.SCIENTIFIC_NAME.isin(T)]
  elif r.dataset_id=='pisco_fish':x=x[x.classcode.eq('SPUL')]
  elif r.dataset_id=='pisco_swath':x=x[x.classcode.isin(['MACPYRAD','PANINT','STRPURAD','MESFRAAD'])]
  frames.append(x)
 x=pd.concat(frames,ignore_index=True);did=r.dataset_id
 if did=='landsat_annual':
  yearcols=[c for c in x if str(c).isdigit() and 1900<=int(c)<=2100]
  v=x[yearcols].apply(pd.to_numeric,errors='coerce')
  metrics.append(dict(dataset_id=did,source_rows=total,screened_rows=len(x),screen_scope='Harmonized polygon annual rows; complete supplier delivery indexed separately',missing_measurements=int(v.isna().sum().sum()),negative_measurements=int((v<0).sum().sum()),explicit_zero_measurements=int(v.eq(0).sum().sum()),exact_repeat_rows=int(x.duplicated().sum()),measurement_columns=';'.join(yearcols)))
  continue
 if did.startswith('pisco') or did=='invertebrate_sizes':
  site='site';code='classcode';count='count';area=None;method='method';program='campus';size='fish_tl' if did=='pisco_fish' else 'size';unit=['site','zone','transect']+(['level'] if did=='pisco_fish' else [])
  dates=pd.to_datetime(x[['year','month','day']],errors='coerce');years=num(x,'year')
 elif did.startswith('mbon'):
  site='site_id';code='proj_taxon_id';count='count';area='area';method='sample_method';program='data_source';size=None;unit=['data_source','sample_method','site_id','subsite_id','transect_id','replicate_id'];dates=pd.to_datetime(x.date,errors='coerce');years=dates.dt.year
 elif did.startswith('lter'):
  site='SITE';code='SP_CODE' if 'SP_CODE' in x else None;count='FRONDS' if did=='lter_kelp' else 'COUNT';area='AREA';method=None;program=None;size={'lter_fish':'SIZE','lter_lobster':'SIZE_MM','lter_kelp':'HLD_DIAM'}[did];unit=[c for c in ['SITE','TRANSECT','QUAD','SIDE','REPLICATE'] if c in x];dates=pd.to_datetime(x.DATE,errors='coerce');years=num(x,'YEAR')
 else:
  site='SiteNumber';code='ScientificName';count='Stipe_Count';area=None;method=None;program=None;size='Diameter_cm';unit=['SiteNumber','Marker'];dates=pd.to_datetime(x.SurveyDate,errors='coerce');years=num(x,'SurveyYear')
 cnt=num(x,count);sz=num(x,size);ar=num(x,area)
 # Exact equality is a screen, not proof of accidental duplication (animal identifiers often absent).
 row=dict(dataset_id=did,source_rows=total,screened_rows=len(x),screen_scope='All delivered focal rows before geographic, method or eligibility filters' if did!='invertebrate_sizes' else 'All delivered size codes before eligibility filters',missing_measurements=int(cnt.isna().sum()),negative_measurements=int((cnt<0).sum()),explicit_zero_measurements=int(cnt.eq(0).sum()),exact_repeat_rows=int(x.duplicated().sum()),missing_complete_dates=int(dates.isna().sum()),calendar_year_disagreements=int((dates.notna()&years.notna()&dates.dt.year.ne(years)).sum()),missing_physical_key_rows=int(x[unit].isna().any(axis=1).sum()),physical_key_fields=';'.join(unit),area_field=area or 'not provided',area_values=';'.join(map(str,sorted(ar.dropna().unique()))) if area else '',nonpositive_or_missing_area=int((ar.isna()|ar.le(0)).sum()) if area else np.nan,positive_count_missing_or_nonpositive_size=int((cnt.gt(0)&(sz.isna()|sz.le(0))).sum()) if size else np.nan,size_screen_meaning={'pisco_swath':'Mixed meanings; use per-code rows. Only MACPYRAD needs size as stipes/plant here.','lter_kelp':'Holdfast diameter; not a required length for stipe-based conversion','kfm_kelp_sizes':'Holdfast diameter; verify conversion requirements'}.get(did,'Length screen, not an automatic exclusion' if size else 'No individual length field'))
 x=x.copy();x['_date']=dates;x['_year']=years;x['_code']=x[code].fillna('Unknown').astype(str) if code else did;x['_method']=x[method].fillna('Unknown').astype(str) if method else did;x['_program']=x[program].fillna('Unknown').astype(str) if program else ('KFM' if did=='kfm_kelp_sizes' else 'LTER');x['_cnt']=cnt;x['_size']=sz;x['_area']=ar
 # Multiple dates behind a reused unit label: an annual sum can pool effort. Not all such units enter the current paper.
 uu=x[unit+['_year','_date']].drop_duplicates();gg=uu.groupby(unit+['_year'],dropna=False)._date.nunique();multi=gg[gg>1].reset_index(name='distinct_dates');multi['dataset_id']=did;visits.append(multi)
 row['unit_years_with_multiple_dates']=len(multi);metrics.append(row)
 for key,g in x.groupby(['_program','_method','_code'],dropna=False):
  protocol.append(dict(dataset_id=did,program=key[0],method=key[1],taxon_code=key[2],rows=len(g),first_year=g._year.min(),last_year=g._year.max(),area_values=';'.join(map(str,sorted(g._area.dropna().unique()))),positive_count_rows=int(g._cnt.gt(0).sum()),missing_counts=int(g._cnt.isna().sum()),negative_counts=int(g._cnt.lt(0).sum()),positive_count_missing_or_nonpositive_size=int((g._cnt.gt(0)&(g._size.isna()|g._size.le(0))).sum()) if size else np.nan,missing_or_nonpositive_area=int((g._area.isna()|g._area.le(0)).sum()) if area else np.nan,lookedfor_values=';'.join(map(str,g.LOOKEDFOR.dropna().unique())) if 'LOOKEDFOR' in g else '',level_values=';'.join(map(str,g.level.dropna().unique())) if 'level' in g else ''))
 print(did,len(x),'repeat rows',row['exact_repeat_rows'],'multiple-date units',len(multi),flush=True)
pd.DataFrame(metrics).to_csv(O/'measurement_design_scan.csv',index=False)
pd.DataFrame(protocol).to_csv(O/'measurement_protocols.csv',index=False)
pd.concat(visits,ignore_index=True).to_csv(O/'repeated_unit_dates.csv',index=False)
for source in ['04_pisco_processing.R','05_kfm_processing.R','06_lter_processing.R']:
 p=D.parent/'code/R'/source;hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
(O/'measurement_source_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
# Fish effort must come from all surveyed taxa, not just sheephead detections.
f=pd.read_csv(D/'PISCO/MLPA_kelpforest_fish_2024.csv',encoding='latin1',low_memory=False)
f=f[f.campus.isin(['UCSB','VRG']) & f.method.isin(['SBTL_FISH_PISCO','SBTL_FISH_VRG']) & f.LOOKEDFOR.eq('yes') & f.level.ne('CAN')].copy()
f['_date']=pd.to_datetime(f[['year','month','day']],errors='coerce')
ev=['campus','method','site','year','month','day','zone','transect']
e=f.groupby(ev,dropna=False).agg(levels=('level',lambda s:';'.join(sorted(s.dropna().unique()))),n_levels=('level','nunique')).reset_index()
q=e.groupby(['site','year'],dropna=False).agg(survey_units=('n_levels','size'),level_combinations=('levels',lambda s:' | '.join(sorted(s.unique()))),min_levels=('n_levels','min'),max_levels=('n_levels','max')).reset_index();q['unequal_level_coverage']=q.min_levels.ne(q.max_levels)
q.to_csv(O/'pisco_fish_level_coverage.csv',index=False)
rep=e.groupby(['campus','method','site','year','zone','transect'],dropna=False).size().reset_index(name='survey_dates');rep=rep[rep.survey_dates>1]
rep.to_csv(O/'pisco_fish_revisited_units.csv',index=False)
print('PISCO standard-fish site-years',len(q),'unequal levels',int(q.unequal_level_coverage.sum()),'revisited units',len(rep),flush=True)
