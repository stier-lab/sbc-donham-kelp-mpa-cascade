"""Read-only source inventory. No abundance estimation, exclusions or model changes."""
from pathlib import Path
import pandas as pd, numpy as np, hashlib,json
D=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data')
O=Path('outputs/dataset_inventory_20261101');O.mkdir(exist_ok=True)
T={'Macrocystis pyrifera':'M. pyrifera','Panulirus interruptus':'P. interruptus','Semicossyphus pulcher':'S. pulcher','Bodianus pulcher':'S. pulcher','Strongylocentrotus purpuratus':'S. purpuratus','Mesocentrotus franciscanus':'M. franciscanus','Strongylocentrotus franciscanus':'M. franciscanus'}
C={'MACPYRAD':'M. pyrifera','MESFRAAD':'M. franciscanus','STRPURAD':'S. purpuratus','MAPY':'M. pyrifera','PINT':'P. interruptus','PANINT':'P. interruptus','SPUL':'S. pulcher','SEMPUL':'S. pulcher','STRPUR':'S. purpuratus','SPUR':'S. purpuratus','MESFRA':'M. franciscanus','SFRAN':'M. franciscanus','STRFRA':'M. franciscanus'}
S=[('pisco_fish','PISCO/MLPA_kelpforest_fish_2024.csv','pisco'),('pisco_swath','PISCO/MLPA_kelpforest_swath_2024.csv','pisco'),('mbon_fish','MBON/SBCMBON_kelp_forest_integrated_fish_20231022.csv','mbon'),('mbon_benthic','MBON/SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv','mbon'),('lter_fish','LTER/Annual_fish_comb_20240307.csv','lter'),('lter_kelp','LTER/Annual_Kelp_All_Years_20240305.csv','lter'),('lter_lobster','LTER/Lobster_Abundance_All_Years_20230831.csv','lter'),('invertebrate_sizes','ALL_sizefreq_2024.csv','size'),('kfm_kelp_sizes','MBON/KFM_Macrocystis_RawData_1984-2023.csv','kfmkelp')]
parts=[];inv=[];hashes={};codes=[];yearrows=[]
for did,rel,kind in S:
 p=D/rel;hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();head=list(pd.read_csv(p,nrows=0).columns)
 n=0;frames=[]
 for x in pd.read_csv(p,chunksize=200000,low_memory=False,encoding='latin1'):
  n+=len(x)
  if kind in ['pisco','size']:
   codes.extend([(did,str(k),int(v)) for k,v in x.classcode.value_counts().items()])
   # Survey opportunity retained independently of whether a focal organism was recorded.
   x=x.copy();x['program']='PISCO' if kind=='pisco' else 'Shared_sizes';x['site_id']=x.site.astype(str);x['taxon']=x.classcode.map(C).fillna('Other / survey opportunity');x['method_key']=x['method'].fillna('Unknown')
   x['date_key']=pd.to_datetime(x[['year','month','day']],errors='coerce');x['count_value']=pd.to_numeric(x['count'],errors='coerce');x['size_value']=pd.to_numeric(x['fish_tl'] if did=='pisco_fish' else x['size'],errors='coerce');x['depth_value']=pd.to_numeric(x.depth,errors='coerce');x['campaign']=x.survey_year.astype(str);x['zone_key']=x.zone.fillna('Unknown')
   x['unit_key']=x['transect'].astype(str)+(('|'+x.level.astype(str)) if 'level' in x else '')
  elif kind=='mbon':
   x=x[x.taxon_name.isin(T)].copy();x['program']=x.data_source.str.upper();x['site_id']=x.site_id.astype(str);x['taxon']=x.taxon_name.map(T);x['method_key']=x.sample_method.fillna('Unknown');x['date_key']=pd.to_datetime(x.date,errors='coerce');x['count_value']=pd.to_numeric(x['count'],errors='coerce');x['size_value']=np.nan;x['depth_value']=np.nan;x['zone_key']='Not supplied';x['campaign']='Not supplied';x['unit_key']=x.subsite_id.astype(str)+'|'+x.transect_id.astype(str)+'|'+x.replicate_id.astype(str)
  elif kind=='lter':
   if did=='lter_fish':x=x[x.SCIENTIFIC_NAME.isin(['Semicossyphus pulcher','Bodianus pulcher'])].copy()
   x=x.copy();x['program']='LTER';x['site_id']=x.SITE.astype(str);x['taxon']={'lter_fish':'S. pulcher','lter_kelp':'M. pyrifera','lter_lobster':'P. interruptus'}[did];x['method_key']=did;x['date_key']=pd.to_datetime(x.DATE,errors='coerce');x['count_value']=pd.to_numeric(x['FRONDS'] if did=='lter_kelp' else x.COUNT,errors='coerce');x['size_value']=pd.to_numeric(x['SIZE_MM'] if did=='lter_lobster' else x['SIZE'] if did=='lter_fish' else x.HLD_DIAM,errors='coerce');x['depth_value']=np.nan;x['zone_key']='Not supplied';x['campaign']='Not supplied';x['unit_key']=x.TRANSECT.astype(str)
  else:
   x=x.copy();x['program']='KFM';x['site_id']=x.SiteNumber.astype(str);x['taxon']='M. pyrifera';x['method_key']='plant stipe/diameter';x['date_key']=pd.to_datetime(x.SurveyDate,errors='coerce');x['count_value']=pd.to_numeric(x.Stipe_Count,errors='coerce');x['size_value']=pd.to_numeric(x.Diameter_cm,errors='coerce');x['depth_value']=np.nan;x['zone_key']='Not supplied';x['campaign']=x.SurveyYear.astype(str);x['unit_key']=x.Marker.astype(str)
  x['unit_key']=x['unit_key'].fillna('Unknown').astype(str);x['site_id']=x['site_id'].fillna('Unknown').astype(str)
  x['year_key']=x.date_key.dt.year;x['month_key']=x.date_key.dt.month
  if kind in ['pisco','size']:
   x['year_key']=x['year_key'].fillna(pd.to_numeric(x.year,errors='coerce'));x['month_key']=x['month_key'].fillna(pd.to_numeric(x.month,errors='coerce').where(pd.to_numeric(x.month,errors='coerce').between(1,12)))
  keep=['program','site_id','taxon','method_key','date_key','year_key','month_key','count_value','size_value','depth_value','zone_key','campaign','unit_key']
  frames.append(x[keep])
 z=pd.concat(frames,ignore_index=True);z['dataset_id']=did
 inv.append(dict(dataset_id=did,source_file=rel,source_sha256=hashes[str(p)],file_bytes=p.stat().st_size,source_rows=n,scanned_rows=len(z),programs=';'.join(sorted(z.program.unique())),first_date=str(z.date_key.min().date()),last_date=str(z.date_key.max().date()),raw_sites=z.site_id.nunique(),invalid_dates=int(z.date_key.isna().sum()),count_missing=int(z.count_value.isna().sum()),count_negative=int(z.count_value.lt(0).sum()),zero_count_rows=int(z.count_value.eq(0).sum()),depth_field_present='depth' in head,positive_depth_rows=int(z.depth_value.gt(0).sum()),size_field_present=kind!='mbon',positive_count_missing_nonpositive_size=int(((z.count_value>0)&(~(z.size_value>0))).sum()) if kind!='mbon' else None,selected_field_duplicate_rows=int(z.duplicated(keep).sum()),duplicate_interpretation='Screen only: selected fields omit observer/taxon-code details; no rows removed',scope='All source rows including nonfocal survey opportunity' if kind in ['pisco','size'] else 'Focal taxon records; not an all-species census',release_label=p.stem))
 # Annual sampling metadata distinguishes changing footprints from long overall coverage.
 for key,g in z.dropna(subset=['year_key']).groupby(['program','site_id','taxon','method_key','year_key']):
  dp=g.depth_value[g.depth_value>0];dd=sorted(set(g.date_key.dropna().dt.strftime('%Y-%m-%d')))
  yearrows.append(dict(dataset_id=did,program=key[0],site_id=key[1],taxon=key[2],method=key[3],year=int(key[4]),dates=';'.join(dd),months=';'.join(map(str,sorted(set(g.month_key.dropna().astype(int))))),zones=';'.join(sorted(set(g.zone_key))),unit_ids=';'.join(sorted(set(g.unit_key))),n_units=g.unit_key.nunique(),depth_min_m=dp.min(),depth_max_m=dp.max(),depth_record_median_m=dp.median(),positive_depth_records=len(dp),rows=len(g),missing_counts=int(g.count_value.isna().sum()),negative_counts=int(g.count_value.lt(0).sum()),zero_rows=int(g.count_value.eq(0).sum()),incomplete_full_dates=int(g.date_key.isna().sum())))
 # Aggregated station/response/method coverage; exact date strings are metadata, never individual animal records.
 for key,g in z.groupby(['program','site_id','taxon','method_key'],dropna=False):
  yrs=sorted(set(g.year_key.dropna().astype(int)));depth=g.depth_value[g.depth_value>0];ds=sorted(set(g.date_key.dropna().dt.strftime('%Y-%m-%d')))
  parts.append(dict(dataset_id=did,program=key[0],site_id=key[1],taxon=key[2],method=key[3],observed_years=';'.join(map(str,yrs)),survey_dates=';'.join(ds),months=';'.join(map(str,sorted(set(g.month_key.dropna().astype(int))))),raw_record_rows=len(g),zero_count_rows=int(g.count_value.eq(0).sum()),missing_count_rows=int(g.count_value.isna().sum()),negative_count_rows=int(g.count_value.lt(0).sum()),positive_count_missing_nonpositive_size=int(((g.count_value>0)&(~(g.size_value>0))).sum()) if kind!='mbon' else None,depth_min_m=depth.min(),depth_max_m=depth.max(),positive_depth_records=len(depth),zones=';'.join(sorted(set(g.zone_key))),survey_units=g.unit_key.nunique(),count_unit='fronds/stipes' if 'kelp' in did else 'organisms / recorded group frequency',size_status='Separate size source required' if kind=='mbon' else 'Recorded size; validity and population match require response-specific checks',coverage_basis='Recorded survey opportunity; not necessarily valid focal measurement' if key[2]=='Other / survey opportunity' else 'Recorded focal rows, including zeros and missing measurements; absence of a row is not coded zero'))
 print(did,n,len(z),flush=True)
pd.DataFrame(inv).to_csv(O/'raw_dataset_scan.csv',index=False);pd.DataFrame(parts).to_csv(O/'raw_station_response_scan.csv',index=False)
c=pd.DataFrame(codes,columns=['dataset_id','classcode','rows']).groupby(['dataset_id','classcode'],as_index=False).rows.sum();c.to_csv(O/'taxon_code_inventory.csv',index=False)
(O/'raw_source_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')

pd.DataFrame(yearrows).to_csv(O/'raw_station_year_audit.csv',index=False)
