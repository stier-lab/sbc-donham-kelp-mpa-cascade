from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,re
O=Path('outputs/seasonal_zone_audit_20261101');D=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data')
a=pd.read_csv('outputs/researcher_scope_20261101/annual_input.csv');mm=pd.read_csv(D/'MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv').rename(columns={'site_id':'site','CA_MPA_Name_Short':'MPA'})
# Match known retained site/month membership, not every station in a broad MPA map.
known=pd.read_csv('tables/table_s_seasonal_site_months.csv');known=known[known.Source.isin(['KFM','LTER'])][['Source','Taxa','MPA','site','year']].drop_duplicates()
# Extend the old seasonal inventory to the current roster using source mappings.
text=(D.parent/'code/R/01_utils.R').read_text()
exc=[]
for name in ['EXCLUDED_KFM_SITES','EXCLUDED_LTER_SITES','EXCLUDED_LTER_SITE_CODES']:
 exc+=re.findall(r'"([^"]+)"',re.search(name+r' <- c\((.*?)\)',text,re.S).group(1))
mm['Source']=mm.data_source.str.upper();fishmm=mm[mm.status.isin(['mpa','reference'])].copy();mm=mm[~mm.site.isin(exc)&mm.status.isin(['mpa','reference'])]
known=mm[['Source','MPA','site']].drop_duplicates().merge(a[a.Source.isin(['KFM','LTER'])][['Source','Taxa','MPA','year']].drop_duplicates(),on=['Source','MPA'])
known_fish=fishmm[['Source','MPA','site']].drop_duplicates().merge(a[(a.Source=='KFM')&(a.Taxa=='S. pulcher')][['Source','Taxa','MPA','year']].drop_duplicates(),on=['Source','MPA'])
files=['MBON/SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv','MBON/SBCMBON_kelp_forest_integrated_fish_20231022.csv'];parts=[];inputs=[]
names={'Macrocystis pyrifera':'M. pyrifera','Panulirus interruptus':'P. interruptus','Semicossyphus pulcher':'S. pulcher','Bodianus pulcher':'S. pulcher','Strongylocentrotus purpuratus':'S. purpuratus','Mesocentrotus franciscanus':'M. franciscanus','Strongylocentrotus franciscanus':'M. franciscanus'}
for file in files:
 p=D/file;inputs.append(p);x=pd.read_csv(p,low_memory=False);x['Source']=x.data_source.str.upper();x['Taxa']=x.taxon_name.map(names);x['site']=x.site_id;x['year']=pd.to_datetime(x.date).dt.year;x['month']=pd.to_datetime(x.date).dt.month
 if 'integrated_fish' in file:x=x[(x.Source=='KFM')&(x.sample_method=='visualfish')&(x.Taxa=='S. pulcher')]
 else:x=x[((x.Source=='KFM')&~((x.Taxa=='M. pyrifera')&(x.proj_taxon_id=='t-k-002')))|((x.Source=='LTER')&(x.sample_method=='quad')&x.Taxa.isin(['S. purpuratus','M. franciscanus']))]
 x=x.merge(known_fish if 'integrated_fish' in file else known,on=['Source','Taxa','site','year']);x['unit']=x.subsite_id.astype(str)+'|'+x.transect_id.astype(str);x['file']=file;parts.append(x[['Source','Taxa','MPA','site','year','month','unit','sample_method','file']].drop_duplicates())
for file,tx in [('Annual_Kelp_All_Years_20240305.csv','M. pyrifera'),('Annual_fish_comb_20240307.csv','S. pulcher'),('Lobster_Abundance_All_Years_20230831.csv','P. interruptus')]:
 p=D/'LTER'/file;inputs.append(p);x=pd.read_csv(p,low_memory=False);x=x.rename(columns={'SITE':'site','YEAR':'year','MONTH':'month'});x['Source']='LTER';x['Taxa']=tx
 if tx=='M. pyrifera':x=x[x.FRONDS!=-99999]
 if tx=='S. pulcher':x=x[(x.SCIENTIFIC_NAME=='Semicossyphus pulcher')&(x.COUNT!=-99999)]
 x=x.merge(known,on=['Source','Taxa','site','year']);x['unit']=x.TRANSECT.astype(str);x['sample_method']='annual_'+tx;x['file']=file;parts.append(x[['Source','Taxa','MPA','site','year','month','unit','sample_method','file']].drop_duplicates())
e=pd.concat(parts);e.to_csv(O/'other_program_survey_units.csv',index=False)
rows=[]
for sid,h in a[a.Source.isin(['KFM','LTER'])].groupby('series_id'):
 r=h.iloc[0];z=e[(e.Source==r.Source)&(e.Taxa==r.Taxa)&(e.MPA==r.MPA)&e.year.isin(h.year)];foot=[set(q.site+'|'+q.unit) for _,q in z.groupby('year')];sites=[set(q.site) for _,q in z.groupby('year')]
 rows.append(dict(series_id=sid,Source=r.Source,MPA=r.MPA,Taxa=r.Taxa,Resp=r.Resp,paired_years=len(h),audited_years=z.year.nunique(),sites=z.site.nunique(),station_set_changes=len(set(tuple(sorted(s)) for s in sites))>1,survey_unit_set_changes=len(set(tuple(sorted(s)) for s in foot))>1,common_units=len(set.intersection(*foot)) if foot else 0,all_units=len(set.union(*foot)) if foot else 0,months=';'.join(map(str,sorted(z.month.unique()))),depth_available=False,scope='Survey-footprint audit; source date metadata only for lobster, not current corrected lengths'))
pd.DataFrame(rows).to_csv(O/'other_program_summary.csv',index=False)
(O/'other_program_hashes.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},indent=2));print(pd.DataFrame(rows).groupby('Source')[['paired_years','audited_years','station_set_changes','survey_unit_set_changes']].sum())
