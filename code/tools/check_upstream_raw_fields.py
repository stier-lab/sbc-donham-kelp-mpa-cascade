from pathlib import Path
import pandas as pd,numpy as np,os,json,hashlib
root=Path.cwd();dr=Path(os.environ.get('MPA_DATA_REPO',str(root.parent/'Donham-Stier-CA-MPA-Data-2026')));out=root/'outputs/upstream_sites_20261101';base=dr/'data'
targets=['Anacapa Island SMCA','Painted Cave SMCA','Point Dume SMCA','Point Dume SMR']
m=pd.read_csv(base/'PISCO/master_site_table_Emilyedit.csv',encoding='latin1')
m=m[m.CA_MPA_Name_Short.isin(targets)]
assert m.groupby(['site','CA_MPA_Name_Short']).site_status.nunique().max()==1
m.to_csv(out/'pisco_source_assignments.csv',index=False)
checks=[];hashes={}
k=pd.read_csv(base/'MBON/SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv',usecols=['site_id','taxon_name','sample_method','count','area'],dtype={'count':'string'})
k=k[k.site_id.isin(['a-k-11','a-k-31','a-k-32','a-k-34'])&k.taxon_name.isin(['Macrocystis pyrifera','Mesocentrotus franciscanus','Panulirus interruptus','Strongylocentrotus purpuratus'])]
k['count']=pd.to_numeric(k['count'],errors='raise')
rows=[]
for key,g in k.groupby(['site_id','taxon_name','sample_method']):
 rows.append(dict(zip(['site_id','taxon_name','sample_method'],key),records=len(g),missing_count=int(g['count'].isna().sum()),negative_count=int((g['count']<0).sum()),nonpositive_area=int((g.area<=0).sum())))
pd.DataFrame(rows).to_csv(out/'kfm_raw_field_checks.csv',index=False)
for f in ['PISCO/MLPA_kelpforest_swath_2024.csv','PISCO/MLPA_kelpforest_fish_2024.csv']:
 d=pd.read_csv(base/f,low_memory=False,encoding='latin1');d=d[d.site.isin(m.site)]
 checks.append(dict(file=f,rows=len(d),missing_count=int(d['count'].isna().sum()),negative_count=int((d['count']<0).sum()),first_year=int(d.year.min()),last_year=int(d.year.max())))
 if 'fish_tl' in d:
  bad=d[(d.classcode=='SPUL')&(d['count']>0)&((~np.isfinite(d.fish_tl))|(d.fish_tl<=0))&(d.level!='CAN')];bad.to_csv(out/'invalid_sheephead_lengths.csv',index=False)
 hashes[str(base/f)]=hashlib.sha256((base/f).read_bytes()).hexdigest()
pd.DataFrame(checks).to_csv(out/'raw_field_checks.csv',index=False)
for f in ['PISCO/master_site_table_Emilyedit.csv','MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv','MBON/SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv','MBON/SBCMBON_kelp_forest_integrated_fish_20231022.csv','MBON/KFM_Macrocystis_RawData_1984-2023.csv','ALL_sizefreq_2024.csv']:
 hashes[str(base/f)]=hashlib.sha256((base/f).read_bytes()).hexdigest()
(out/'raw_input_hashes.json').write_text(json.dumps(hashes,indent=2));print(checks);print('Invalid positive-count sheephead lengths:',len(bad))
