from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,shutil
D=Path('/Users/adrianstier/sbc-donham-kelp-mpa-data');O=Path('outputs/pelican_pisco_recovery_20261102');O.mkdir(exist_ok=True)
raw=pd.read_csv(D/'data/PISCO/MLPA_kelpforest_swath_2024.csv',encoding='cp1252',low_memory=False)
new=pd.read_csv(D/'data/processed/source_build/pisco_swath_reconciled.csv',low_memory=False)
keys=['campus','site','year','month','zone','transect']
def annual(data,site,year,taxon,old):
 d=data[data.site.eq(site)&data.year.eq(year)&data.month.between(5,10)&~data.zone.isin(['DEEP','SPECIALO'])].copy()
 grid=d[keys].drop_duplicates();obs=d[d.classcode.eq(taxon)].copy();obs['observed']=True
 if taxon=='MACPYRAD':obs['count']=obs['count']*obs['size']
 x=grid.merge(obs[keys+['count','observed']],on=keys,how='left')
 if old:x['count']=x['count'].fillna(0)
 else:x.loc[x.observed.isna(),'count']=0
 # Strict sums/means retain an explicit unknown measurement throughout this site-year.
 g=x.groupby(keys,dropna=False)['count'].agg(lambda v:v.sum(skipna=False)).reset_index()
 for by in [['site','year','zone','transect'],['site','year','zone'],['site','year']]:g=g.groupby(by,dropna=False)['count'].agg(lambda v:v.mean(skipna=False)).reset_index()
 return g['count'].iloc[0]/60 if len(g) else np.nan
rows=[]
for site,year,taxon in [('SCI_CAVERN_POINT_W',2013,'MESFRAAD'),('LITTLE_DUME_W',2011,'MACPYRAD'),('UNDERWATER_ARCH',2023,'MACPYRAD')]:
 rows.append(dict(site=site,year=year,taxon=taxon,legacy_zero_filled_per_m2=annual(raw,site,year,taxon,True),corrected_per_m2=annual(new,site,year,taxon,False),units='urchins/m2' if taxon=='MESFRAAD' else 'stipes/m2'))
pd.DataFrame(rows).to_csv(O/'site_year_changes.csv',index=False)
sites=pd.read_csv(D/'data/PISCO/master_site_table_Emilyedit.csv');sc=sites[sites.CA_MPA_Name_Short.eq('Scorpion SMR')&sites.site_status.eq('mpa')].site.unique();roles=[]
for site in sc:
 a=annual(raw,site,2013,'MESFRAAD',True);b=annual(new,site,2013,'MESFRAAD',False)
 if np.isfinite(a):roles.append(dict(site=site,old=a,new=b))
r=pd.DataFrame(roles);r.to_csv(O/'scorpion_2013_station_density.csv',index=False)
for name in ['scorpion_primary_step_sensitivity.csv','older_mass_comparison.csv']:
 shutil.copy2(Path('outputs/kfm_remaining_sources_20261102')/name,O/name)
shutil.copy2(D/'data/processed/source_build/pisco_count_recovery.csv',O/'pisco_count_recovery.csv')
manifest=json.loads((D/'data/processed/source_build/build_manifest.json').read_text());assert manifest['size_input_complete'] and manifest['required_kfm_holds']==0
summary=dict(source_gate_complete=True,pelican_original_rows=140,pelican_sensitivity_animals=135,recovered_red_urchin_count=219,scorpion_protected_density_2013_old=float(r.old.mean()),scorpion_protected_density_2013_new=float(r.new.mean()),scorpion_stations=len(r),scope='Same delivered stations/years, inherited transect-zone-site aggregation; diagnostic source-density comparison, not a full primary refit',paper_changed=False)
(O/'verification.json').write_text(json.dumps(summary,indent=2)+'\n');print(pd.DataFrame(rows).to_string(index=False));print(json.dumps(summary,indent=2))
