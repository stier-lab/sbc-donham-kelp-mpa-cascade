"""Review-only same-cohort model refit, current versus Emily archived size pool."""
from pathlib import Path
import os,shutil,subprocess,sys,json
import pandas as pd
R=Path(__file__).resolve().parents[2];M=R/'outputs/integrated_review_20261103';O=R/'outputs/emily_biomass_bridge_20261103'
arm=sys.argv[1];assert arm in ['current','archived'];name='emily_matched_'+arm;V=M/'sensitivities'/name;P=V/'private';P.mkdir(parents=True,exist_ok=True)
for n in ['comparison_year_roster.csv','station_assignments.csv','visit_biomass_support.csv','size_pool_inventory.csv','input_hashes.json']:shutil.copy2(M/n,V/n)
for n in ['measurement_uses.parquet','size_pool_members.parquet','lobster_imputation_members.csv']:shutil.copy2(M/'private'/n,P/n)
a=pd.read_csv(O/'paired_annual_values.csv');matched=set(map(tuple,a.loc[a.matched,['series_id','year']].itertuples(index=False,name=None)))
x=pd.read_parquet(M/'private/supported_measurement_uses.parquet');target=x.Source.eq('PISCO')&x.Resp.eq('Bio')&x.Taxa.isin(['M. franciscanus','S. purpuratus'])
keys=pd.MultiIndex.from_frame(x[['series_id','year']]);keep=keys.isin(matched);x.loc[target&~keep,'biomass']=float('nan')
if arm=='archived':
 pools=pd.read_csv(O/'emily_style_pools.csv');pools['pool_id']='archive_pool|'+pools.series_id+'|'+pools.status+'|'+pools.year.astype(str)
 x=x.merge(pools[['series_id','status','year','pool_id','legacy_mean_mass']],on=['series_id','status','year'],how='left',validate='many_to_one')
 m=target&keep;positive=m&x.density.gt(0);zero=m&x.density.eq(0);assert x.loc[positive,'legacy_mean_mass'].notna().all()
 x.loc[positive,'biomass']=x.loc[positive,'density']*x.loc[positive,'legacy_mean_mass'];x.loc[positive,'donor_pool']=x.loc[positive,'pool_id'].map(lambda s:json.dumps([s]));x.loc[positive,'size_basis']='emily_archive_pooled_size';x.loc[zero,'biomass']=0.;x.loc[zero,'donor_pool']='';x.loc[zero,'size_basis']='observed_zero';x=x.drop(columns=['pool_id','legacy_mean_mass'])
 inv=pd.read_csv(V/'size_pool_inventory.csv');added=pools.assign(Source='PISCO',site='archive_pooled_comparison',date='archive_pooled_year',zone='archive_pooled',mean_mass=pools.legacy_mean_mass,source_rows=pools.source_rows,singleton=False,zero_empirical_variance=False)[['Source','site','year','date','zone','Taxa','pool_id','mean_mass','total_mass','n_sized','source_rows','singleton','zero_empirical_variance']];pd.concat([inv,added],ignore_index=True).to_csv(V/'size_pool_inventory.csv',index=False)
 mem=pd.read_parquet(P/'size_pool_members.parquet');am=pd.read_csv(O/'archived_donor_members.csv');pd.concat([mem,am],ignore_index=True).to_parquet(P/'size_pool_members.parquet',index=False)
 v=pd.read_csv(V/'visit_biomass_support.csv');v.loc[v.series_id.isin(a.series_id.unique()),['prediction_error_sd','prediction_shared_sd']]=0.;v.to_csv(V/'visit_biomass_support.csv',index=False)
x.to_parquet(P/'supported_measurement_uses.parquet',index=False)
env=dict(os.environ,INTEGRATED_REVIEW_VARIANT=name)
for script in ['estimate.py','joint_models.py']:subprocess.run(['python3',str(R/'code/integrated_review'/script)],cwd=R,env=env,check=True)
print(name,'complete')
