"""Review-only direct-size sensitivity, preserving existing density and survey roster."""
from pathlib import Path
import os,subprocess,shutil
import pandas as pd
R=Path(__file__).resolve().parents[2];M=R/'outputs/integrated_review_20261103';V=M/'sensitivities/direct_sizes_only';P=V/'private';P.mkdir(parents=True,exist_ok=True)
for name in ['comparison_year_roster.csv','station_assignments.csv','visit_biomass_support.csv','size_pool_inventory.csv','input_hashes.json']:
 shutil.copy2(M/name,V/name)
for name in ['measurement_uses.parquet','size_pool_members.parquet','lobster_imputation_members.csv']:
 shutil.copy2(M/'private'/name,P/name)
x=pd.read_parquet(M/'private/supported_measurement_uses.parquet');v=pd.read_csv(M/'visit_biomass_support.csv')
b=v[v.basis.str.startswith('calibrated_',na=False)][['series_id','site','year','date','zone']].drop_duplicates();assert len(b)
z=x.merge(b.assign(borrowed=True),on=['series_id','site','year','date','zone'],how='left',validate='many_to_one');m=z.borrowed.fillna(False)&z.Resp.eq('Bio');z.loc[m,'biomass']=float('nan');z.loc[m,'size_basis']='borrowed_sizes_withheld';z.loc[m,'donor_pool']='';z=z.drop(columns='borrowed');z.to_parquet(P/'supported_measurement_uses.parquet',index=False)
v.loc[v.basis.str.startswith('calibrated_',na=False),'biomass_supported']=False
v.to_csv(V/'visit_biomass_support.csv',index=False)
env=dict(os.environ,INTEGRATED_REVIEW_VARIANT='direct_sizes_only')
for script in ['estimate.py','joint_models.py']:
 subprocess.run(['python3',str(R/'code/integrated_review'/script)],cwd=R,env=env,check=True)
base=pd.read_csv(M/'annual_values.csv');alt=pd.read_csv(V/'annual_values.csv');assert len(base)==len(alt);assert (base[['series_id','year']]==alt[['series_id','year']]).all().all()
z=base[['series_id','year','Source','Taxa','Resp','response_supported']].merge(alt[['series_id','year','response_supported']],on=['series_id','year'],suffixes=('_current','_direct'),validate='one_to_one');assert not (z.response_supported_direct&~z.response_supported_current).any();z.to_csv(V/'coverage_changes.csv',index=False)
print('borrowed visits',len(b),'lost annual supported years',((z.response_supported_current)&(~z.response_supported_direct)).sum())
