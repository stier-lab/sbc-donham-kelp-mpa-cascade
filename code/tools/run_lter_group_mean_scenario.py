"""Illustrative group-size CV sensitivity; assumed CV is not observed LTER data."""
from pathlib import Path
import os,shutil,subprocess,sys
import pandas as pd
R=Path(__file__).resolve().parents[2];M=R/'outputs/integrated_review_20261103';cv=float(sys.argv[1]);assert cv in [.2,.3]
name='lter_group_cv'+str(int(100*cv));V=M/'sensitivities'/name;P=V/'private';P.mkdir(parents=True,exist_ok=True)
for n in ['comparison_year_roster.csv','station_assignments.csv','visit_biomass_support.csv','size_pool_inventory.csv','input_hashes.json']:shutil.copy2(M/n,V/n)
for n in ['measurement_uses.parquet','size_pool_members.parquet','lobster_imputation_members.csv']:shutil.copy2(M/'private'/n,P/n)
x=pd.read_parquet(M/'private/supported_measurement_uses.parquet');z=x[(x.Source=='LTER')&(x.Resp=='Bio')&(x.Taxa.isin(['M. franciscanus','S. purpuratus']))];assert z.size_basis.eq('same_quadrat_size_class_mean').all()
for t,b in [('M. franciscanus',2.917),('S. purpuratus',2.870)]:
 f=(1+cv*cv)**(b*(b-1)/2);k=(x.Source=='LTER')&(x.Resp=='Bio')&x.Taxa.eq(t)&x.biomass.gt(0);x.loc[k,'biomass']*=f
x.to_parquet(P/'supported_measurement_uses.parquet',index=False)
env=dict(os.environ,INTEGRATED_REVIEW_VARIANT=name)
for script in ['estimate.py','joint_models.py']:subprocess.run(['python3',str(R/'code/integrated_review'/script)],cwd=R,env=env,check=True)
print(name,'complete')
