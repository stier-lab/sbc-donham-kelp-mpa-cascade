"""Run the approved coverage and baseline sensitivities without changing main data."""
from common import *
import subprocess,sys,shutil
for variant in ['constant_stations','constant_zones','four_before','sparse_size_variance','exact_KFM_dates','second_seed']:
 target=BASE_OUT/'sensitivities'/variant;target.mkdir(parents=True,exist_ok=True);(target/'private').mkdir(exist_ok=True)
 for name in (['comparison_year_roster.csv','station_assignments.csv'] if variant=='exact_KFM_dates' else ['comparison_year_roster.csv','visit_biomass_support.csv','size_pool_inventory.csv','station_assignments.csv']):
  dst=target/name
  if not dst.exists():dst.symlink_to(BASE_OUT/name)
 for name in (['measurement_uses.parquet','lobster_imputation_members.csv'] if variant=='exact_KFM_dates' else ['supported_measurement_uses.parquet','size_pool_members.parquet']):
  dst=target/'private'/name
  if not dst.exists():dst.symlink_to(BASE_OUT/'private'/name)
 env=os.environ.copy();env['INTEGRATED_REVIEW_VARIANT']=variant;env['OPENBLAS_NUM_THREADS']='1'
 for script in (['biomass','estimate','joint_models'] if variant=='exact_KFM_dates' else ['estimate','joint_models']):
  with (target/(script+'.log')).open('w') as f:subprocess.run([sys.executable,str(ROOT/'code/integrated_review'/f'{script}.py')],cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
 print(variant,flush=True)
parts=[]
for variant in ['', 'constant_stations','constant_zones','four_before','sparse_size_variance','exact_KFM_dates','second_seed']:
 folder=BASE_OUT if not variant else BASE_OUT/'sensitivities'/variant
 p=pd.read_csv(folder/'pooled_results.csv');p['coverage_variant']=variant or 'available';parts.append(p)
pd.concat(parts,ignore_index=True).to_csv(BASE_OUT/'pooled_coverage_sensitivities.csv',index=False)
