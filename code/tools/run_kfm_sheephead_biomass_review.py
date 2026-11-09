"""Rebuild the restricted KFM review, its sensitivities, report and guide tables."""
from pathlib import Path
import os,subprocess,sys
A=Path(__file__).resolve().parents[2]
env=dict(os.environ,OPENBLAS_NUM_THREADS='1')
base=[sys.executable,'code/tools/kfm_sheephead_biomass_review.py']
variants=[[],['--variant','empirical','--count-model','empirical'],['--variant','pool_strength4','--pool-strength','4'],['--variant','count_prior01','--count-prior','.1'],['--variant','dependence_positive','--dependence','.5'],['--variant','dependence_negative','--dependence','-.5']]
for flags in variants:subprocess.run(base+flags,cwd=A,env=env,check=True)
for script in ['report_kfm_sheephead_biomass.py','check_kfm_sheephead_biomass.py','report_dataset_inventory.py','check_current_decisions.py','check_documentation_consolidation.py']:
 subprocess.run([sys.executable,'code/tools/'+script],cwd=A,env=env,check=True)
