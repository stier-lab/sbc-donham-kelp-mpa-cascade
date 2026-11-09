"""Rebuild scoped review evidence; canonical analysis and manuscript are untouched."""
from pathlib import Path
import subprocess,sys,os
ROOT=Path(__file__).resolve().parents[2]
os.chdir(ROOT)
O=Path('outputs/biomass_next_steps_20261102')
if (O/'completed.json').exists():(O/'completed.json').unlink()
env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
def run(*args):subprocess.run(list(args),check=True,env=env)
py=sys.executable
run(py,'code/tools/rebuild_density_baseline.py','--output',str(O/'density'))
run(py,'code/tools/test_density_baseline.py','--output',str(O/'density'))
run(py,'code/tools/rebuild_borrowing_support.py','--repo','.', '--data-root','/Users/adrianstier/sbc-donham-kelp-mpa-data','--output',str(O/'borrowing'))
run(py,'code/tools/test_borrowing_support.py')
base=[py,'code/tools/biomass_observation_model.py','--visits',str(O/'borrowing/visit_support.csv'),'--count-units',str(O/'density/biomass_roster_count_units.csv'),'--tolerance','7','--draws','800','--parameter-bootstraps','40','--trend-bootstraps','500']
for name,extra in [('uncertainty',[]),('uncertainty_pool10',['--pool-df','10']),('uncertainty_seed2',['--seed','20261103'])]:run(*base,'--output',str(O/name),*extra)
run(py,'code/tools/test_biomass_observation_model.py')
run('Rscript','code/tools/verify_direct_biomass_process.R')
run(py,'code/tools/report_biomass_next_steps.py')
run(py,'code/tools/report_dataset_inventory.py')
run(py,'code/tools/check_current_decisions.py')
run(py,'code/tools/verify_consolidated_inventory.py')
