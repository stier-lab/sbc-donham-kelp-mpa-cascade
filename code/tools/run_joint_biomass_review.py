"""Rebuild the isolated joint-model review, checks and consolidated guide."""
from pathlib import Path
import subprocess,sys,os
os.chdir(Path(__file__).resolve().parents[2]);env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1');py=sys.executable
O=Path('outputs/joint_borrowed_biomass_20261102')
for f in ['report_verified.json','independent_checks.json']:
    if (O/f).exists():(O/f).unlink()
def run(*a):subprocess.run(a,check=True,env=env)
run(py,'code/tools/joint_borrowed_biomass.py')
for name,extra in [('seed2',['--seed','20261105']),('pool10',['--pool-df','10','--bootstraps','499']),('strict0',['--tolerance','0','--bootstraps','499'])]:run(py,'code/tools/joint_borrowed_biomass.py','--output',str(O/name),'--rhos','0',*extra)
run(py,'code/tools/check_joint_biomass_diagnostics.py')
run(py,'code/tools/test_joint_borrowed_biomass.py')
run(py,'code/tools/finalize_joint_biomass_review.py')
run(py,'code/tools/report_dataset_inventory.py')
run(py,'code/tools/check_current_decisions.py')
run(py,'code/tools/verify_consolidated_inventory.py')
