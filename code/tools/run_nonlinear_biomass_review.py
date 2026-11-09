"""Rebuild conditional nonlinear calibration and documentation; no paper changes.
The Scott receipt is a dated verified input, not an email query repeated by this runner.
"""
from pathlib import Path
import os,sys,subprocess
os.chdir(Path(__file__).resolve().parents[2]);env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1');py=sys.executable
for script in ['nonlinear_biomass_calibration.py','test_nonlinear_biomass_calibration.py','report_nonlinear_biomass_calibration.py','report_dataset_inventory.py','check_current_decisions.py','verify_consolidated_inventory.py']:
    subprocess.run([py,'code/tools/'+script],check=True,env=env)
