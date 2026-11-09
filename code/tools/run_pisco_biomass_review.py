"""Single entry point for the isolated, non-production PISCO biomass review."""
from pathlib import Path
import subprocess,sys
root=Path(__file__).resolve().parents[2]
for cmd in [[sys.executable,'code/tools/propagate_pisco_borrowing.py'],['Rscript','code/tools/verify_pisco_borrowing.R'],[sys.executable,'code/tests/check_pisco_borrowing.py'],[sys.executable,'code/tools/report_pisco_borrowing.py']]:subprocess.run(cmd,cwd=root,check=True)
