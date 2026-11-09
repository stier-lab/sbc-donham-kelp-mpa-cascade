"""Rebuild corrected whole-dataset results and the separate KFM fish review.

The driver writes review outputs only. It never runs the manuscript's legacy
`run_all.R` or edits frozen tables, figures, manuscript, or supplement.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser()
parser.add_argument("--reuse-integrated", action="store_true", help="Use a completed, hash-verified integrated build")
parser.add_argument("--reuse-kfm", action="store_true", help="Use the separately verified KFM sheephead review")
args = parser.parse_args()
env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")


def run(path: str) -> None:
    print(f"Running {path}", flush=True)
    subprocess.run([sys.executable, path], cwd=ROOT, env=env, check=True)


if not args.reuse_integrated:
    run("code/integrated_review/run.py")
if not args.reuse_kfm:
    run("code/tools/run_kfm_sheephead_biomass_review.py")
else:
    folder = ROOT / "outputs/kfm_sheephead_biomass_20261103"
    checks = json.loads((folder / "verification.json").read_text())
    assert checks["input_hashes_verified"] and checks["eligible_comparison_years"] == 52
    for record in json.loads((folder / "inputs.json").read_text()):
        path = Path(record["path"])
        assert path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest() == record["sha256"], path
run("code/tools/report_whole_dataset_biomass_update.py")
run("code/tests/check_whole_dataset_biomass_review.py")
