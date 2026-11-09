"""One review entry point. Never promotes results to the manuscript."""
from pathlib import Path
import argparse,subprocess,sys,os,json,hashlib,datetime
A=Path(__file__).resolve().parents[2];D=A.parent/'sbc-donham-kelp-mpa-data';O=A/'outputs/integrated_review_20261103';O.mkdir(exist_ok=True);os.chdir(A)
ap=argparse.ArgumentParser();ap.add_argument('--analysis-only',action='store_true',help='Skip isolated legacy reconstruction; still rebuild and test registered raw sources.');args=ap.parse_args()
(O/'completed.json').unlink(missing_ok=True);(O/'input_hashes.json').unlink(missing_ok=True);env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
def run(cmd,cwd=A,name='phase'):
 print(name,flush=True)
 with (O/(name+'.log')).open('w') as f:subprocess.run(cmd,cwd=cwd,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
code_before={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in (A/'code/integrated_review').glob('*.py')}
run(['make','check'],D,'source_checks')
run([sys.executable,'code/integrated_review/test_uncertainty.py'],name='uncertainty_tests')
if not args.analysis_only:
 run([sys.executable,'code/integrated_review/stage_sources.py'],name='stage_sources')
 run(['Rscript','integrated_source_run.R'],O/'stage','legacy_reconstruction')
for phase in ['raw_measurements','coverage','biomass','estimate','joint_models','sensitivities','check','report']:
 run([sys.executable,'code/integrated_review/'+phase+'.py'],name=phase)
assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in code_before.items()),'Pipeline code changed during run; repeat before certification'
(O/'completed.json').write_text(json.dumps(dict(completed=datetime.datetime.now(datetime.timezone.utc).isoformat(),paper_updated=False,source_code_hashes={str(p.relative_to(A)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (A/'code/integrated_review').glob('*.py')},policy_sha256=hashlib.sha256((A/'config/integrated_review_policy.json').read_bytes()).hexdigest()),indent=2))
