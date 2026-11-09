"""Collect isolated agent evidence without adopting recommendations as policy."""
from pathlib import Path
import csv,hashlib,json,shutil
O=Path('outputs/parallel_resolution_20261102')
S=Path('/tmp/donham_swarm_20261102')
rows=[];manifest=[]
prior={}
if (O/'agent_findings.csv').exists():
    with (O/'agent_findings.csv').open() as f:
        prior={(r['audit_area'],r['finding_id']):r['coordinator_disposition'] for r in csv.DictReader(f)}
for key in ['kfm','support','uncertainty','coverage']:
    src=S/key
    if not (src/'report.md').exists() or not (src/'findings.csv').exists():
        raise SystemExit(f'Missing deliverable: {key}')
    dst=O/key;dst.mkdir(exist_ok=True)
    # Copy only small auditable reports/scripts/tables, never bulk private source files.
    for p in src.iterdir():
        if p.is_file() and p.suffix in {'.md','.csv','.py','.R','.json','.txt','.log'} and p.stat().st_size<2_000_000:
            shutil.copy2(p,dst/p.name)
            manifest.append(dict(path=str(dst/p.name),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),agent=key))
    with (src/'findings.csv').open() as f:
        for r in csv.DictReader(f):
            r['audit_area']=key;r['coordinator_disposition']=prior.get((key,r['finding_id']),'Pending independent verification')
            rows.append(r)
cols=list(dict.fromkeys(k for r in rows for k in r))
with (O/'agent_findings.csv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=cols);w.writeheader();w.writerows(rows)
(O/'agent_evidence_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(f'Collected {len(rows)} findings; recommendations remain unadopted.')
