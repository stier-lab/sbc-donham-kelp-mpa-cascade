#!/usr/bin/env python3
"""Guard matching, primary reconstruction and reported sensitivity conclusions."""
import csv,hashlib,json,math
from pathlib import Path
r=Path(__file__).resolve().parents[2]
def rows(name):return list(csv.DictReader((r/'tables'/name).open()))
m=json.loads((r/'tables/kfm_continuity_manifest.json').read_text())
for path,digest in m['inputs_sha256'].items():assert hashlib.sha256((r/path).read_bytes()).hexdigest()==digest,path
assert len(m['after_only_in_pbacips'])==6
assert all(x['Taxa']=='S. pulcher' for x in m['after_only_in_pbacips'])
a=rows('table_s_kfm_continuity_annual.csv')
assert len(a)==74 and all(int(x['year'])<=2022 for x in a)
assert all(int(x['mpa_n_persistent'])==int(x['reference_n_persistent'])==1 for x in a)
assert all(float(x['offset'])>0 for x in a)
c=rows('table_s_kfm_continuity_contrasts.csv')
for mpa in {x['MPA'] for x in c}:
    x,y=[x for x in c if x['MPA']==mpa]
    assert x['years']==y['years'] and x['n_years']==y['n_years']=='37'
    assert float(x['estimate'])*float(y['estimate'])>0
s=next(x for x in c if x['MPA']=='Scorpion SMR' and x['scenario']=='persistent_sites')
assert float(s['HAC_ci_low'])<0<float(s['HAC_ci_high'])
p=rows('table_s_kfm_continuity_pooled.csv')
base={(x['Taxa'],x['Resp']):x for x in p if x['scenario']=='primary'}
for x in p:
    b=base[x['Taxa'],x['Resp']]
    assert float(x['estimate'])*float(b['estimate'])>0
    assert (float(x['q'])<.05)==(float(b['q'])<.05)
for x in rows('table_s_kfm_continuity_inputs.csv'):
    assert hashlib.md5((r/x['path']).read_bytes()).hexdigest()==x['md5'],x['path']
print('PASS: matched calendars, observed-zero retention, current inputs, original reconstruction and nine pooled conclusions')
