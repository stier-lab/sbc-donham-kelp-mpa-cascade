#!/usr/bin/env python3
"""Independent scope, annual-value, period and primary-route checks for the PISCO rule."""
import csv,json,math,hashlib
from pathlib import Path
R=Path(__file__).resolve().parents[2]
def read(p):return list(csv.DictReader((R/p).open()))
policy=json.loads((R/'data/harmonized/pisco_baseline_policy.json').read_text())
assert policy['minimum_before_years']==3 and policy['sensitivity_minimum']==4
annual=read('tables/table_s_pisco_fixed_annual.csv')
fixed={(x['MPA'],x['Resp'],x['year']):x for x in annual if x['scenario']=='strict_min3' and x['membership']=='fixed'}
assert len(fixed)==84
for name,old,taxcol in [('harmonized_response_ratios.csv','ratios_before.csv','y'),('harmonized_raw_responses.csv','raw_before.csv','taxon_name')]:
 current=read('data/harmonized/'+name);before=read('data/external/pisco_policy_20261031/'+old)
 target=lambda x:x['source']=='PISCO' and x[taxcol]=='Semicossyphus pulcher' and x['CA_MPA_Name_Short'] in policy['MPAs']
 assert [x for x in current if not target(x)]==[x for x in before if not target(x)]
 assert hashlib.md5((R/'data/harmonized'/name).read_bytes()).hexdigest()==policy['files'][name]['after_md5']
 for x in filter(target,current):
  f=fixed[x['CA_MPA_Name_Short'],x['resp'],x['year']]
  assert x['BA']==('Before' if int(x['year'])<2012 else 'After')
  if taxcol=='y':
   assert abs(float(x['lnDiff'])-math.log((float(f['mpa'])+float(f['correction']))/(float(f['reference'])+float(f['correction']))))<1e-10
   assert abs(float(x['Diff'])-math.exp(float(x['lnDiff'])))<1e-10
  else:assert abs(float(x['value'])-float(f[x['status']]))<1e-8
primary=read('data/sumstats_final.csv');old=read('data/external/pisco_policy_20261031/primary_before.csv')
key=lambda x:tuple(x[k] for k in ['MPA','Taxa','Source','Resp'])
historical={key(x):x for x in old}
for x in primary:
 target=x['Source']=='PISCO' and x['Taxa']=='S. pulcher' and x['MPA'] in policy['MPAs']
 if target:
  assert x['AnalysisType']=='pBACIPS'
  assert x['Model']==('Linear' if x['MPA']=='Point Vicente SMCA' else 'Step')
  yy=[int(k[2]) for k in fixed if k[:2]==(x['MPA'],x['Resp'])]
  assert len(yy)==int(x['N']) and sum(y<2012 for y in yy)>=3 and sum(y>=2012 for y in yy)>=4
 else:
  for c in ['Mean','SE','N','Model','AnalysisType']:assert x[c]==historical[key(x)][c],(key(x),c)
assert sum(x['AnalysisType']=='pBACIPS' for x in primary)==42
assert sum(int(x['N']) for x in primary)==2815
print('PASS: only six histories changed; source values, strict prior years, 42 before/after contrasts and 2815 annual rows agree.')
