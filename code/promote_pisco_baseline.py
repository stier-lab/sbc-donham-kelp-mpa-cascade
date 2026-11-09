#!/usr/bin/env python3
"""Rebuild the authorized six fixed-station PISCO inputs after module 59.

Default verifies the derived values without changing files. --apply writes the two
harmonized tables and refreshes only this policy marker's output fingerprints.
"""
import csv,json,math,hashlib,io,argparse
from pathlib import Path
R=Path(__file__).resolve().parents[1]
read=lambda p:list(csv.DictReader(p.open()))
ref=R/'data/external/pisco_policy_20261031'
p=R/'data/harmonized/pisco_baseline_policy.json';policy=json.loads(p.read_text())
fixed=[x for x in read(R/'tables/table_s_pisco_fixed_annual.csv') if x['scenario']=='strict_min3' and x['membership']=='fixed']
lookup={(x['MPA'],x['Resp'],x['year']):x for x in fixed}
assert len(lookup)==84 and policy['minimum_before_years']==3
raw=read(ref/'raw_before.csv')
maximum={}
for x in raw:
 if x['source']=='PISCO' and x['taxon_name']=='Semicossyphus pulcher':
  k=x['CA_MPA_Name_Short'],x['resp'];maximum[k]=max(maximum.get(k,0),float(x['value']))
apply=argparse.ArgumentParser();apply.add_argument('--apply',action='store_true');args=apply.parse_args()
for name,prior,taxon in [('harmonized_response_ratios.csv','ratios_before.csv','y'),('harmonized_raw_responses.csv','raw_before.csv','taxon_name')]:
 source=ref/prior;rows=read(source);out=[]
 assert hashlib.md5(source.read_bytes()).hexdigest()==policy['files'][name]['before_md5']
 for row in rows:
  x=dict(row)
  if x['source']=='PISCO' and x[taxon]=='Semicossyphus pulcher' and x['CA_MPA_Name_Short'] in policy['MPAs']:
   k=x['CA_MPA_Name_Short'],x['resp'],x['year']
   if k not in lookup:continue
   v=lookup[k];x['BA']='Before' if int(x['year'])<2012 else 'After'
   if taxon=='taxon_name':x['value']=v[x['status']]
   else:
    a=float(v['mpa'])+float(v['correction']);b=float(v['reference'])+float(v['correction']);scale=maximum[k[:2]]
    x.update(mpa=str(a/scale),reference=str(b/scale),lnDiff=v['lnDiff'],Diff=str(math.exp(float(v['lnDiff']))))
  out.append(x)
 path=R/'data/harmonized'/name
 if args.apply:
  with path.open('w',newline='')as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(out)
  policy['files'][name]['after_md5']=hashlib.md5(path.read_bytes()).hexdigest()
 else:
  current=read(path);assert len(out)==len(current)
  for x,y in zip(out,current):
   for k in x:
    if x[k]==y[k]:continue
    try:assert math.isclose(float(x[k]),float(y[k]),rel_tol=1e-12,abs_tol=1e-12),(name,k,x[k],y[k])
    except ValueError:raise AssertionError((name,k,x[k],y[k]))
 print(('Wrote' if args.apply else 'Verified'),name,len(out),'rows')
if args.apply:p.write_text(json.dumps(policy,indent=2)+'\n')
