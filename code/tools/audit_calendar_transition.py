"""Compare the surveyed-calendar primary policy with its frozen preceding fit."""
import csv, json, math, hashlib
from pathlib import Path
R=Path(__file__).resolve().parents[2]; B=R/'data/external/calendar_policy_20261030'
read=lambda p:list(csv.DictReader(p.open()))
key=lambda x:tuple(x[c] for c in ('MPA','Taxa','Source','Resp'))
a=read(R/'data/sumstats_final.csv');o={key(x):x for x in read(B/'primary_before.csv')}
assert len(a)==len(o)==145 and {key(x) for x in a}==set(o)
rows=[]
for r in a:
 old=o[key(r)]
 assert all(r[c]==old[c] for c in ('Model','N','AnalysisType')),(key(r),'Changed membership or retained model')
 rows.append({**{c:r[c] for c in ('MPA','Taxa','Source','Resp','Model','N')},
              'previous_estimate':old['Mean'],'estimate':r['Mean'],
              'previous_se':old['SE'],'se':r['SE'],
              'changed':abs(float(old['Mean'])-float(r['Mean']))>1e-7})
def write(name,rows):
 with (R/'tables'/name).open('w',newline='')as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
write('table_s_calendar_transition_histories.csv',rows)
old={ (r['Taxa'],r['Response']):r for r in read(B/'pooled_before.csv')}
pooled=[]
for r in read(R/'tables/table_02_meta_analysis.csv'):
 p=old[(r['Taxa'],r['Response'])]
 pooled.append({**{c:r[c]for c in ['Taxa','Response']},**{'previous_'+c:p[c]for c in ['Estimate','CI_lower','CI_upper','pval_fdr']},
                **{c:r[c]for c in ['Estimate','CI_lower','CI_upper','pval_fdr']},
                'direction_changed':math.copysign(1,float(r['Estimate']))!=math.copysign(1,float(p['Estimate'])),
                'adjusted_decision_changed':(float(p['pval_fdr'])<.05)!=(float(r['pval_fdr'])<.05)})
write('table_s_calendar_transition_pooled.csv',pooled)
print(json.dumps({'same_histories':len(a),'changed_estimates':sum(r['changed']for r in rows),
 'changed_adjusted_conclusions':[{c:r[c]for c in ['Taxa','Response']}for r in pooled if r['adjusted_decision_changed']],
 'changed_directions':sum(r['direction_changed']for r in pooled)},indent=2))
