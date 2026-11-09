"""Independent source, aggregation and result checks; no primary changes."""
import csv, hashlib, math
from pathlib import Path
from collections import defaultdict
root=Path(__file__).resolve().parents[2]
def read(p):
 with p.open() as f:return list(csv.DictReader(f))
def tab(n):return read(root/f'tables/table_s_lobster_missing_{n}.csv')
for x in tab('inputs'):
 p=Path(x['path']);p=p if p.is_absolute() else root/p
 assert hashlib.md5(p.read_bytes()).hexdigest()==x['md5'],p
src=read(root/'data/LTER/Lobster_Abundance_All_Years_20230831.csv')
raw=[x for x in src if x['SITE'] in ('IVEE','NAPL')]
missing=[x for x in raw if float(x['COUNT'])>0 and float(x['SIZE_MM'])<=0]
assert len(missing)==5 and all(x['NUM_AO']==x['COUNT']=='1' for x in missing)
for site,year in [('IVEE','2017'),('NAPL','2014')]:
 donors=[x for x in src if x['SITE']==site and x['YEAR']==year and float(x['COUNT'])>0 and float(x['SIZE_MM'])>0]
 targets=[x for x in missing if x['SITE']==site and x['YEAR']==year]
 assert {x['DATE'] for x in donors}=={x['DATE'] for x in targets}
 # Independent nested survey arithmetic for the change in the affected protected site.
 def biomass(use_mean):
  mean=sum(float(x['COUNT'])*.001352821*float(x['SIZE_MM'])**2.913963113 for x in donors)/sum(float(x['COUNT']) for x in donors)
  groups=defaultdict(float)
  for x in src:
   if x['SITE']!=site or x['YEAR']!=year:continue
   n=float(x['COUNT']);s=float(x['SIZE_MM'])
   b=n*.001352821*s**2.913963113 if n>0 and s>0 else (n*mean if n>0 and use_mean else 0)
   groups[x['TRANSECT'],x['REPLICATE']]+=b
  ts=defaultdict(list)
  for (t,r),b in groups.items():ts[t].append(b)
  return sum(sum(v)/len(v) for v in ts.values())/len(ts)/300
 expected=100*(biomass(True)/biomass(False)-1)
 observed=next(x for x in tab('annual') if x['scenario']=='site_year_mean_mass' and x['year']==year and x['status']=='mpa' and (x['MPA'].startswith('Campus') if site=='IVEE' else x['MPA'].startswith('Naples')))
 assert math.isclose(expected,float(observed['pct_change']),abs_tol=1e-10)
pooled=tab('pooled');base={(x['Taxa'],x['Resp']):x for x in pooled if x['scenario']=='current_zero'}
for x in pooled:
 b=base[x['Taxa'],x['Resp']]
 assert float(x['estimate'])*float(b['estimate'])>0
 assert (float(x['q'])<.05)==(float(b['q'])<.05)
primary_perm=next(x for x in read(root/'tables/table_s_prereview_crosstaxa_permutation.csv') if x['Relationship']=='P. interruptus biomass -> S. purpuratus biomass')
for x in tab('pathway'):
 assert x['k']=='4' and float(x['exact_p'])>.05
 assert abs(float(x['exact_p'])-float(primary_perm['p_permutation']))<.00005
tr=tab('temporal');b={x['Species']:x for x in tr if x['scenario']=='current_zero'}
for x in tr:
 assert float(x['estimate'])*float(b[x['Species']]['estimate'])>0
 assert (float(x['p'])<.05)==(float(b[x['Species']]['p'])<.05)
assert all(x['singular']=='FALSE' and x['convergence']=='' for x in tab('fit_checks'))
print('PASS: five antenna-only records; same-date donors; independent hierarchical changes; all six scenarios preserve pooled/temporal decisions and the same inconclusive exact pathway test. Source fingerprints current.')
