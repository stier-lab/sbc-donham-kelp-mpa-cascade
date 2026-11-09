"""Isolated correction of proven excluded-donor leakage, using original <=2011 sizes."""
from pathlib import Path
import pandas as pd,numpy as np,ast
from scipy.stats import t as student
O=Path('outputs/biomass_exclusion_reconciliation_20261102');D=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data')
tree=ast.parse(Path('code/tools/audit_station_coverage_refits.py').read_text());exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['fit','slope']],type_ignores=[]),'fit','exec'))
a=pd.read_csv('outputs/biomass_resolution_20261102/review_annual_input.csv');inv=pd.read_csv('outputs/biomass_support_20261102/animal_pool_inventory.csv').set_index(['series_id','year','status']);off=pd.read_csv('outputs/all_program_station_refits_20261102/offsets.csv').set_index('series_id').offset
p=pd.read_csv('outputs/biomass_support_20261102/private/pool_members.csv',low_memory=False);l=pd.read_csv(D/'source_archive/kfm_size_feasibility_2026-10-29/legacy/Invertebrate Size Frequency.csv');affected=pd.read_csv(O/'affected_histories.csv');changes=[];effects=[];variant=a.copy()
for item in affected.itertuples():
 sid=item.series_id;site=item.delivered_site;code=11005 if 'M. franciscanus' in sid else 11006;power=2.917 if code==11005 else 2.870;st='reference' if site=='a-k-01' else 'mpa';z=a[a.series_id.eq(sid)&a.t.gt(0)].copy();z['new']=z.lnDiff;reconciled=[]
 for year in z.loc[z.year.le(2011),'year']:
  pool=p[p.series_id.eq(sid)&p.year.eq(year)&p.status.eq(st)];target=pool[pool.campus.eq('KFM')&pool.site.eq(site)]
  if target.empty:continue
  exact=l[l.Year.eq(year)&l.Species.eq(code)&l.Site.isin([1,21] if site=='a-k-01' else [37])].copy();replacement=l[l.Year.eq(year)&l.Species.eq(code)&l.Site.eq(1 if site=='a-k-01' else 14)].copy()
  for q in [exact,replacement]:q['mm']=(q.Size*10).round(6)
  if '|PISCO|' in sid:exact=exact[exact.mm.ge(25)];replacement=replacement[replacement.mm.ge(25)]
  if exact.empty or replacement.empty:continue
  hist=target.groupby('size')['count'].sum().subtract(exact.groupby('mm').size(),fill_value=0);assert hist.abs().max()<1e-9
  keep=pool[~(pool.campus.eq('KFM')&pool.site.eq(site))]
  oldmean=np.average(.00059*pool['size']**power,weights=pool['count']);newmean=((.00059*keep['size']**power*keep['count']).sum()+(.00059*replacement.mm**power).sum())/(keep['count'].sum()+len(replacement));factor=newmean/oldmean
  vv={stt:inv.loc[(sid,year,stt),'value'] for stt in ['mpa','reference']};vv[st]*=factor;new=np.log((vv['mpa']+off[sid])/(vv['reference']+off[sid]));z.loc[z.year.eq(year),'new']=new;variant.loc[variant.series_id.eq(sid)&variant.year.eq(year),'lnDiff']=new;reconciled.append(year)
  changes.append(dict(series_id=sid,year=year,status=st,delivered_site=site,removed_original_records=int(target['count'].sum()),allowed_replacement_records=len(replacement),old_mean_mass=oldmean,new_mean_mass=newmean,mass_factor=factor,old_lnDiff=float(z.loc[z.year.eq(year),'lnDiff'].iloc[0]),new_lnDiff=new))
 for scope,zz in [('all_observed_post_years',z),('reconciled_years_only',z[z.year.isin(reconciled)])]:
  for label,col in [('current','lnDiff'),('allowed_original_sizes','new')]:
   f=slope(zz,col)
   if f:effects.append(dict(series_id=sid,scope=scope,scenario=label,n_reconciled_years=len(reconciled),**f))
pd.DataFrame(changes).to_csv(O/'annual_donor_corrections.csv',index=False);pd.DataFrame(effects).to_csv(O/'effects.csv',index=False);variant.to_csv(O/'partial_reconciled_review_variant.csv',index=False)
print(pd.DataFrame(effects).query("scope=='all_observed_post_years'")[['series_id','scenario','estimate','low','high','n_reconciled_years']].to_string(index=False))
