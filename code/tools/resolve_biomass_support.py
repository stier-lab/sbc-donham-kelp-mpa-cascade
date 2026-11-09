"""Review-only size-support sensitivities. Never writes manuscript inputs.
Shared source-row Bayesian weights preserve dependence across reused animals.
Draw ranges describe size-driven variation in slope estimates, not confidence intervals.
"""
from pathlib import Path
import ast,json,hashlib
import numpy as np,pandas as pd
from scipy import sparse
from scipy.stats import t as student
O=Path('outputs/biomass_resolution_20261102');O.mkdir(exist_ok=True)
B=Path('outputs/biomass_support_20261102');A=Path('outputs/pisco_fish_rules_20261101')
# Reuse the previously independently verified calendar-gap HAC implementation only.
tree=ast.parse(Path('code/tools/audit_station_coverage_refits.py').read_text())
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['fit','slope']],type_ignores=[]),'verified_fit','exec'))
a=pd.read_csv(A/'approved_annual_input.csv');p=pd.read_csv(B/'private/pool_members.csv',low_memory=False)
i=pd.read_csv(B/'animal_pool_inventory.csv');off=pd.read_csv('outputs/all_program_station_refits_20261102/offsets.csv').set_index('series_id').offset
p['mass']=np.select([p.Taxa.eq('M. franciscanus'),p.Taxa.eq('S. purpuratus')],[.00059*p['size']**2.917,.00059*p['size']**2.870],default=.001352821*p['size']**2.913963113)
keys=['series_id','year','status'];pools=p[keys].drop_duplicates().reset_index(drop=True);pools['pool_id']=np.arange(len(pools));p=p.merge(pools,on=keys,validate='many_to_one')
# Each raw record's frequency denotes animals; keep original row identity across uses.
raw=p[['source_row','count','mass','survey_group']].drop_duplicates();assert not raw.source_row.duplicated().any()
raw=raw.sort_values('source_row').reset_index(drop=True);rowidx=pd.Series(np.arange(len(raw)),index=raw.source_row);p['j']=p.source_row.map(rowidx)
W=sparse.csr_matrix((p['count'],(p.pool_id,p.j)),shape=(len(pools),len(raw)))
N=np.asarray(W.sum(axis=1)).ravel();mass=raw.mass.to_numpy();means=np.asarray(W@mass).ravel()/N
pools['mean_mass']=means;pools['n_animals']=N
lookup={(r.series_id,r.year,r.status):r.pool_id for r in pools.itertuples()}
inv=i.set_index(keys);effects=[];draws=[];losses=[];annual=[]
# Conditional Bayesian bootstrap; singleton variability is deliberately not invented.
rng=np.random.default_rng(20261102);weights=rng.gamma(raw['count'].to_numpy()[:,None],1,size=(len(raw),500))/raw['count'].to_numpy()[:,None]
boot=(W@(weights*mass[:,None]))/(W@weights);del weights
ratio=boot/means[:,None]
# A separate group-weight stress check retains shared group IDs across histories.
# Missing identifiers can merge visits, so this is not certified dive-level resampling.
groups,labels=pd.factorize(raw.survey_group);gw=rng.exponential(size=(len(labels),500))[groups]
group_ratio=((W@(gw*mass[:,None]))/(W@gw))/means[:,None];del gw
# Borrowing alternatives retain only supported paired years in BOTH comparison fits.
restricted={}
for flag in ['same_program','same_count_site']:
 q=p[p[flag].astype(str).str.lower().eq('true')];s=q.groupby('pool_id').apply(lambda z:np.average(z.mass,weights=z['count']),include_groups=False)
 restricted[flag]=s.reindex(pools.pool_id).to_numpy()/means
for sid,z in a[a.t>0].groupby('series_id'):
 if sid not in set(pools.series_id):continue
 z=z.sort_values('year').copy();base=slope(z,'lnDiff');assert base
 ids={st:np.array([lookup.get((sid,y,st),-1) for y in z.year]) for st in ['mpa','reference']}
 vals={st:np.array([inv.loc[(sid,y,st),'value'] for y in z.year]) for st in ids}
 def get_ratio(r,st):
  v=np.ones((len(z),)+r.shape[1:]);ok=ids[st]>=0;v[ok]=r[ids[st][ok]];assert np.all(vals[st][~ok]==0);return v
 def response(r):
  x=get_ratio(r,'mpa');y=get_ratio(r,'reference');ext=(slice(None),)+(None,)*(r.ndim-1)
  return np.log((vals['mpa'][ext]*x+off[sid])/(vals['reference'][ext]*y+off[sid]))
 assert np.max(np.abs(response(np.ones(len(pools)))-z.lnDiff))<1e-8
 X=np.column_stack([np.ones(len(z)),z.t]);bb=(np.linalg.pinv(X)@response(ratio))[1]
 draws.append(dict(series_id=sid,scheme='individual',estimate=base['estimate'],size_draw_low=np.quantile(bb,.025),size_draw_high=np.quantile(bb,.975),draws=500,n_years=len(z)))
 gb=(np.linalg.pinv(X)@response(group_ratio))[1]
 draws.append(dict(series_id=sid,scheme='recorded_group',estimate=base['estimate'],size_draw_low=np.quantile(gb,.025),size_draw_high=np.quantile(gb,.975),draws=500,n_years=len(z)))
 for label,r in [('current',np.ones(len(pools)))]+list(restricted.items()):
  zz=z.copy();zz['new']=response(r);zz=zz[np.isfinite(zz.new)];before=slope(zz,'lnDiff');after=slope(zz,'new')
  losses.append(dict(series_id=sid,scenario=label,n_original=len(z),n_supported=len(zz),lost_years=';'.join(map(str,sorted(set(z.year)-set(zz.year))))))
  if before and after:effects.append(dict(series_id=sid,scenario=label,matched_original=before['estimate'],**after))
 # Explicit stress test: independently halve/double singleton pool mean mass by role.
 for direction in [-1,1]:
  r=np.ones(len(pools));mask=pools.n_animals.eq(1);r[mask & pools.status.eq('mpa')]=2.**direction;r[mask & pools.status.eq('reference')]=2.**(-direction)
  zz=z.copy();zz['new']=response(r);f=slope(zz,'new');effects.append(dict(series_id=sid,scenario='singleton_stress_'+str(direction),matched_original=base['estimate'],**f))
# Correct unknown kelp biomass in review input by withholding paired year, not a station.
bad=pd.read_csv(B/'unresolved_direct_measurements.csv')
# Data-driven exclusion joins; the same rule applies to every flagged biomass station-year.
badkeys=bad[['series_id','year']].drop_duplicates().assign(unknown_measurement=True)
marked=a.merge(badkeys,on=['series_id','year'],how='left',validate='one_to_one')
corrected=marked[marked.unknown_measurement.isna()].drop(columns='unknown_measurement')
assert len(a)-len(corrected)==len(badkeys)
corrected.to_csv(O/'review_annual_input.csv',index=False)
badkeys['reason']='Observed organisms lack size measurements; biomass unknown. Omit paired biomass year to preserve reference membership; retain raw counts.'
badkeys.to_csv(O/'measurement_exclusions.csv',index=False)
for name,obj in [('conditional_size_draws',draws),('matched_donor_effects',effects),('donor_coverage_losses',losses)]:pd.DataFrame(obj).to_csv(O/(name+'.csv'),index=False)
# Red-urchin audit preserves program-specific sampled populations.
red=i[i.Taxa.eq('M. franciscanus')].copy();red['size_rule']=np.where(red.Source.eq('PISCO'),'Test diameter >=25 mm','Unfiltered measured diameter; count size classes need source confirmation')
red['handling']='Retain density; biomass provisional where donor representativeness is unresolved'
red.to_csv(O/'red_urchin_role_years.csv',index=False)
r=p[p.Taxa.eq('M. franciscanus')];rows=[]
for sid,g in r.groupby('series_id'):
 rows.append(dict(series_id=sid,n_source_rows=g.source_row.nunique(),n_role_years=g.pool_id.nunique(),min_diameter_mm=g['size'].min(),max_diameter_mm=g['size'].max(),fraction_animals_below_25=np.average(g['size']<25,weights=g['count']),fraction_same_count_site=np.average(g.same_count_site,weights=g['count']),fraction_complete_dates=np.average(g.full_date,weights=g['count'])))
pd.DataFrame(rows).to_csv(O/'red_urchin_summary.csv',index=False)
# Exact deterministic verification and independent resampling seed diagnostic.
old=pd.read_csv(A/'effects.csv');old=old[(old.scenario=='Total__After')&(old.metric=='Post_slope')]
e=pd.DataFrame(effects);check=e[e.scenario.eq('current')].merge(old,on='series_id',suffixes=('_new','_old'))
assert len(check)*2==len(draws);assert max(abs(check.estimate_new-check.estimate_old))<1e-8;assert max(abs(check.se_new-check.se_old))<1e-8
protected=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text())
for path,h in protected.items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==h
(O/'verification.json').write_text(json.dumps(dict(histories_checked=len(check),shared_source_rows=len(raw),pool_occurrences=len(p),conditional_draws=500,corrected_annual_rows=len(corrected),protected_files_unchanged=len(protected),red_histories=len(rows)),indent=2))
print((O/'verification.json').read_text())
