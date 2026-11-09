"""Review pipeline: all-month density and empirical predictive biomass draws.
Outer-site holdouts for donor-error calibration; no cross-protection borrowing.
Intervals describe conditional measurement/prediction variation, NOT total CIs.
"""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,ast
from scipy.stats import t as student
O=Path('outputs/pisco_biomass_borrowing_20261102');B=Path('outputs/pisco_all_months_20261102');V=Path('outputs/pisco_timing_donors_20261102');N=1000;rng=np.random.default_rng(20261102)
s=pd.read_csv(V/'size_sample_inventory.csv');s['date']=pd.to_datetime(s.date)
v=pd.read_csv(B/'visit_values.csv');v=v[v.scenario.eq('allmonths_sharedzones_visit_sizes')].copy();v['date']=pd.to_datetime(v.date)
a=pd.read_csv(B/'annual_values.csv');a=a[a.scenario.eq('allmonths_sharedzones_visit_sizes')];off=pd.read_csv('outputs/all_program_station_refits_20261102/offsets.csv').set_index('series_id').offset
# The count-led review table preserves all 63 post-protection comparison-years.
density=a[a.series_id.str.endswith('|Den')].copy();density.to_csv(O/'selected_density_annual.csv',index=False);assert (density.t>0).sum()==63
T=ast.parse(Path('code/tools/audit_station_coverage_refits.py').read_text());exec(compile(ast.Module(body=[n for n in T.body if isinstance(n,ast.FunctionDef) and n.name in ['fit','slope']],type_ignores=[]),'fit','exec'))
def predict(t,kind,outer=None):
 d=s[(s.series_id==t.series_id)&(s.year==t.year)&(s.zone==t.zone)&(s.status==t.status)]
 if outer is not None:d=d[d.site!=outer]
 if kind=='other_site':d=d[d.site!=t.site]
 else:d=d[(d.site==t.site)&(d.date!=t.date)]
 if d.empty:return None
 d=d.copy();d['gap']=(d.date-t.date).abs().dt.days;d=d[d.gap.eq(d.groupby('site').gap.transform('min'))]
 g=d.groupby('site').agg(total_mass=('total_mass','sum'),n_sized=('n_sized','sum'));pred=(g.total_mass/g.n_sized).mean()
 return dict(predicted=pred,max_gap=int(d.gap.max()),donors=';'.join(sorted(d.site.unique())))
# Recompute inner predictions without the outer target site anywhere in training.
cal={};calrows=[]
for (sid,outer,status),_ in v.groupby(['series_id','site','status']):
 for kind in ['other_site','other_visit']:
  records=[]
  for t in s[(s.series_id==sid)&(s.status==status)&(s.site!=outer)].itertuples():
   pred=predict(t,kind,outer)
   if pred:
    row=dict(series_id=sid,outer_site=outer,status=status,kind=kind,cal_site=t.site,cal_year=t.year,cal_zone=t.zone,cal_date=str(t.date.date()),observed=t.mean_mass,predicted=pred['predicted'],residual=np.log(t.mean_mass/pred['predicted']),gap_days=pred['max_gap'],donor_sites=pred['donors'],cluster=str(t.site)+'|'+str(t.year));records.append(row);calrows.append(row)
  z=pd.DataFrame(records)
  if len(z):
   assert not z.cal_site.eq(outer).any();assert all(outer not in d.split(';') for d in z.donor_sites)
   cal[(sid,outer,status,kind)]=z
pd.DataFrame(calrows).to_csv(O/'outer_site_calibration.csv',index=False)
def weighted_quantile(z,probs):
 n=z.groupby('cluster').size();weights=1/z.cluster.map(n).to_numpy();order=np.argsort(z.residual.to_numpy());r=z.residual.to_numpy()[order];w=weights[order];q=(np.cumsum(w)-.5*w)/w.sum();return np.interp(probs,q,r)
def calibrated(t,kind):
 pred=predict(t,kind);z=cal.get((t.series_id,t.site,t.status,kind))
 if pred is None:return None,'no donors in same role/zone/year'
 if z is None or z.cluster.nunique()<2:return None,'fewer than two calibration site-years'
 if pred['max_gap']>z.gap_days.max():return None,'date gap exceeds calibration support'
 q=weighted_quantile(z,[.025,.5,.975]);return dict(**pred,point=pred['predicted']*np.exp(q[1]),low=pred['predicted']*np.exp(q[0]),high=pred['predicted']*np.exp(q[2]),cal_clusters=z.cluster.nunique(),cal_samples=len(z)),''
# Honest outer-site interval assessment against known local, positive-count samples.
targets=s.merge(v[['series_id','site','year','zone','date','density']],on=['series_id','site','year','zone','date'],validate='one_to_one');targets=targets[targets.density>0]
cv=[]
for t in targets.itertuples():
 for kind in ['other_site','other_visit']:
  pred,reason=calibrated(t,kind);row=dict(series_id=t.series_id,Taxa=t.Taxa,site=t.site,status=t.status,year=t.year,zone=t.zone,date=str(t.date.date()),target_n=t.n_sized,observed=t.mean_mass,kind=kind,available=pred is not None,reason=reason)
  if pred:row.update(pred,covered=pred['low']<=t.mean_mass<=pred['high'])
  cv.append(row)
pd.DataFrame(cv).to_csv(O/'heldout_interval_validation.csv',index=False)
# Frequency-aware local-size bootstrap. One raw sample retains its draws across uses.
r=pd.read_csv('/Users/adrianstier/sbc-donham-kelp-mpa-data/data/ALL_sizefreq_2024.csv',low_memory=False);r=r[r.campus.isin(['VRG','UCSB'])&r.classcode.isin(['MESFRA','STRPUR','PANINT'])].copy();r['Taxa']=r.classcode.map({'MESFRA':'M. franciscanus','STRPUR':'S. purpuratus','PANINT':'P. interruptus'});r=r[r['count'].gt(0)&r['size'].gt(0)&(r.Taxa.eq('P. interruptus')|r['size'].ge(25))];r['date']=pd.to_datetime(r[['year','month','day']],errors='coerce');r['mass']=np.select([r.Taxa.eq('M. franciscanus'),r.Taxa.eq('S. purpuratus')],[.00059*r['size']**2.917,.00059*r['size']**2.870],default=.001352821*r['size']**2.913963113)
physical=['site','year','zone','date','Taxa'];r=r.merge(s[physical].drop_duplicates(),on=physical,validate='many_to_one');local={};means={}
for key,g in r.groupby(physical):
 w=rng.gamma(g['count'].to_numpy()[:,None],1,size=(len(g),N));local[key]=(g.mass.to_numpy()@w)/w.sum(axis=0);means[key]=np.average(g.mass,weights=g['count'])
# Shared empirical error within target site/species/year. Calibration clusters
# receive equal weight; each selected cluster contributes one observed error.
errcache={}
def error_draws(t,kind):
 key=(t.series_id,t.site,t.year,kind)
 if key not in errcache:
  z=cal[(t.series_id,t.site,t.status,kind)];clusters=[g.residual.to_numpy() for _,g in z.groupby('cluster')];ii=rng.integers(len(clusters),size=N);errcache[key]=np.array([rng.choice(clusters[j]) for j in ii])
 return errcache[key]
ledger=[];draws=[];points=[]
for t in v.itertuples():
 tx=t.series_id.split('|')[1];key=(t.site,t.year,t.zone,t.date,tx);row=dict(series_id=t.series_id,site=t.site,status=t.status,year=t.year,zone=t.zone,date=str(t.date.date()),density=t.density)
 if t.density==0:row.update(method='zero_count',mass=0.,reason='');d=np.zeros(N);point=0.
 elif key in local:
  row.update(method='measured_same_visit',mass=means[key],reason='');d=t.density*local[key];point=t.density*means[key]
 else:
  pred,reason=calibrated(t,'other_visit');kind='other_visit'
  if pred is None:pred,reason=calibrated(t,'other_site');kind='other_site'
  if pred:
   row.update(method='borrowed_'+kind,mass=pred['point'],reason='',**{k:val for k,val in pred.items() if k!='point'});point=t.density*pred['point'];d=t.density*pred['predicted']*np.exp(error_draws(t,kind))
  else:row.update(method='unsupported',mass=np.nan,reason=reason);point=np.nan;d=np.full(N,np.nan)
 row['biomass_point']=point;ledger.append(row);draws.append(d);points.append(point)
ld=pd.DataFrame(ledger);ld['row_id']=np.arange(len(ld));ld.to_csv(O/'visit_biomass_support.csv',index=False);draws=np.array(draws)
checks=ld[['series_id','site','year','zone','date','method','density']].copy()
for d in range(3):checks['draw_'+str(d)]=draws[:,d]
checks.to_csv(O/'private/visit_draw_checks.csv',index=False)
# No dropping unsized visits or stations within a role; carry missingness to paired year.
annual=[];effects=[];saved=[];coverage=[]
for sid,g in ld.groupby('series_id'):
 role={};rolepoint={}
 for (year,status),z in g.groupby(['year','status']):
  station=[];stationp=[]
  for site,ss in z.groupby('site'):
   zone=[];zonep=[]
   for _,zz in ss.groupby('zone'):
    zone.append(draws[zz.row_id].mean(axis=0));zonep.append(zz.biomass_point.to_numpy().mean())
   station.append(np.mean(zone,axis=0));stationp.append(np.mean(zonep))
  role[(year,status)]=np.mean(station,axis=0);rolepoint[(year,status)]=np.mean(stationp)
 base=a[a.series_id==sid].sort_values('year');zrows=[];yd=[]
 for t in base.itertuples():
  pp=rolepoint[(t.year,'mpa')];rr=rolepoint[(t.year,'reference')];new=np.log((pp+off[sid])/(rr+off[sid]));y=np.log((role[(t.year,'mpa')]+off[sid])/(role[(t.year,'reference')]+off[sid]));ok=np.isfinite(new) and np.isfinite(y).all()
  annual.append(dict(series_id=sid,year=t.year,t=t.t,protected_biomass=pp,reference_biomass=rr,biomass_lnratio=new,local_only_lnratio=t.revised,original_lnratio=t.original,supported=ok))
  if ok and t.t>0:zrows.append(dict(year=t.year,t=t.t,new=new,old=t.original));yd.append(y)
 post=base[base.t>0];coverage.append(dict(series_id=sid,original_post_years=len(post),local_only_post_years=int(post.revised.notna().sum()),borrowing_post_years=len(zrows),unsupported_years=';'.join(map(str,sorted(set(post.year)-{q['year'] for q in zrows})))))
 if len(zrows)>=3:
  z=pd.DataFrame(zrows);fitnew=slope(z,'new');fitold=slope(z,'old');bb=(np.linalg.pinv(np.column_stack([np.ones(len(z)),z.t]))@np.array(yd))[1];effects.append(dict(series_id=sid,current_same_year_slope=fitold['estimate'],**fitnew,prediction_draw_low=np.quantile(bb,.025),prediction_draw_high=np.quantile(bb,.975),draw_mean=bb.mean(),draws=N))
  for d in range(3):
   for row,val in zip(zrows,np.array(yd)[:,d]):saved.append(dict(series_id=sid,draw=d,year=row['year'],t=row['t'],lnDiff=val,expected_slope=bb[d]))
pd.DataFrame(annual).to_csv(O/'biomass_annual.csv',index=False);pd.DataFrame(coverage).to_csv(O/'biomass_coverage.csv',index=False);pd.DataFrame(effects).to_csv(O/'biomass_effects.csv',index=False);pd.DataFrame(saved).to_csv(O/'private/verification_draws.csv',index=False)
policy=dict(scope='Six affected PISCO animal comparisons; existing retained sites/years',density='All valid months, equal visits, constant shared INNER/MID zones',biomass='Review sensitivity only: measured sizes, then calibrated same-role/zone borrowing; unsupported remains missing',uncertainty='Empirical site-year prediction error plus local size resampling; fixed counts, calibration equations and offsets; not total confidence intervals',paper_update=False,draws=N,seed=20261102)
(O/'review_analysis_policy.json').write_text(json.dumps(policy,indent=2))
protected=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text());assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in protected.items());info=dict(density_post_years=63,local_only_biomass_post_years=25,borrow_supported_post_years=sum(q['borrowing_post_years'] for q in coverage),conditional_biomass_fits=len(effects),draws=N,outer_validation_attempts=len(cv),calibration_rows=len(calrows),protected_unchanged=len(protected));(O/'verification.json').write_text(json.dumps(info,indent=2));print(json.dumps(info,indent=2));print(pd.DataFrame(coverage).to_string(index=False));print(ld.groupby('method').size())
