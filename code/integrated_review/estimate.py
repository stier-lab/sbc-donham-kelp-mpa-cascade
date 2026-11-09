"""Exact measurement draws + calendar-time working covariance, separate estimands.
Measurement covariance is propagated jointly. Residual process covariance is
estimated AFTER accounting for that measurement covariance (no variance addition
to an already total empirical standard error). Shared references are retained
through the same physical measurement draws and reported covariance matrix.
"""
from common import *
import sys
sys.path.insert(0,str(ROOT/'code/tools'))
from linear_mixed import fit
from scipy.stats import t as student,shapiro
B=POLICY.get('measurement_draws',1200)
x=pd.read_parquet(PRIVATE/'supported_measurement_uses.parquet');roster=read_tracked(OUT/'comparison_year_roster.csv');v=read_tracked(OUT/'visit_biomass_support.csv');p=read_tracked(OUT/'size_pool_inventory.csv');mem=pd.read_parquet(PRIVATE/'size_pool_members.parquet')
off=read_tracked(ROOT/'outputs/all_program_station_refits_20261102/offsets.csv').set_index('series_id').offset
if VARIANT in ['constant_stations','constant_zones']:
 keep=[]
 for series,g in x.groupby('series_id'):
  if VARIANT=='constant_stations':
   all_years=set(g.year);sets=g.groupby('site').year.agg(set);sites=set(sets[sets.map(lambda s:s==all_years)].index);z=g[g.site.isin(sites)]
  else:
   sets=g.groupby('year').zone.agg(set);zones=set.intersection(*sets);z=g[g.zone.isin(zones)]
  valid=z.groupby('year').status.agg(set);years=set(valid[valid.map(lambda s:s=={'mpa','reference'})].index);keep.append(z[z.year.isin(years)])
 x=pd.concat(keep,ignore_index=True)
pool_cv={}
if VARIANT=='sparse_size_variance':
 allowed=x[['Source','site','year','Taxa']].drop_duplicates();used=p.merge(allowed,on=['Source','site','year','Taxa'],how='inner');cvrows=[]
 for rr in used.itertuples():
  gg=mem[mem.pool_id.eq(rr.pool_id)];nn=gg.frequency.sum();mm=rr.mean_mass
  if nn>1 and len(gg)>1 and mm>0:cvrows.append(dict(Source=rr.Source,Taxa=rr.Taxa,ss=float((gg.frequency*(gg.mass/mm-1)**2).sum()),df=float(nn-1)))
 cv=pd.DataFrame(cvrows).groupby(['Source','Taxa'])[['ss','df']].sum();pool_cv=(cv.ss/cv.df).to_dict()
weight_cache={};pool_cache={};mu=p.set_index('pool_id').mean_mass.to_dict();visit_index=v.set_index(['series_id','site','year','date','zone'])
def random(key):return np.random.default_rng(int(hashlib.sha256((('second_seed|' if VARIANT=='second_seed' else '')+str(key)).encode()).hexdigest()[:15],16))
def weights(key):
 if key not in weight_cache:weight_cache[key]=random(key).exponential(size=B)
 return weight_cache[key]
def pool_draw(key):
 if key not in pool_cache:
  g=mem[mem.pool_id.eq(key)];w=np.array([random(r.source_record_id).gamma(shape=r.frequency,scale=1,size=B) for r in g.itertuples()]);pool_cache[key]=g.mass.to_numpy()@w/w.sum(axis=0)
  if VARIANT=='sparse_size_variance':
   pp=p[p.pool_id.eq(key)].iloc[0];n=pp.n_sized;cv=pool_cv.get((pp.Source,pp.Taxa),np.nan)
   if np.isfinite(cv):
    extra=2*cv/((max(n-1,0)+2)*max(n,1));lv=np.log1p(extra);pool_cache[key]*=np.exp(random('sparse|'+key).normal(-.5*lv,np.sqrt(lv),B))
 return pool_cache[key]
def aggregate_units(g,column):
 # Average quadrats within LTER transects before the station/visit average.
 nested=g.unit.str.split('|').str[0] if g.Source.iloc[0]=='LTER' else pd.Series('all',index=g.index)
 values=[];draws=[];outer_weights=[]
 for nest,q in g.groupby(nested):
  a=q[column].to_numpy(float)
  if not np.isfinite(a).all():return np.nan,np.full(B,np.nan)
  w=np.array([weights(k) for k in q.unit_id]);d=np.tile(a[:,None],(1,B))
  if column=='biomass':
   for j,r in enumerate(q.itertuples()):
    if r.imputed_density>0:
     d[j]+=r.imputed_density*(pool_draw(r.imputed_pool)-mu[r.imputed_pool])
    if isinstance(r.donor_pool,str) and r.donor_pool:
     keys=json.loads(r.donor_pool);base=np.mean([mu[k] for k in keys]);dd=np.mean([pool_draw(k) for k in keys],axis=0)
     if base>0:d[j]*=dd/base
     vv=visit_index.loc[(r.series_id,r.site,r.year,r.date,r.zone)]
     if vv.prediction_error_sd>0 or vv.prediction_shared_sd>0:
      key=vv.prediction_error_key;sd=vv.prediction_error_sd
      shared=float(vv.prediction_shared_sd);d[j]*=np.exp(random(key).normal(-.5*shared*shared,shared,size=B)+random(vv.prediction_target_key).normal(-.5*sd*sd,sd,size=B))
  values.append(a.mean());draws.append(np.sum(w*d,axis=0)/w.sum(axis=0));outer_weights.append(weights('LTER_permanent_transect|'+str(q.site.iloc[0])+'|'+str(nest)) if q.Source.iloc[0]=='LTER' else np.ones(B))
 outer_weights=np.array(outer_weights);return np.mean(values),np.sum(outer_weights*np.array(draws),axis=0)/outer_weights.sum(axis=0)

annual=[];annual_draws={};station=[]
for (series,yr),g in x.groupby(['series_id','year']):
 col='density' if g.Resp.iloc[0]=='Den' else 'biomass';roles={};rd={}
 for role,z in g.groupby('status'):
  sites=[];sd=[]
  for site,zz in z.groupby('site'):
   zones=[];zd=[]
   for zone,zzz in zz.groupby('zone'):
    dates=[];dd=[]
    # Method is part of a PISCO survey event; repeated visits receive equal weight.
    zz4=zzz.copy();zz4['_method']=zz4.unit.str.split('|').str[0] if zz4.Source.iloc[0]=='PISCO' else 'recorded_method'
    for _,ev in zz4.groupby(['date','_method']):
     value,dr=aggregate_units(ev,col);dates.append(value);dd.append(dr)
    zones.append(np.mean(dates));zd.append(np.mean(dd,axis=0))
   value=np.mean(zones);dr=np.mean(zd,axis=0);sites.append(value);sd.append(dr)
   station.append(dict(series_id=series,year=yr,site=site,status=role,value=value,units=len(zz),visits=zz.date.nunique(),zones=';'.join(sorted(zz.zone.unique()))))
  roles[role]=np.mean(sites);rd[role]=np.mean(sd,axis=0)
 meta=roster[(roster.series_id==series)&roster.year.eq(yr)].iloc[0].to_dict();meta.update(protected=roles.get('mpa',np.nan),reference=roles.get('reference',np.nan));annual.append(meta);annual_draws[(series,yr)]=rd
an=pd.DataFrame(annual);st=pd.DataFrame(station);st.to_csv(OUT/'station_annual_values.csv',index=False)
# Fixed legacy offsets keep source corrections separate from zero handling.
# New histories use the same half-minimum-positive rule and are labeled explicitly.
for series,g in an.groupby('series_id'):
 positives=g[['protected','reference']].to_numpy().ravel();positives=positives[np.isfinite(positives)&(positives>0)]
 o=float(off[series]) if series in off.index else .5*min(positives) if len(positives) else np.nan
 an.loc[g.index,'offset']=o;an.loc[g.index,'offset_basis']='frozen_existing' if series in off.index else 'half_minimum_positive'
 an.loc[g.index,'ln_ratio']=np.log((g.protected+o)/(g.reference+o))
an['response_supported']=np.isfinite(an.ln_ratio);an['annual_index']=np.arange(len(an));
all_draws=[]
for r in an.itertuples():
 rd=annual_draws[(r.series_id,r.year)];all_draws.append(np.log((rd['mpa']+r.offset)/(rd['reference']+r.offset)))
np.savez_compressed(PRIVATE/'annual_measurement_draws.npz',draws=np.array(all_draws));an.to_csv(OUT/'annual_values.csv',index=False)
results=[];diagnostics=[];elig=[];influence=[];coefs=[];drawcoeff=[]
for series,allg in an.groupby('series_id'):
 g=allg[allg.response_supported].sort_values('year');pre=int(g.t.lt(0).sum());post=int(g.t.ge(0).sum());ba=pre>=POLICY['minimum_before_years'];status='before_after_and_post_trend' if ba and post>=3 else 'post_trend_only' if post>=3 else 'insufficient_post_years'
 elig.append(dict(series_id=series,before_years=pre,after_years=post,eligible_model=status,eligible_four_before=pre>=4,first_year=g.year.min(),last_year=g.year.max(),years=';'.join(map(str,g.year)),missing_response_years=';'.join(map(str,allg.loc[~allg.response_supported,'year']))))
 if post<3:continue
 models=['post_trend']+(['flat_before_after','pretrend_sensitivity'] if ba else [])
 for model in models:
  h=g[g.t>=0].copy() if model=='post_trend' else g.copy();time=h.t.to_numpy(float);after=(time>=0).astype(float);recovery=np.maximum(time,0);years=h.year.to_numpy(int);y=h.ln_ratio.to_numpy(float);n=len(y)
  names=['intercept','recovery_slope'] if model=='post_trend' else ['baseline','level_change','recovery_slope'] if model=='flat_before_after' else ['baseline','pre_slope','level_change','slope_change']
  X=np.column_stack([np.ones(n),time]) if model=='post_trend' else np.column_stack([np.ones(n),after,recovery]) if model=='flat_before_after' else np.column_stack([np.ones(n),time,after,recovery])
  if np.linalg.matrix_rank(X)<X.shape[1] or n<=X.shape[1]:continue
  d=np.array([np.log((annual_draws[(series,yr)]['mpa']+o)/(annual_draws[(series,yr)]['reference']+o)) for yr,o in zip(h.year,h.offset)])
  if not np.isfinite(d).all():raise ValueError('Nonfinite measurement draws '+series)
  S=measurement_covariance(d,years,h.Source.iloc[0],h.Resp.iloc[0],h.Taxa.iloc[0])
  for rho in [0,.3,.6]:
   K=rho**np.abs(years[:,None]-years[None,:]);f=fit(y,X,S,[K]);b=f['beta'];cov=f['bc'];res=y-X@b;df=n-X.shape[1];crit=student.ppf(.975,df)
   # Keep observations in their real years; missing calendar years are not compressed.
   lag_pairs=[(i,j) for i in range(n) for j in range(i) if years[i]-years[j]==1]
   lag=float(np.corrcoef([res[i] for i,j in lag_pairs],[res[j] for i,j in lag_pairs])[0,1]) if len(lag_pairs)>2 and np.std(res)>0 else np.nan
   for j,name in enumerate(names):
    se=np.sqrt(max(cov[j,j],0));results.append(dict(series_id=series,model=model,rho=rho,term=name,estimate=b[j],se=se,low=b[j]-crit*se,high=b[j]+crit*se,n_years=n,before_years=pre,after_years=post,df=df,process_variance=f['q'][0],boundary=f['q'][0]<1e-8))
   diagnostics.append(dict(series_id=series,model=model,rho=rho,n=n,rank=np.linalg.matrix_rank(X),condition=np.linalg.cond(X),calendar_lag1=lag,process_variance=f['q'][0],boundary=f['q'][0]<1e-8,upper_boundary=f['upper_boundary'],measurement_rank=np.linalg.matrix_rank(S),zero_measurement_years=int((np.diag(S)==0).sum()),uncertainty='conditional on allometry, source support, observed units and donor calibration'))
   if rho==.3:
    Cinv=np.linalg.pinv(f['C']);H=np.linalg.solve(X.T@Cinv@X,X.T@Cinv);bd=H@(d-d.mean(axis=1,keepdims=True))
    for j,name in enumerate(names):coefs.append(dict(series_id=series,model=model,term=name,estimate=b[j],conditional_total_variance=cov[j,j]));drawcoeff.append(bd[j])
    for i in range(n):
     ix=np.arange(n)!=i
     if np.linalg.matrix_rank(X[ix])==X.shape[1]:
      beta=np.linalg.lstsq(X[ix],y[ix],rcond=None)[0];influence.append(dict(series_id=series,model=model,omitted_year=int(years[i]),unweighted_recovery_coefficient=beta[-1]))
pd.DataFrame(results).to_csv(OUT/'model_results.csv',index=False);pd.DataFrame(diagnostics).to_csv(OUT/'model_diagnostics.csv',index=False);pd.DataFrame(elig).to_csv(OUT/'model_eligibility.csv',index=False);pd.DataFrame(influence).to_csv(OUT/'year_influence.csv',index=False)
cf=pd.DataFrame(coefs);cf.to_csv(OUT/'coefficient_index.csv',index=False);np.savez_compressed(PRIVATE/'raw_measurement_coefficient_draws.npz',measurement_deviations=np.array(drawcoeff));np.savez_compressed(PRIVATE/'raw_monte_carlo_coefficient_covariance.npz',covariance=np.cov(np.array(drawcoeff)))
# Paired source change table; no comparison of different estimands is labeled an effect change.
old=read_tracked(ROOT/'outputs/pisco_fish_rules_20261101/approved_annual_input.csv');cmp=an.merge(old[['series_id','year','lnDiff']],on=['series_id','year'],how='outer',indicator=True);cmp['difference']=cmp.ln_ratio-cmp.lnDiff;cmp.to_csv(OUT/'annual_change_comparison.csv',index=False)
write_inputs();print('Annual responses',len(an),'supported',an.response_supported.sum(),'fitted coefficients',len(results))

# Unstructured finite-draw exports are diagnostics, never the fitted covariance.
for old_name in ['coefficient_draws.npz','coefficient_measurement_covariance.npz']:
 (PRIVATE/old_name).unlink(missing_ok=True)
(PRIVATE/'auxiliary_draw_export.json').write_text(json.dumps(dict(scope='Unstructured finite Monte Carlo coefficient draws for diagnostics only; not joint-model covariance',main_fitted_covariance='joint_coefficient_covariance.npz'),indent=2))
