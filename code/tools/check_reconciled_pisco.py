from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
from scipy.stats import t
out=Path('outputs/pisco_source_reconciliation_20261101');e=pd.read_csv(out/'effects.csv');o=pd.read_csv(out/'model_rows.csv');p=pd.read_csv(out/'pooled.csv');errors=[];errs=[]
assert not e.duplicated(['scenario','metric','series_id']).any()
assert not o.duplicated(['scenario','metric','series_id','year']).any()
lookup=e.set_index(['scenario','metric','series_id'])
for key,d in o.groupby(['scenario','metric','series_id']):
 scenario,metric,id=key;r=lookup.loc[key];d=d.sort_values('year');yr=d.year.to_numpy();tm=d.t.to_numpy();post=d.post.to_numpy();y=d.lnDiff.to_numpy();ones=np.ones(len(d))
 if metric in ['Post_slope','BA_post_slope','No_baseline_post_slope']:X=np.column_stack([ones,tm]);g=np.array([0,1])
 elif metric=='Change':X=np.column_stack([ones,post,post*tm]);g=np.array([0,1,tm[post==1].mean()])
 else:X=np.column_stack([ones,tm,post,post*tm]);g=np.array([0,0,0,1]) if metric=='Slope_change' else np.array([0,0,1,tm[post==1].mean()])
 b=np.linalg.lstsq(X,y,rcond=None)[0];u=y-X@b;B=np.linalg.inv(X.T@X);n,k=X.shape;lag=min(int(d.lag.iloc[0]),int(yr.max()-yr.min()));W=np.maximum(0,1-abs(yr[:,None]-yr[None,:])/(lag+1));S=X*u[:,None];V=B@(S.T@W@S)@B*n/(n-k)
 err=max(abs(g@b-r.estimate),abs(g@V@g-r.vi));assert err<1e-7,(key,err);errors.append(err)
 assert (tm!=0).all()
 if 'Nine_year' in scenario:assert (tm[post==1]<=9).all()
for _,r in p.iterrows():
 x=e[(e.scenario==r.scenario)&(e.metric==r.metric)&(e.Taxa==r.Taxa)&(e.Resp==r.Resp)]
 if np.isfinite(r.sigma_mpa):
  mp=x.MPA.to_numpy();V=np.diag(x.vi.to_numpy())+r.sigma_mpa*(mp[:,None]==mp[None,:]);
  if np.isfinite(r.sigma_series):V+=np.eye(len(x))*r.sigma_series
  w=np.linalg.solve(V,np.ones(len(x)));est=w@x.estimate.to_numpy()/w.sum();err=abs(est-r.estimate);assert err<1e-7;errs.append(err)
 if np.isfinite(r.low):assert abs(r.low-(r.estimate-t.ppf(.975,r.df)*r.se))<1e-7
for _,z in p.groupby(['scenario','metric']):
 z=z[z.p.notna()].sort_values('p');q=np.minimum(1,np.minimum.accumulate((z.p.to_numpy()*9/np.arange(1,len(z)+1))[::-1])[::-1]);assert np.allclose(q,z.q)
a=pd.read_csv(out/'annual_input.csv');base=a[a.scenario=='Previous_review'].set_index(['series_id','year']);targets=set(pd.read_csv(out/'reconstruction.csv').series_id)
for sc in ['Complete_zones','Fixed_valid_zones','Supported_corrections']:
 b=a[a.scenario==sc].set_index(['series_id','year']);assert b.index.equals(base.index);other=~b.index.get_level_values('series_id').isin(targets);assert np.array_equal(b.loc[other,'lnDiff'],base.loc[other,'lnDiff'])
 for mt in e.metric.unique():assert set(e[(e.scenario==sc)&(e.metric==mt)].series_id)==set(e[(e.scenario=='Previous_review')&(e.metric==mt)].series_id)
 before=e[(e.scenario=='Previous_review')&(e.metric=='Change')].set_index('series_id');after=e[(e.scenario==sc)&(e.metric=='Change')].set_index('series_id');assert np.allclose(before.estimate.sort_index(),after.estimate.sort_index()) and np.allclose(before.vi.sort_index(),after.vi.sort_index())
z=pd.read_csv(out/'scenario_zone_values.csv');an=pd.read_csv(out/'source_anomalies.csv');roles=pd.read_csv(out/'role_values.csv');patch=pd.read_csv(out/'patches.csv');co=pd.read_csv(out/'reconstruction.csv').set_index('series_id')
for (sid,sc),g in z.groupby(['series_id','scenario']):
 if sc!='Legacy_reproduced':
  assert g.valid.all()
  bad=an[an.series_id==sid].iloc[0];assert not ((g.site==bad.site)&(g.year==bad.year)&(g.zone==bad.zone)).any()
 if sc=='Fixed_valid_zones':assert len(set(tuple(sorted(q.site+'|'+q.zone)) for _,q in g.groupby('year')))==1
 # Independent year-role reconstruction, explicit equal station then zone weighting.
 for (year,status),gg in g.groupby(['year','status']):
  value=np.mean([np.mean(q.value.to_numpy()) for _,q in gg.groupby('site')]);r=roles[(roles.series_id==sid)&(roles.scenario==sc)&(roles.year==year)&(roles.status==status)].iloc[0];assert abs(r.value-value)<1e-9
 rr=roles[(roles.series_id==sid)&(roles.scenario==sc)].pivot(index='year',columns='status',values='value');offset=co.loc[sid,'offset'];pp=patch[(patch.series_id==sid)&(patch.scenario==sc)].set_index('year');assert np.allclose(np.log((rr.mpa+offset)/(rr.reference+offset)),pp.loc[rr.index,'lnDiff'])
# Frozen size-based companion history already excludes the flawed Scorpion zone.
f=pd.read_csv('outputs/seasonal_zone_audit_20261101/fixed_footprints.csv');f=f[(f.series_id=='Scorpion SMR|M. franciscanus|PISCO|Bio')&(f.scenario=='Fixed_station_zones')];assert not (f.site=='SCI_CAVERN_POINT_W').any()
h=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text())
for path,digest in h.items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest
for path,digest in json.loads((out/'source_hashes.json').read_text()).items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest
# Preserve every prior review artifact, including its scripts (no overwritten comparison).
for path,digest in json.loads(Path('outputs/seasonal_zone_audit_20261101/manifest.json').read_text()).items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest
r={'effects_checked':len(errors),'max_effect_error':max(errors),'pooled_points_checked':len(errs),'max_GLS_error':max(errs),'legacy_role_year_values_verified':int(co.raw_rows.sum()),'source_hashes_and_previous_review_unchanged':True,'protected_files_unchanged':len(h),'same_model_membership_and_years':True,'55_before_after_histories_unchanged':True,'missing_units_excluded_without_imputation':True,'Scorpion_biomass_companion_already_excludes_bad_station':True}
(out/'verification.json').write_text(json.dumps(r,indent=2));print(r)
dep=pd.read_csv(out/'point_vicente_2013_depth_audit.csv');assert set(dep.loc[dep.season_retained,'zone'])=={'DEEP'};assert set(dep.loc[~dep.season_retained,'month'])=={4}
v=z[(z.series_id=='Point Vicente SMCA|M. pyrifera|PISCO|Bio')&(z.scenario=='Fixed_valid_zones')&(z.status=='mpa')];assert set(v.zone)=={'DEEP'} and (v[v.year>=2017].value==0).all()
target='Point Vicente SMCA|M. pyrifera|PISCO|Bio';b=a[(a.series_id==target)&(a.scenario=='Supported_corrections')].set_index('year');c=a[(a.series_id==target)&(a.scenario=='Complete_zones')].set_index('year');assert np.array_equal(b.lnDiff,c.lnDiff)
r['DEEP_only_Point_Vicente_withheld_from_supported_corrections']=True;(out/'verification.json').write_text(json.dumps(r,indent=2))
