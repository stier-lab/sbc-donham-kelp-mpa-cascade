from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
from scipy.stats import t
out=Path('outputs/full_followup_recovery_20261101');e=pd.read_csv(out/'effects.csv');o=pd.read_csv(out/'model_rows.csv');p=pd.read_csv(out/'pooled.csv');errors=[];errs=[]
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
for metric in e.metric.unique():
 a=e[(e.scenario=='Nine_year_matched')&(e.metric==metric)];b=e[(e.scenario=='Full_matched')&(e.metric==metric)];assert set(a.series_id)==set(b.series_id)
for scenario in e.scenario.unique():
 def ids(mt):return set(e[(e.scenario==scenario)&(e.metric==mt)].series_id)
 assert not ids('BA_post_slope')&ids('No_baseline_post_slope')
 assert ids('BA_post_slope')|ids('No_baseline_post_slope')==ids('Post_slope')
h=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text())
for f,digest in h.items():assert hashlib.sha256(Path(f).read_bytes()).hexdigest()==digest
result={'individual_effects_checked':len(errors),'max_individual_error':max(errors),'pooled_points_checked':len(errs),'max_pooled_error':max(errs),'protected_files_unchanged':len(h),'matched_cohorts_and_disjoint_groups':True,'BH_and_interval_arithmetic':True}
(out/'verification.json').write_text(json.dumps(result,indent=2));print(result)
