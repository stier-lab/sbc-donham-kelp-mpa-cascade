#!/usr/bin/env python3
from pathlib import Path
import numpy as np,pandas as pd,json,hashlib
from scipy.stats import t,false_discovery_control
import sys
out=Path(sys.argv[1] if len(sys.argv)>1 else 'outputs/pooled_common_estimand_20261101')
for p,h in json.loads((out/'freeze_manifest.json').read_text()).items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
rows=pd.read_csv(out/'model_rows.csv');e=pd.read_csv(out/'effects.csv');lookup=e.set_index(['series_id','scenario','metric']);maxerr=0
for key,x in rows.groupby(['series_id','scenario','metric'],sort=False):
 x=x.sort_values('year');y=x.lnDiff.to_numpy();tt=x.t.to_numpy();post=x.post.to_numpy();pt=x.pt.to_numpy();n=len(x);metric=key[2];scenario=key[1]
 if metric=='Post_difference':X=np.column_stack([np.ones(n),tt]);g=X.mean(axis=0)
 elif scenario=='Trend':X=np.column_stack([np.ones(n),tt,post,pt]);g=np.array([0,0,1,pt[post==1].mean()])
 else:X=np.column_stack([np.ones(n),post,pt]);g=np.array([0,1,pt[post==1].mean()])
 b=np.linalg.lstsq(X,y,rcond=None)[0];u=y-X@b;bread=np.linalg.inv(X.T@X);S=X*u[:,None];yr=x.year.to_numpy();lag=min(x.lag.iloc[0],yr.max()-yr.min())
 W=np.maximum(0,1-abs(yr[:,None]-yr[None,:])/(lag+1));V=bread@(S.T@W@S)@bread*n/(n-X.shape[1])
 est=g@b;se=np.sqrt(g@V@g);r=lookup.loc[key];err=max(abs(est-r.estimate),abs(se-r.se));maxerr=max(maxerr,err);assert err<1e-7,(key,err)
 if metric=='Post_difference':assert np.isclose(est,y.mean())
 elif scenario!='Trend':assert np.isclose(est,y[post==1].mean()-y[post==0].mean())
protected=json.loads((out/'protected_hashes.json').read_text());assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in protected.items())
ans=dict(individual_effects_checked=len(e),max_abs_error=maxerr,protected_unchanged=len(protected))
if (out/'pooled.csv').exists():
 p=pd.read_csv(out/'pooled.csv');me=0
 assert not p.duplicated(['scenario','metric','Taxa','Resp','view']).any()
 for _,r in p.iterrows():
  if r.k<2 or r.n_mpa<2:continue
  x=e[(e.scenario==r.scenario)&(e.metric==r.metric)&(e.Taxa==r.Taxa)&(e.Resp==r.Resp)]
  if r['view']=='Divers':x=x[x.Source!='Landsat']
  if r['view']=='Satellite':x=x[x.Source=='Landsat']
  V=np.diag(x.vi.to_numpy())+r.sigma_mpa*(x.MPA.to_numpy()[:,None]==x.MPA.to_numpy()[None,:])
  if pd.notna(r.sigma_series):V+=r.sigma_series*np.eye(len(x))
  w=np.linalg.solve(V,np.ones(len(x)));est=w@x.estimate/w.sum();me=max(me,abs(est-r.estimate));assert np.isclose(est,r.estimate,atol=1e-7)
  if pd.notna(r.low):
   assert np.isclose(r.low,r.estimate-t.ppf(.975,r.df)*r.se)
   assert np.isclose(r.p,2*t.sf(abs(r.estimate/r.se),r.df))
 for _,g in p.groupby(['scenario','metric','view']):
  x=g.p.dropna();adj=false_discovery_control(np.r_[x.to_numpy(),np.ones(max(0,9-len(x)))])[:len(x)]
  assert np.allclose(adj,g.loc[x.index,'q'])
 ans.update(pooled_rows=len(p),max_pooled_GLS_error=me)
(out/'verification.json').write_text(json.dumps(ans,indent=2));print(json.dumps(ans,indent=2))
