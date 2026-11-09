#!/usr/bin/env python3
"""Independent NumPy regression and calendar-HAC reconstruction; no primary writes."""
from pathlib import Path
import hashlib,json,sys
import numpy as np
import pandas as pd
from scipy.stats import t as tdist
out=Path(sys.argv[1] if len(sys.argv)>1 else 'outputs/satellite_baseline_comparison_20261101')
cohort=pd.read_csv(out/'cohort.csv')
if (out/'freeze_manifest.json').exists():
 for p,h in json.loads((out/'freeze_manifest.json').read_text())['files'].items():
  assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h, p
d=pd.read_csv(out/'model_input_years.csv');r=pd.read_csv(out/'estimates.csv');c=pd.read_csv(out/'coefficients.csv')
maxerr=0
for (site,window),x in d.groupby(['MPA','window']):
 x=x.sort_values('year');n=len(x);yr=x.year.to_numpy();post=x.post.to_numpy();tt=x.t.to_numpy();pt=x.post_time.to_numpy();y=x.lnDiff.to_numpy()
 for model in ['Stable','Changing']:
  X=np.column_stack([np.ones(n),post,pt]) if model=='Stable' else np.column_stack([np.ones(n),tt,post,pt])
  names=['(Intercept)','post','post_time'] if model=='Stable' else ['(Intercept)','t','post','post_time']
  b=np.linalg.lstsq(X,y,rcond=None)[0];res=y-X@b;S=X*res[:,None];bread=np.linalg.inv(X.T@X);k=X.shape[1]
  cc=c[(c.MPA==site)&(c.window==window)&(c.model==model)].set_index('term').loc[names].coefficient.to_numpy();assert np.allclose(b,cc,atol=1e-9)
  # Separate before/after fits must reproduce the fully interacted OLS model.
  before=post==0;after=post==1
  xb=np.column_stack([np.ones(before.sum()),tt[before]]) if model=='Changing' else np.ones((before.sum(),1))
  bb=np.linalg.lstsq(xb,y[before],rcond=None)[0]
  assert np.allclose(bb,b[:2] if model=='Changing' else b[:1])
  if model=='Stable':
   gmean=np.zeros(k);gmean[-2:]=[1,pt[after].mean()]
   assert np.allclose(gmean@b,y[after].mean()-y[before].mean())
  xa=np.column_stack([np.ones(after.sum()),tt[after]]);ba=np.linalg.lstsq(xa,y[after],rcond=None)[0]
  assert np.allclose(xa@ba,(X@b)[after])
  for lag in [1,3,5]:
   weights=np.maximum(0,1-np.abs(yr[:,None]-yr[None,:])/(lag+1))
   V=bread@(S.T@weights@S)@bread*n/(n-k)
   rr=r[(r.MPA==site)&(r.window==window)&(r.model==model)&(r.lag==lag)]
   for _,row in rr.iterrows():
    g=np.zeros(k)
    if row.quantity=='Mean_after':g[-2:]=[1,pt[after].mean()]
    if row.quantity=='Last_after':g[-2:]=[1,pt[-1]]
    if row.quantity=='Level':g[-2]=1
    if row.quantity=='Slope':g[-1]=1
    est=g@b;se=np.sqrt(g@V@g);lo=est-tdist.ppf(.975,n-k)*se;hi=est+tdist.ppf(.975,n-k)*se
    err=np.max(np.abs(np.array([est,se,lo,hi])-row[['estimate','se','low','high']].to_numpy(float)));maxerr=max(maxerr,err)
    assert err<1e-7,(site,window,model,lag,row.quantity,err)
p=pd.read_csv(out/'predictions.csv');post=p[p.post==1].pivot(index=['MPA','year'],columns='model',values='fitted');assert np.allclose(post.Stable,post.Changing)
assert len(r)==len(cohort)*3*2*3*4
hashes=json.loads((out/'protected_hashes.json').read_text());changed=[p for p,h in hashes.items() if hashlib.sha256(Path(p).read_bytes()).hexdigest()!=h];assert not changed,changed
result={'contrasts_checked':len(r),'max_absolute_error':maxerr,'protected_files_unchanged':len(hashes),'split_regression_check':'pass','same_post_fitted_values':'pass','same_model_years':f'shared input per site/window, {len(cohort)*3} inputs'}
(out/'verification.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
