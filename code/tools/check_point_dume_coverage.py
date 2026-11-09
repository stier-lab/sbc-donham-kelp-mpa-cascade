from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
from scipy.stats import t
out=Path('outputs/point_dume_coverage_20261101');e=pd.read_csv(out/'effects.csv');p=pd.read_csv(out/'pooled.csv');original=pd.read_csv('outputs/full_followup_recovery_20261101/effects.csv').query("scenario=='Full'").set_index(['metric','series_id']);err=[]
assert not e.duplicated(['scenario','metric','series_id']).any()
for _,r in e.iterrows():
 a=original.loc[(r.metric,r.series_id)];assert abs(a.estimate-r.estimate)<1e-12 and abs(a.vi-r.vi)<1e-12
for _,r in p.iterrows():
 x=e[(e.scenario==r.scenario)&(e.metric==r.metric)&(e.Taxa==r.Taxa)&(e.Resp==r.Resp)]
 if np.isfinite(r.sigma_mpa):
  mp=x.MPA.to_numpy();V=np.diag(x.vi.to_numpy())+r.sigma_mpa*(mp[:,None]==mp[None,:])
  if np.isfinite(r.sigma_series):V+=np.eye(len(x))*r.sigma_series
  w=np.linalg.solve(V,np.ones(len(x)));est=w@x.estimate.to_numpy()/w.sum();er=abs(est-r.estimate);assert er<1e-7;err.append(er)
 if np.isfinite(r.low):assert abs(r.low-(r.estimate-t.ppf(.975,r.df)*r.se))<1e-7
for _,z in p.groupby(['scenario','metric']):
 z=z[z.p.notna()].sort_values('p');q=np.minimum(1,np.minimum.accumulate((z.p.to_numpy()*9/np.arange(1,len(z)+1))[::-1])[::-1]);assert np.allclose(q,z.q)
c=pd.read_csv(out/'coverage_categories.csv');assert len(c)==169 and c.series_id.nunique()==169
for _,r in pd.read_csv(out/'rule_membership.csv').iterrows():
 n=int(r.rule.split('_')[0][1:]);span=int(r.rule.split('_')[1][4:]);early=r.rule.split('_')[2][5:];expect=(r.na>=n and r.span>=span and (early=='any' or r.first_age<=int(early)));assert expect==r.included
ba=e.query("scenario=='All_three_year' and metric=='Change'").set_index('series_id');new=e.query("scenario=='n5_span5_early3' and metric=='Change'").set_index('series_id');pd.testing.assert_frame_equal(ba[['estimate','vi']].sort_index(),new[['estimate','vi']].sort_index())
f=pd.read_csv(out/'footprint_ratios.csv');fs=pd.read_csv(out/'footprint_slopes.csv');ferr=[]
for _,r in fs.iterrows():
 z=f[(f.scenario==r.scenario)&(f.Taxa==r.Taxa)];a=np.polyfit(z.year-z.year.mean(),z.ln_ratio,1)[0];ferr.append(abs(a-r.slope));assert abs(a-r.slope)<1e-10
h=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text())
for fn,digest in h.items():assert hashlib.sha256(Path(fn).read_bytes()).hexdigest()==digest
v={'inherited_effect_rows_verified':len(e),'pooled_points_verified':len(err),'max_GLS_error':max(err),'footprint_slopes_verified':len(ferr),'max_slope_error':max(ferr),'before_after_55_unchanged':True,'coverage_rules_verified':True,'protected_files_unchanged':len(h)};(out/'verification.json').write_text(json.dumps(v,indent=2));print(v)
