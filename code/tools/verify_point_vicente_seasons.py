from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
from scipy.stats import t
out=Path('outputs/point_vicente_seasons_20261101');e=pd.read_csv(out/'effects.csv');o=pd.read_csv(out/'model_rows.csv');p=pd.read_csv(out/'pooled.csv');errors=[];errs=[]
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
a=pd.read_csv(out/'annual_input.csv');base=a[a.scenario=='Previous_review'].set_index(['series_id','year']);local=pd.read_csv(out/'local_annual.csv');sl=pd.read_csv(out/'local_slopes.csv');changed=set(local.series_id);assert len(changed)==4
for sc,g in a.groupby('scenario'):
 gg=g.set_index(['series_id','year']);other=~gg.index.get_level_values('series_id').isin(changed);pd.testing.assert_series_equal(gg.loc[other,'lnDiff'].sort_index(),base.loc[~base.index.get_level_values('series_id').isin(changed),'lnDiff'].sort_index())
 for mt in e.metric.unique():assert set(e[(e.scenario==sc)&(e.metric==mt)].series_id)==set(e[(e.scenario=='Previous_review')&(e.metric==mt)].series_id)
 if sc=='Previous_review':continue
 for sid,h in g[g.series_id.isin(changed)].groupby('series_id'):
  prev=base.loc[sid];added=set(h.year)-set(prev.index);removed=set(prev.index)-set(h.year);assert not removed;assert added==({2014} if 'full11' in sc else set())
 before=e[(e.scenario=='Previous_review')&(e.metric=='Change')].set_index('series_id');after=e[(e.scenario==sc)&(e.metric=='Change')].set_index('series_id');assert len(after)==55;assert np.allclose(before.estimate.sort_index(),after.estimate.sort_index());assert np.allclose(before.vi.sort_index(),after.vi.sort_index())
# Independent local point estimates and calendar HAC3; use original source data rows for the benchmark.
ll=pd.concat([local,a[(a.scenario=='Previous_review')&(a.MPA=='Point Vicente SMCA')&(a.Source=='PISCO')&(a.t>0)]],ignore_index=True);le=[]
for r in sl.itertuples():
 q=ll[(ll.scenario==r.scenario)&(ll.series_id==r.series_id)].sort_values('year');X=np.column_stack([np.ones(len(q)),q.t]);y=q.lnDiff.to_numpy();b=np.linalg.lstsq(X,y,rcond=None)[0];u=y-X@b;B=np.linalg.inv(X.T@X);yr=q.year.to_numpy();W=np.maximum(0,1-abs(yr[:,None]-yr[None,:])/4);scores=X*u[:,None];V=B@(scores.T@W@scores)@B*len(q)/(len(q)-2);err=max(abs(b[1]-r.estimate),abs(np.sqrt(V[1,1])-r.se));assert err<1e-8;le.append(err)
# Same anchored late-period observations must be identical between Apr-Oct and May-Oct.
for zn in ['INNER_MID','INNER_MID_OUTER']:
 x=local[local.scenario==f'AprOct_{zn}_late9_anchored'].set_index(['series_id','year']);y=local[local.scenario==f'MayOct_{zn}_late9_anchored'].set_index(['series_id','year']);assert np.allclose(x.lnDiff.sort_index(),y.lnDiff.sort_index())
# No DEEP zones or missing values; constant identical zone sets on both roles for all anchor/date-rule scenarios.
foot=pd.concat([pd.read_csv(out/'anchored_footprints.csv'),pd.read_csv(out/'date_matched_footprints.csv')]);roles=pd.concat([pd.read_csv(out/'anchored_role_values.csv'),pd.read_csv(out/'date_matched_role_values.csv')]);offset=pd.read_csv(out/'reconstruction.csv').set_index(['Taxa','Resp']).offset
for sc,g in foot.groupby('scenario'):
 assert not g.zone.isin(['DEEP','SPECIALO']).any();assert g.value.notna().all();sets=[set(q.zone) for _,q in g.groupby(['site','year'])];assert all(x==sets[0] for x in sets)
 assert len(set(tuple(sorted(q[['site','zone']].drop_duplicates().apply(tuple,axis=1))) for _,q in g.groupby('year')))==1
 for (year,tx,resp,status),q in g.groupby(['year','Taxa','Resp','status']):
  estimate=np.mean([np.mean(v.value) for _,v in q.groupby('site')]);saved=roles[(roles.scenario==sc)&(roles.year==year)&(roles.Taxa==tx)&(roles.Resp==resp)&(roles.status==status)].iloc[0];assert abs(estimate-saved.value)<1e-9
 for (tx,resp),q in roles[roles.scenario==sc].groupby(['Taxa','Resp']):
  v=q.pivot(index='year',columns='status',values='value');r=local[(local.scenario==sc)&(local.Taxa==tx)&(local.Resp==resp)].set_index('year');assert np.allclose(np.log((v.mpa+offset[(tx,resp)])/(v.reference+offset[(tx,resp)])),r.loc[v.index,'lnDiff'])
# Date rules recomputed from per-station gap ledger.
gaps=pd.read_csv(out/'station_pair_date_gaps.csv');rules=pd.read_csv(out/'date_rule_grid.csv')
for r in rules.itertuples():
 q=gaps[(gaps.season==r.season)&(gaps.zones==r.zones)].groupby('reference').max_gap_days.max();ok=set(q[q<=r.max_gap_rule_days].index);assert len(ok)==r.n_references;assert ok==(set(r.references.split(';')) if isinstance(r.references,str) else set())
h=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text())
for file,digest in h.items():assert hashlib.sha256(Path(file).read_bytes()).hexdigest()==digest
for file,digest in json.loads((out/'source_hashes.json').read_text()).items():assert hashlib.sha256(Path(file).read_bytes()).hexdigest()==digest
for directory in ['seasonal_zone_audit_20261101','pisco_source_reconciliation_20261101']:
 for file,digest in json.loads((Path('outputs')/directory/'manifest.json').read_text()).items():assert hashlib.sha256(Path(file).read_bytes()).hexdigest()==digest
r={'individual_contrasts_checked':len(errors),'max_effect_error':max(errors),'pooled_points_checked':len(errs),'max_pooled_error':max(errs),'local_slopes_checked':len(le),'max_local_error':max(le),'same169_post_histories':True,'all55_before_after_histories_unchanged':True,'only_four_Point_Vicente_responses_changed':True,'only_added_year_is_observed_2014':True,'same_shallow_zones_both_roles_constant_station_sets':True,'identical_late9_AprOct_MayOct_values':True,'source_and_previous_review_hashes_match':True,'protected_files_unchanged':len(h)}
(out/'verification.json').write_text(json.dumps(r,indent=2));print(r)
