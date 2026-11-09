from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy.stats import t
O=Path('outputs/cross_site_seasons_20261101')
f=pd.read_csv(O/'local_effects.csv');u=pd.read_csv(O/'model_rows.csv');effects=pd.read_csv(O/'effects.csv');p=pd.read_csv(O/'pooled.csv');lookup=f.set_index(['scenario','variant','series_id']);err=[]
assert not f.duplicated(['scenario','variant','series_id']).any()
for key,z in u.groupby(['scenario','variant','series_id']):
 r=lookup.loc[key];z=z.sort_values('year');X=np.column_stack([np.ones(len(z)),z.t]);y=z.lnDiff.to_numpy();b=np.linalg.lstsq(X,y,rcond=None)[0];res=y-X@b;B=np.linalg.inv(X.T@X);yr=z.year.to_numpy();W=np.maximum(0,1-abs(yr[:,None]-yr[None,:])/(min(3,int(yr.max()-yr.min()))+1));S=X*res[:,None];V=B@(S.T@W@S)@B*len(z)/(len(z)-2);e=max(abs(b[1]-r.estimate),abs(V[1,1]-r.vi));assert e<1e-7,(key,e);err.append(e);assert (z.t>0).all() and not z.year.duplicated().any();assert abs(r.low-(r.estimate-t.ppf(.975,r.df)*r.se))<1e-7
errs=[]
for r in p.itertuples():
 x=effects[(effects.scenario==r.scenario)&(effects.Taxa==r.Taxa)&(effects.Resp==r.Resp)];assert len(x)==r.k
 if np.isfinite(r.sigma_mpa):
  m=x.MPA.to_numpy();V=np.diag(x.vi)+r.sigma_mpa*(m[:,None]==m[None,:]);
  if np.isfinite(r.sigma_series):V+=np.eye(len(x))*r.sigma_series
  w=np.linalg.solve(V,np.ones(len(x)));er=abs(w@x.estimate/w.sum()-r.estimate);assert er<1e-7;errs.append(er)
 if np.isfinite(r.low):assert abs(r.low-(r.estimate-t.ppf(.975,r.df)*r.se))<1e-7
for _,z in p.groupby('scenario'):
 z=z[z.p.notna()].sort_values('p');qq=np.minimum(1,np.minimum.accumulate((z.p.to_numpy()*9/np.arange(1,len(z)+1))[::-1])[::-1]);assert np.allclose(qq,z.q)
c=pd.read_csv(O/'selected_coverage.csv').fillna({'years':'','lost_original':'','added_years':''});pairs=pd.read_csv(O/'pair_coverage.csv').fillna({'years':''});pg=pd.read_csv(O/'pair_year_gaps.csv');dates=pd.read_csv(O/'survey_dates_depths.csv');dates['date']=pd.to_datetime(dates[['year','month','day']]);windows={'MayOct':(5,10),'AprOct':(4,10),'AprNov':(4,11),'MarNov':(3,11),'AllMonths':(1,12)}
# Independently recompute physical date distances, not just trust a saved match flag.
ng=0
for (season,ps,rs,year,zn),z in pg[pg.max_gap_days.notna()].groupby(['season','protected','reference','year','zones']):
 lo,hi=windows[season];q=dates[(dates.year==year)&dates.month.between(lo,hi)];gg=[]
 for zone in (['INNER','MID'] if zn=='IM' else ['INNER','MID','OUTER']):
  pp=q[(q.site==ps)&(q.zone==zone)].date.to_numpy();rr=q[(q.site==rs)&(q.zone==zone)].date.to_numpy();assert len(pp) and len(rr);gg.extend(abs((pp[:,None]-rr[None,:])/np.timedelta64(1,'D')).ravel())
 assert np.allclose(z.max_gap_days,max(gg));ng+=1
for r in c.itertuples():
 q=pairs[(pairs.MPA==r.MPA)&(pairs.season==r.season)&(pairs.zones==r.zones)&(pairs.limit==r.limit)]
 if not r.scenario.startswith('Anchor_'):
  best=q.sort_values(['n','span','first','protected','reference'],ascending=[False,False,True,True,True]).iloc[0];assert best.protected==r.protected and best.reference==r.reference and best.years==r.years
 else:
  anchor=c[(c.MPA==r.MPA)&(c.scenario=='AprNov_IM_d90')].iloc[0];assert anchor.protected==r.protected and anchor.reference==r.reference
 gaps=pg[(pg.MPA==r.MPA)&(pg.season==r.season)&(pg.zones==r.zones)&(pg.protected==r.protected)&(pg.reference==r.reference)]
 good=gaps[gaps.protected_complete & gaps.reference_complete & ((gaps.max_gap_days<=r.limit)|(r.limit==0))];assert set(good.year)==set(int(x) for x in r.years.split(';') if x)
foot=pd.read_csv(O/'selected_footprints.csv');roles=pd.read_csv(O/'role_values.csv');annual=pd.read_csv(O/'local_annual.csv');offs=pd.read_csv(O/'reconstruction.csv').set_index(['MPA','Taxa','Resp']).offset
for (sc,m),z in foot.groupby(['scenario','MPA']):
 rule=c[(c.scenario==sc)&(c.MPA==m)].iloc[0];assert set(z.site)=={rule.protected,rule.reference};assert z.value.notna().all();assert not z.zone.isin(['DEEP','SPECIALO']).any()
 sets=z.groupby(['site','year']).zone.apply(lambda x:tuple(sorted(set(x))));assert len(set(sets))==1
 rr=roles[(roles.scenario==sc)&(roles.MPA==m)]
 rebuilt=z.groupby(['status','year','Taxa','Resp']).value.mean().sort_index();saved=rr.set_index(['status','year','Taxa','Resp']).value.sort_index();assert np.allclose(rebuilt,saved)
 for (tx,resp),q in rr.groupby(['Taxa','Resp']):
  v=q.pivot(index='year',columns='status',values='value');a=annual[(annual.scenario==sc)&(annual.MPA==m)&(annual.Taxa==tx)&(annual.Resp==resp)].set_index('year');off=offs[(m,tx,resp)];assert np.allclose(np.log((v.mpa+off)/(v.reference+off)),a.loc[v.index,'lnDiff'])
base=effects[effects.variant=='Previous_review'].set_index('series_id');changed=set(annual.series_id);assert len(changed)==20
for sc,z in effects.groupby('scenario'):
 other=z[~z.series_id.isin(changed)].set_index('series_id');b=base.loc[~base.index.isin(changed)];assert set(other.index)==set(b.index);assert np.allclose(other.estimate.sort_index(),b.estimate.sort_index());assert np.allclose(other.vi.sort_index(),b.vi.sort_index())
 if z.variant.iloc[0]=='Retention_only':assert np.allclose(z.set_index('series_id').estimate,base.loc[z.series_id].estimate)
old=pd.read_csv('outputs/pisco_source_reconciliation_20261101/effects.csv');old=old[(old.scenario=='Supported_corrections')&(old.metric=='Change')].set_index('series_id');ba=pd.read_csv(O/'unchanged_before_after.csv').set_index('series_id');assert len(ba)==55 and set(old.index)==set(ba.index);assert np.allclose(ba.estimate,old.loc[ba.index].estimate);assert np.allclose(ba.vi,old.loc[ba.index].vi);assert not set(ba.index)&changed
h=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text())
for file,digest in h.items():assert hashlib.sha256(Path(file).read_bytes()).hexdigest()==digest
for file,digest in json.loads((O/'source_hashes.json').read_text()).items():assert hashlib.sha256(Path(file).read_bytes()).hexdigest()==digest
for folder in ['point_vicente_seasons_20261101','pisco_source_reconciliation_20261101']:
 for file,digest in json.loads((Path('outputs')/folder/'manifest.json').read_text()).items():assert hashlib.sha256(Path(file).read_bytes()).hexdigest()==digest,(file,'changed')
res=dict(local_fits=len(err),max_local_error=max(err),pooled_points=len(errs),max_pooled_error=max(errs),physical_date_gap_checks=ng,coverage_rules=len(c),before_after_unchanged=55,protected_files_unchanged=len(h),source_hashes_match=True,previous_audits_unchanged=True,unchanged_non_target_effects=True)
(O/'verification.json').write_text(json.dumps(res,indent=2)+'\n');print(res)
