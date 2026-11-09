from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
from scipy.stats import t
p=Path('outputs/change_attribution_20261101');e=pd.read_csv(p/'effects.csv');o=pd.read_csv(p/'observations.csv');z=pd.read_csv(p/'pooled.csv');old=pd.read_csv('data/sumstats_final.csv');old['series_id']=old[['MPA','Taxa','Source','Resp']].agg('|'.join,axis=1);orig=old.set_index('series_id')
assert not e.duplicated(['stage','metric','series_id']).any()
assert not o.duplicated(['stage','metric','series_id','year']).any()
errors=[];n=0
for (stage,metric,id),d in o.groupby(['stage','metric','series_id']):
 target=e[(e.stage==stage)&(e.metric==metric)&(e.series_id==id)]
 if stage=='membership_fixed':target=e[(e.stage=='M2_additions_only')&(e.metric==metric)&(e.series_id==id)]
 if target.empty:continue
 r=target.iloc[0];d=d.sort_values('year');y=d.lnDiff.to_numpy();yr=d.year.to_numpy();post=d.post.to_numpy();tm=d.t.to_numpy();nb=sum(post==0)
 X=np.column_stack([np.ones(len(d)),tm]) if nb==0 else np.column_stack([np.ones(len(d)),post,post*tm])
 b=np.linalg.lstsq(X,y,rcond=None)[0];u=y-X@b;B=np.linalg.inv(X.T@X);nn,pp=X.shape
 g=X.mean(axis=0) if nb==0 else np.array([0 if metric=='Change' else 1,1,tm[post==1].mean()])
 lag=min(3,int(yr.max()-yr.min()));W=np.maximum(0,1-abs(yr[:,None]-yr[None,:])/(lag+1));S=X*u[:,None];V=B@(S.T@W@S)@B*nn/(nn-pp)
 er=max(abs(g@b-r.estimate),abs(g@V@g-r.vi));assert er<1e-7,(stage,id,er);errors.append(er);n+=1
 if stage=='B1_common_definition_same_observations':assert len(d)==orig.loc[id,'N']
# Each B model step holds history membership fixed. Designation/year-window steps explicitly alter years.
for metric in ['Change','Post_difference']:
 sets=[set(e[(e.stage==s)&(e.metric==metric)].series_id) for s in ['B0_matched_legacy','B1_common_definition_same_observations','B2_calendar_designation_rule','B3_nine_year_window','B4_target_observations']]
 assert all(x==sets[0] for x in sets)
# Pooled GLS points independently reconstructed from reported variance components.
pe=[]
for _,r in z.iterrows():
 x=e[(e.stage==r.stage)&(e.metric==r.metric)&(e.Taxa==r.Taxa)&(e.Resp==r.Resp)]
 if pd.isna(r.sigma):continue
 sig=[float(v) for v in str(r.sigma).split(';')];mp=x.MPA.to_numpy();V=np.diag(x.vi.to_numpy())+sig[0]*(mp[:,None]==mp[None,:]);
 if len(sig)>1:V+=np.eye(len(x))*sig[1]
 one=np.ones(len(x));w=np.linalg.solve(V,one);est=w@x.estimate.to_numpy()/(w@one);er=abs(est-r.estimate);assert er<1e-7;pe.append(er)
 if np.isfinite(r.low):assert abs(r.low-(r.estimate-t.ppf(.975,r.df)*r.se))<1e-7
# Inference-only contrasts preserve estimates and effect inputs exactly.
for a,b in [('A1_pool_structure','A2_cluster_uncertainty'),('U1_same_data_IID_model','U2_same_data_IID_CR2'),('U3_same_data_HAC_model','U4_same_data_HAC_CR2')]:
 x=z[z.stage==a].set_index(['metric','Taxa','Resp']);y=z[z.stage==b].set_index(['metric','Taxa','Resp']);assert np.allclose(x.estimate,y.estimate,atol=1e-10)
for metric in ['Change','Post_difference']:
 base=e[(e.stage=='M2_additions_only')&(e.metric==metric)].set_index('series_id')
 for s in ['M0_old_roster','M1_removals_only','M3_final_roster']:
  x=e[(e.stage==s)&(e.metric==metric)].set_index('series_id');assert np.allclose(x[['estimate','vi']],base.loc[x.index,['estimate','vi']])
for _,d in z.groupby(['stage','metric']):
 valid=d.p.notna();ps=d.loc[valid,'p'].to_numpy();order=np.argsort(ps);q=np.minimum(1,np.minimum.accumulate((ps[order]*9/np.arange(1,len(ps)+1))[::-1])[::-1]);assert np.allclose(q,d.loc[valid,'q'].to_numpy()[order])
h=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text());print('hash schema',type(h).__name__)
if isinstance(h,dict):
 for name,digest in h.items():assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==digest,name
report=dict(individual_contrasts=n,max_individual_error=max(errors),pooled_points=len(pe),max_pooled_error=max(pe),protected_files=len(h),matched_membership=True,unchanged_membership_inputs=True,bh_verified=True)
(p/'verification.json').write_text(json.dumps(report,indent=2));print(report)
