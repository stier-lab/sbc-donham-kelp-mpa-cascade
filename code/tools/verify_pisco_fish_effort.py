from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
O=Path('outputs/pisco_fish_effort_20261101');e=pd.read_csv(O/'effects.csv');u=pd.read_csv(O/'model_rows.csv');d=pd.read_csv(O/'model_input.csv')
r=pd.read_csv(O/'reconstruction_checks.csv');assert len(r)==34 and r.verified.all() and r.max_ratio_error.max()<1e-7
lookup={k:g for k,g in u.groupby(['series_id','scenario','metric'])};errors=[]
for v in e.itertuples():
 z=lookup[(v.series_id,v.scenario,v.metric)].sort_values('year');t=z.t.to_numpy();post=(t>0).astype(float);y=z.lnDiff.to_numpy();years=z.year.to_numpy(dtype=int)
 assert len(set(years))==len(years)
 X=np.column_stack([np.ones(len(t)),t]) if v.metric=='Post_slope' else np.column_stack([np.ones(len(t)),post,post*t])
 g=np.array([0,1]) if v.metric=='Post_slope' else np.array([0,1,t[post==1].mean()]);b=np.linalg.lstsq(X,y,rcond=None)[0];n,k=X.shape;res=y-X@b;score=X*res[:,None]
 padded=np.zeros((years.max()-years.min()+1,k));padded[years-years.min()]=score;meat=padded.T@padded;lag=min(3,len(padded)-1)
 for h in range(1,lag+1):
  cross=padded[h:].T@padded[:-h];meat+=(1-h/(lag+1))*(cross+cross.T)
 inv=np.linalg.inv(X.T@X);V=inv@meat@inv*n/(n-k);se=np.sqrt(g@V@g);est=g@b
 errors.append((abs(est-v.estimate),abs(se-v.se)));assert abs(est-v.estimate)<1e-8 and abs(se-v.se)<1e-8
# Confirm identical observations per before/after contrast; treatment values may differ.
for (c,sid),z in d.groupby(['contrast','series_id']):
 l=z[z.arm=='Before'];r=z[z.arm=='After'];assert set(l.year)==set(r.year)
 assert set(l.cluster)==set(r.cluster)
# Other programs unchanged in each paired contrast.
x=d[d.Source!='PISCO'];j=x[x.arm=='Before'].merge(x[x.arm=='After'],on=['contrast','series_id','year']);assert np.max(abs(j.lnDiff_x-j.lnDiff_y))==0
# Same event keys for complete-level full-water vs bottom-only sensitivity.
ev=pd.read_csv(O/'event_values.csv');keys=['site','year','date','method','zone','transect'];a=ev[ev.scenario=='Complete_noncan'];b=ev[ev.scenario=='Complete_BOT'];assert set(map(tuple,a[keys].to_numpy()))==set(map(tuple,b[keys].to_numpy()))
assert not ev.duplicated(keys+['scenario']).any()
# Averaging an un-revisited unit must reproduce its un-normalized value.
sv=pd.read_csv(O/'site_year_values.csv');f=pd.read_csv(O/'survey_units.csv');re=f.groupby(['site','year','zone','transect']).date.nunique().groupby(['site','year']).max();no=set(re[re==1].index)
l=sv[sv.scenario=='Legacy'].merge(sv[sv.scenario=='Visit_mean'],on=['site','year']);l=l[[k in no for k in zip(l.site,l.year)]]
for resp in ['Den','Bio']:
 z=l[[resp+'_x',resp+'_y']].dropna();assert np.max(abs(z.iloc[:,0]-z.iloc[:,1]))<1e-8
protected=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text())
for p,h in protected.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
sources=json.loads((O/'input_hashes.json').read_text())
for p,h in sources.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
pools=pd.read_csv(O/'pool_checks.csv');assert pools.status.eq('fitted').all()
res=dict(legacy_histories_reproduced=34,independently_checked_effects=len(e),max_estimate_error=max(v[0] for v in errors),max_se_error=max(v[1] for v in errors),matched_contrast_years=True,same_event_keys_bottom_sensitivity=True,other_programs_unchanged=True,protected_files_unchanged=len(protected),hashed_inputs=len(sources),pooled_models=len(pools),pooled_optimizer_fallbacks=int(pools.solver.ne('nlminb').sum()),limits='Pooled REML not independently reimplemented; constrained source support is not validation of habitat equivalence or missing-pass meaning.')
(O/'verification.json').write_text(json.dumps(res,indent=2)+'\n');print(res)
