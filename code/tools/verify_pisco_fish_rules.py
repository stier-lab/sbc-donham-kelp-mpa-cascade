from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
O=Path('outputs/pisco_fish_rules_20261101');e=pd.read_csv(O/'effects.csv');u=pd.read_csv(O/'model_rows.csv');d=pd.read_csv(O/'model_input.csv')
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


for (contrast,sid),z in d.groupby(['contrast','series_id']):
 a=z[z.arm=='Before'];b=z[z.arm=='After']
 if contrast=='Same_years':assert set(a.year)==set(b.year)
 else:assert set(b.year)<=set(a.year)
 if not (a.Source.iloc[0]=='PISCO' and a.Taxa.iloc[0]=='S. pulcher'):
  j=a.merge(b,on='year');assert len(a)==len(b) and (j.lnDiff_x==j.lnDiff_y).all()
old=pd.read_csv('outputs/full_followup_recovery_20261101/effects.csv');old=old[(old.scenario=='Full')&old.metric.isin(['Change','Post_slope'])]
oldcheck=e[e.scenario=='Total__Before'].merge(old,on=['series_id','metric'],suffixes=('_new','_old'))
assert len(oldcheck)==len(old) and (oldcheck.estimate_new-oldcheck.estimate_old).abs().max()<1e-8 and (oldcheck.se_new-oldcheck.se_old).abs().max()<1e-8
protected=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text())
for p,h in protected.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
for name in ['pisco_pass_crane_20261101','pisco_fish_effort_20261101']:
 for p,h in json.loads((Path('outputs')/name/'manifest.json').read_text()).items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
for p,h in json.loads((O/'input_hashes.json').read_text()).items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
approved=pd.read_csv(O/'approved_annual_input.csv');assert approved.series_id.nunique()==175 and len(approved)==3225
pools=pd.read_csv(O/'pool_checks.csv');assert pools.status.eq('fitted').all()
j=e[e.scenario=='Total__Before'].merge(e[e.scenario=='Total__After'],on=['series_id','metric'],suffixes=('_before','_after'));unchanged=j[~(j.Source_before.eq('PISCO')&j.Taxa_before.eq('S. pulcher'))];assert (unchanged.estimate_before-unchanged.estimate_after).abs().max()<1e-10
res=dict(individual_effects_independently_checked=len(e),max_estimate_error=max(v[0] for v in errors),max_se_error=max(v[1] for v in errors),previous_full_followup_effects_reproduced=len(oldcheck),approved_histories=175,approved_annual_rows=3225,protected_files_unchanged=len(protected),preceding_audits_unchanged=True,nonfish_effects_unchanged=True,pooled_fits=len(pools),optimizer_fallbacks=int(pools.solver.ne('nlminb').sum()),limits='Pooled REML not independently reimplemented; pathways/resilience not refitted; paper not updated.')
(O/'verification.json').write_text(json.dumps(res,indent=2)+'\n');print(res)
