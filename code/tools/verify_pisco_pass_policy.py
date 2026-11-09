from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
O=Path('outputs/pisco_pass_crane_20261101');e=pd.read_csv(O/'effects.csv');u=pd.read_csv(O/'model_rows.csv');d=pd.read_csv(O/'model_input.csv')
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
 if contrast=='Incomplete_passes':assert set(a.year)==set(b.year)
 else:assert set(a.year)<=set(b.year)
 if a.Source.iloc[0]!='PISCO':
  j=a.merge(b,on='year');assert (j.lnDiff_x==j.lnDiff_y).all()
protected=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text())
for p,h in protected.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
previous=json.loads(Path('outputs/pisco_fish_effort_20261101/manifest.json').read_text())
for p,h in previous.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
for p,h in json.loads((O/'input_hashes.json').read_text()).items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
single=pd.read_csv(O/'single_pass_comparison_links.csv');singleden=single[single.Resp=='Den'];assert len(singleden)==67 and singleden.MPA.nunique()==9
assert singleden.classification.eq('Confirmed omitted bottom pass').sum()==3
crane=pd.read_csv(O/'crane_comparison_links.csv');keys=['site','year','date','method','zone','transect'];cu=crane.drop_duplicates(keys);assert len(cu)==144 and set(cu.levels)<={'BOT;MID','BOT;CNMD'}
z=pd.read_csv(O/'proposed_paired_values.csv');assert len(z)==420 and z.complete_support.sum()==417
pool=pd.read_csv(O/'pool_checks.csv');assert pool.status.eq('fitted').all()
res=dict(effect_estimates_and_SEs_checked=len(e),max_estimate_error=max(v[0] for v in errors),max_se_error=max(v[1] for v in errors),single_pass_units_in_comparison_windows=67,confirmed_omitted_bottom=3,unresolved_single_pass_units=64,retained_crane_units=144,proposed_valid_history_year_rows=417,protected_files_unchanged=len(protected),previous_audit_artifacts_unchanged=len(previous),pools=len(pool),optimizer_fallbacks=int(pool.solver.ne('nlminb').sum()),limits='CRANE coverage contrast deliberately adds2004; incomplete-pass contrast uses identical years. Pooled fits not independently reimplemented.')
(O/'verification.json').write_text(json.dumps(res,indent=2)+'\n');print(res)
