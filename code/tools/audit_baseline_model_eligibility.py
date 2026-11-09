from pathlib import Path
import numpy as np,pandas as pd,json,hashlib
from scipy.stats import t,pearsonr,levene
root=Path.cwd();base=root/'outputs/researcher_scope_20261101';out=root/'outputs/model_eligibility_20261101';out.mkdir(exist_ok=True)
a=pd.read_csv(base/'annual_input.csv');effects=pd.read_csv(base/'effects.csv');eligible=set(effects.query("scenario=='Primary' and metric=='Change'").series_id)
def fit(y,year,X):
 n,p=X.shape;b=np.linalg.lstsq(X,y,rcond=None)[0];u=y-X@b;B=np.linalg.inv(X.T@X);h=np.sum((X@B)*X,axis=1);v=[B*(u@u)/(n-p),B@(X.T@np.diag((u/np.maximum(1-h,1e-8))**2)@X)@B]
 for lag in [1,3,5]:
  W=np.maximum(0,1-abs(year[:,None]-year[None,:])/(lag+1));S=X*u[:,None];v.append(B@(S.T@W@S)@B*n/(n-p))
 se=np.sqrt(np.maximum(0,np.array([np.diag(x) for x in v])));return b,u,se.max(axis=0),se
rows=[];pred=[]
for id,d in a.groupby('series_id',sort=False):
 d=d.sort_values('year');b=d[d.t<0];post=d[d.t.between(1,9)];r=dict(series_id=id,MPA=d.MPA.iloc[0],Taxa=d.Taxa.iloc[0],Resp=d.Resp.iloc[0],Source=d.Source.iloc[0],site_included=True,baseline_design_verified=bool(d.baseline_design_verified.all()),n_before=len(b),before_span=int(b.year.max()-b.year.min()) if len(b) else 0,n_after=len(post),prior_change_fit=id in eligible,before_years=';'.join(map(str,b.year)),p_linear=np.nan,p_curve=np.nan,p_shift=np.nan,slope=np.nan,slope_se=np.nan,drift9=np.nan,drift9_low90=np.nan,drift9_high90=np.nan,loo_drift9_min=np.nan,loo_drift9_max=np.nan,constant=False)
 for factor in [1.10,1.25,1.50,2.00]:r[f'compatible_{factor:.2f}']=False
 if len(b)>=3:
  yr=b.year.to_numpy(float);x=yr-yr.mean();y=b.lnDiff.to_numpy();r['constant']=np.ptp(y)<1e-10
  if not r['constant']:
   X=np.column_stack([np.ones(len(x)),x]);coef,u,se,allse=fit(y,yr,X);r.update(slope=coef[1],slope_se=se[1],p_linear=2*t.sf(abs(coef[1]/se[1]),len(y)-2),drift9=9*coef[1],drift9_low90=9*(coef[1]-t.ppf(.95,len(y)-2)*se[1]),drift9_high90=9*(coef[1]+t.ppf(.95,len(y)-2)*se[1]),iid_p_linear=2*t.sf(abs(coef[1]/allse[0,1]),len(y)-2),hac3_p_linear=2*t.sf(abs(coef[1]/allse[3,1]),len(y)-2))
   if len(b)>=4:
    loo=[9*np.linalg.lstsq(np.delete(X,i,axis=0),np.delete(y,i),rcond=None)[0][1] for i in range(len(y))];r['loo_drift9_min']=min(loo);r['loo_drift9_max']=max(loo)
   if len(b)>=6 and r['before_span']>=5:
    c,_,cs,_=fit(y,yr,np.column_stack([X,x*x]));r['p_curve']=2*t.sf(abs(c[2]/cs[2]),len(y)-3) if cs[2]>0 else np.nan
   if len(b)>=8:
    half=np.arange(len(y))>=len(y)//2;c,_,cs,_=fit(y,yr,np.column_stack([np.ones(len(y)),half]));r['p_shift']=2*t.sf(abs(c[1]/cs[1]),len(y)-2) if cs[1]>0 else np.nan
    r['variance_shift_p']=levene(u[~half],u[half],center='median').pvalue
   pairs=np.flatnonzero(np.diff(yr)==1)
   if len(pairs)>=4 and np.std(u[pairs])>0 and np.std(u[pairs+1])>0:r['residual_adjacent_r'],r['residual_adjacent_p']=pearsonr(u[pairs],u[pairs+1])
   for j in range(len(y)):pred.append(dict(series_id=id,year=yr[j],observed=y[j],fitted=(X@coef)[j],residual=u[j]))
 if not r['prior_change_fit']:status='Not eligible by coverage/design/estimability'
 elif r['constant']:status='Unresolved: constant baseline'
 elif r['p_linear']<.05:status='Detected baseline trend'
 elif r['p_curve']<.05 or r['p_shift']<.05:status='Nonlinear or level-shift flag'
 elif r['n_before']<6 or r['before_span']<5:status='Unresolved: sparse baseline'
 else:status='No detected departure; assess precision'
 r['status']=status
 for factor in [1.10,1.25,1.50,2.00]:
  margin=np.log(factor);r[f'compatible_{factor:.2f}']=status=='No detected departure; assess precision' and r['drift9_low90']>-margin and r['drift9_high90']<margin and r['loo_drift9_min']>-margin and r['loo_drift9_max']<margin
 r['permissive_no_linear_flag']=r['prior_change_fit'] and not r['constant'] and np.isfinite(r['p_linear']) and r['p_linear']>=.05
 # Illustrative relaxed sample-size gate, still requiring precision/LOO; absence of shape test is flagged.
 r['relaxed_n4_50pct']=r['prior_change_fit'] and not r['constant'] and len(b)>=4 and r['before_span']>=3 and r['p_linear']>=.05 and not(r['p_curve']<.05 or r['p_shift']<.05) and r['drift9_low90']>-np.log(1.5) and r['drift9_high90']<np.log(1.5) and r['loo_drift9_min']>-np.log(1.5) and r['loo_drift9_max']<np.log(1.5)
 rows.append(r)
ledger=pd.DataFrame(rows);ledger.to_csv(out/'model_eligibility.csv',index=False);pd.DataFrame(pred).to_csv(out/'baseline_observations.csv',index=False)
sets={'All_before_after':ledger.prior_change_fit,'No_detected_linear_trend':ledger.permissive_no_linear_flag,'Adequate_no_detected_departure':ledger.status.eq('No detected departure; assess precision'),'Tolerance_10':ledger['compatible_1.10'],'Tolerance_25':ledger['compatible_1.25'],'Tolerance_50':ledger['compatible_1.50'],'Tolerance_100':ledger['compatible_2.00'],'Relaxed_n4_tolerance50':ledger.relaxed_n4_50pct};parts=[];membership=[]
for scenario,mask in sets.items():
 ids=set(ledger.loc[mask,'series_id']);e=effects[(effects.scenario=='Primary')&((effects.metric=='Post_difference')|effects.series_id.isin(ids))].copy();e['scenario']=scenario;parts.append(e)
 for id in sorted(ids):membership.append(dict(scenario=scenario,series_id=id))
pd.concat(parts).to_csv(out/'effects.csv',index=False);pd.DataFrame(membership).to_csv(out/'model_membership.csv',index=False)
files=[base/'annual_input.csv',base/'effects.csv',out/'PLAN.md',out/'model_eligibility.csv',out/'effects.csv'];(out/'freeze_manifest.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},indent=2))
print(ledger[ledger.prior_change_fit].status.value_counts().to_string());print('Compatible counts:',{s:int(m.sum()) for s,m in sets.items()})

orig_rows=pd.read_csv(base/'model_rows.csv');orig_rows=orig_rows[orig_rows.scenario=='Primary'];modelrows=[]
for scenario,mask in sets.items():
 ids=set(ledger.loc[mask,'series_id']);z=orig_rows[(orig_rows.metric=='Post_difference')|orig_rows.series_id.isin(ids)].copy();z['scenario']=scenario;modelrows.append(z)
pd.concat(modelrows).to_csv(out/'model_rows.csv',index=False)
(out/'protected_hashes.json').write_text((base/'protected_hashes.json').read_text())
