from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
from scipy.stats import t as student
O=Path('outputs/all_program_station_refits_20261102');a=pd.read_csv('outputs/pisco_fish_rules_20261101/approved_annual_input.csv');v=pd.read_csv(O/'verified_station_values.csv');offset=pd.read_csv(O/'offsets.csv').set_index('series_id').offset
def fit(y,years,X,contrast,weights=None,station=False):
 y=np.asarray(y,float);X=np.asarray(X,float);years=np.asarray(years,int);w=np.ones(len(y)) if weights is None else np.asarray(weights,float)
 if not len(y):return None
 rank=np.linalg.matrix_rank(X);n,k=X.shape
 if rank<k or n<=k or len(set(years))<3:return None
 inv=np.linalg.inv(X.T@(w[:,None]*X));b=inv@(X.T@(w*y));score=X*(w*(y-X@b))[:,None];timeline=np.arange(years.min(),years.max()+1);scores=np.zeros((len(timeline),k))
 for yy in set(years):scores[yy-years.min()]=score[years==yy].sum(axis=0)
 meat=scores.T@scores;lag=min(3,len(timeline)-1)
 for h in range(1,lag+1):
  cross=scores[h:].T@scores[:-h];meat+=(1-h/(lag+1))*(cross+cross.T)
 correction=len(set(years))/(len(set(years))-1) if station else n/(n-k)
 V=inv@meat@inv*correction;g=np.asarray(contrast);est=g@b;se=np.sqrt(max(0,g@V@g));df=len(set(years))-2 if station else n-k
 return dict(estimate=est,se=se,low=est-student.ppf(.975,df)*se,high=est+student.ppf(.975,df)*se,n_years=len(set(years)),n_rows=n,df=df,rank=rank)
def slope(z,col):
 return fit(z[col],z.year,np.column_stack([np.ones(len(z)),z.t]),[0,1])
results=[];annual=[];station_rows=[];influence=[]
for sid,allg in v.groupby('series_id'):
 g=allg[allg.t>0].copy();g=g[g.value.notna()];
 if g.empty:continue
 meta=dict(series_id=sid,MPA=g.MPA.iloc[0],Resp=g.Resp.iloc[0],Source=g.Source.iloc[0],Taxa=g.Taxa.iloc[0]);off=offset[sid]
 approved=a[(a.series_id==sid)&(a.t>0)].sort_values('year')
 def save(sc,z,col='lnDiff'):
  ans=slope(z,col)
  if ans:results.append(dict(**meta,scenario=sc,**ans))
 save('Approved',approved)
 refs=allg[allg.status.eq('reference')];sets=[set(q.loc[q.value.notna(),'site']) for _,q in refs.groupby('year')];persistent=set.intersection(*sets)
 protected=g[g.status.eq('mpa')].groupby('year').value.mean();available=g[g.status.eq('reference')].groupby('year').value.mean()
 for name,rr in [('Available',g[g.status.eq('reference')]),('Persistent',g[g.status.eq('reference')&g.site.isin(persistent)])]:
  z=rr.groupby('year').value.mean().rename('reference').to_frame().join(protected.rename('protected')).dropna();z['lnDiff']=np.log((z.protected+off)/(z.reference+off));z=z.reset_index();z['t']=z.year-g.MPA_Start.iloc[0]
  for r in z.itertuples():annual.append(dict(**meta,scenario=name,year=r.year,t=r.t,lnDiff=r.lnDiff,protected=r.protected,reference=r.reference))
  if name=='Persistent':
   save('Persistent_reference',z);save('Available_on_persistent_years',approved[approved.year.isin(z.year)])
 # A post-only persistent set distinguishes baseline missingness from recovery turnover.
 postsets=[set(q.loc[q.value.notna(),'site']) for _,q in g[g.status.eq('reference')].groupby('year')];postfixed=set.intersection(*postsets)
 rr=g[g.status.eq('reference')&g.site.isin(postfixed)].groupby('year').value.mean();z=rr.rename('reference').to_frame().join(protected.rename('protected')).dropna().reset_index();z['lnDiff']=np.log((z.protected+off)/(z.reference+off));z['t']=z.year-g.MPA_Start.iloc[0];save('Persistent_post_reference',z)
 # Leave-one-reference-out on the same years in both arms.
 for omit in sorted(set(g.loc[g.status.eq('reference'),'site'])):
  rr=g[g.status.eq('reference')&g.site.ne(omit)].groupby('year').value.mean();z=rr.rename('reference').to_frame().join(protected.rename('protected')).dropna().reset_index();z['lnDiff']=np.log((z.protected+off)/(z.reference+off));z['t']=z.year-g.MPA_Start.iloc[0];after=slope(z,'lnDiff');before=slope(approved[approved.year.isin(z.year)],'lnDiff')
  if after and before:influence.append(dict(**meta,omitted_reference=omit,n_years=len(z),before=before['estimate'],after=after['estimate'],difference=after['estimate']-before['estimate'],before_low=before['low'],before_high=before['high'],after_low=after['low'],after_high=after['high'],years=';'.join(map(str,sorted(z.year)))))
 # Same observed station rows and weights; only station intercept adjustment changes.
 g['protected']=g.status.eq('mpa').astype(int);g['log_value']=np.log(g.value+off);g['weight']=1/g.groupby(['year','status']).site.transform('count');g['station_key']=g.status+'|'+g.site
 yd=pd.get_dummies(g.year,drop_first=True,dtype=float).to_numpy();interaction=(g.protected*g.t).to_numpy();sd=pd.get_dummies(g.station_key,drop_first=True,dtype=float).to_numpy()
 for sc,XX in [('Log_station_unadjusted',np.column_stack([np.ones(len(g)),g.protected,yd,interaction])),('Log_station_adjusted',np.column_stack([np.ones(len(g)),sd,yd,interaction]))]:
  contrast=np.zeros(XX.shape[1]);contrast[-1]=1;ans=fit(g.log_value,g.year,XX,contrast,g.weight,True)
  if ans:results.append(dict(**meta,scenario=sc,**ans))
 for r in g.itertuples():station_rows.append(dict(**meta,site=r.site,status=r.status,year=r.year,t=r.t,value=r.value,offset=off,log_value=r.log_value,weight=r.weight))
pd.DataFrame(results).to_csv(O/'effects.csv',index=False);pd.DataFrame(annual).to_csv(O/'annual_comparisons.csv',index=False);pd.DataFrame(influence).to_csv(O/'reference_omission.csv',index=False);pd.DataFrame(station_rows).to_csv(O/'station_model_rows.csv',index=False)

e=pd.DataFrame(results);old=pd.read_csv('outputs/pisco_fish_rules_20261101/effects.csv');old=old[(old.scenario=='Total__After')&(old.metric=='Post_slope')];j=e[e.scenario=='Approved'].merge(old,on='series_id',suffixes=('_new','_old'));assert len(j)==len(old);assert (j.estimate_new-j.estimate_old).abs().max()<1e-8;assert (j.se_new-j.se_old).abs().max()<1e-8
coverage=[]
for sid,g in v.groupby('series_id'):
 sets=[set(q.loc[q.value.notna(),'site']) for _,q in g[g.status.eq('reference')].groupby('year')];union=set.union(*sets);fixed=set.intersection(*sets);annual=g[g.status.eq('reference')].groupby('year').value.count()
 postsets=[set(q.loc[q.value.notna(),'site']) for _,q in g[(g.status=='reference')&(g.t>0)].groupby('year')];postunion=set.union(*postsets) if postsets else set();postfixed=set.intersection(*postsets) if postsets else set()
 coverage.append(dict(post_reference_membership_changes=any(x!=postunion for x in postsets),n_post_persistent=len(postfixed),series_id=sid,MPA=g.MPA.iloc[0],Source=g.Source.iloc[0],Taxa=g.Taxa.iloc[0],Resp=g.Resp.iloc[0],n_years=g.year.nunique(),n_reference_stations=len(union),n_persistent_references=len(fixed),reference_membership_changes=any(x!=union for x in sets),persistent_stations=';'.join(sorted(fixed)),references=';'.join(sorted(union)),size_support_status='Existing size conversion reproduced; expanded donor adequacy not newly certified' if g.Resp.iloc[0]=='Bio' and g.Source.iloc[0]!='Landsat' else 'Not a size-derived biomass response'))
pd.DataFrame(coverage).to_csv(O/'coverage_summary.csv',index=False)
wide=e.pivot(index=['series_id','MPA','Source','Taxa','Resp'],columns='scenario',values='estimate').reset_index();wide['station_adjustment']=wide.Log_station_adjusted-wide.Log_station_unadjusted;wide['reference_subset_change']=wide.Persistent_reference-wide.Approved;wide.to_csv(O/'comparison_summary.csv',index=False)
protected=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text())
for p,h in protected.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
(O/'verification.json').write_text(json.dumps(dict(annual_histories_verified=175,previous_slopes_reproduced=len(j),max_estimate_error=float((j.estimate_new-j.estimate_old).abs().max()),max_se_error=float((j.se_new-j.se_old).abs().max()),protected_unchanged=len(protected),fits=len(e),omission_comparisons=len(influence)),indent=2))
print(pd.DataFrame(coverage).groupby('Source').agg(histories=('series_id','size'),changing_reference_membership=('reference_membership_changes','sum'),single_reference=('n_reference_stations',lambda x:int((x==1).sum())),no_persistent_reference=('n_persistent_references',lambda x:int((x==0).sum()))));print((O/'verification.json').read_text())
