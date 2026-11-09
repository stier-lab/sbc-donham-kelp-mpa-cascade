"""Joint comparison fits with shared physical-station and temporal process terms.
Known measurement covariance comes from the same nonlinear draws across uses.
Separate before-after changes and post-protection trends are never pooled.
"""
from common import *
import sys
sys.path.insert(0,str(ROOT/'code/tools'))
from linear_mixed import fit
from scipy.linalg import block_diag
from scipy.stats import t as student,shapiro,spearmanr
from scipy.optimize import minimize_scalar
overlap=read_tracked(ROOT/'data/external/source_metadata_2026-10-16/satellite_reference_overlap.csv');crosswalk=read_tracked(ROOT/'data/external/source_metadata_2026-10-16/satellite_site_crosswalk.csv').set_index('Site_#').MPA.to_dict();matches=read_tracked(ROOT/'data/external/tom_bell_2026-10-20/current_input_matches.csv');matches=matches[matches.exact];geom=[]
for rr in overlap.itertuples():
 aa=matches[(matches.site==rr.site_a)&(matches.group==rr.group_a)];bb=matches[(matches.site==rr.site_b)&(matches.group==rr.group_b)]
 if len(aa)==1 and len(bb)==1:
  ka='Landsat|'+str(crosswalk[rr.site_a])+'|reference|'+str(int(aa.rep.iloc[0]));kb='Landsat|'+str(crosswalk[rr.site_b])+'|reference|'+str(int(bb.rep.iloc[0]));geom.append(dict(a=ka,b=kb,correlation=float(np.sqrt(rr.fraction_a*rr.fraction_b)),shared_pixels=rr.n_shared_pixels,basis='verified mask-to-delivery crosswalk; iid-pixel overlap working kernel'))
pd.DataFrame(geom).to_csv(OUT/'satellite_shared_reference_kernel.csv',index=False)
an=read_tracked(OUT/'annual_values.csv');an=an[an.response_supported];draws=np.load(PRIVATE/'annual_measurement_draws.npz')['draws'];st=read_tracked(OUT/'station_annual_values.csv');elig=read_tracked(OUT/'model_eligibility.csv').set_index('series_id');out=[];diags=[];indexed=[];covblocks=[]
for (source,taxon,response),allg in an.groupby(['Source','Taxa','Resp']):
 for model in ['post_trend','flat_before_after','pretrend_sensitivity']:
  parts=[];designs=[];terms=[];refs=[];scales=[]
  for series,g in allg.groupby('series_id'):
   if series not in elig.index or elig.loc[series,'after_years']<3:continue
   if model!='post_trend' and elig.loc[series,'before_years']<POLICY['minimum_before_years']:continue
   g=g[g.t.ge(0)] if model=='post_trend' else g;g=g.sort_values('year');tm=g.t.to_numpy(float);n=len(g);post=(tm>=0).astype(float);rec=np.maximum(tm,0)
   names=['intercept','recovery_slope'] if model=='post_trend' else ['baseline','level_change','recovery_slope'] if model=='flat_before_after' else ['baseline','pre_slope','level_change','slope_change']
   X=np.column_stack([np.ones(n),tm]) if model=='post_trend' else np.column_stack([np.ones(n),post,rec]) if model=='flat_before_after' else np.column_stack([np.ones(n),tm,post,rec])
   if np.linalg.matrix_rank(X)<len(names) or n<=len(names):continue
   scale=max(float(g.ln_ratio.std()),.1);parts.append(g);designs.append(X);scales.extend([scale]*n)
   terms.extend([dict(series_id=series,Source=source,Taxa=taxon,Resp=response,MPA=g.MPA.iloc[0],model=model,term=name,n_years=n,before_years=int(g.t.lt(0).sum()),after_years=int(g.t.ge(0).sum()),mean_post_time=float(g.loc[g.t.ge(0),'t'].mean())) for name in names])
   for r in g.itertuples():
    z=st[(st.series_id==series)&st.year.eq(r.year)];rweights={}
    for role,sign in [('mpa',1.),('reference',-1.)]:
     sites=sorted(z.loc[z.status.eq(role),'site'].unique())
     for site in sites:rweights[source+'|'+site]=sign/len(sites)
    refs.append(rweights)
  if not parts:continue
  h=pd.concat(parts,ignore_index=True);X=block_diag(*designs);y=h.ln_ratio.to_numpy();sc=np.array(scales);N=len(h);ncoef=X.shape[1];d=draws[h.annual_index.to_numpy(int)];S=measurement_covariance(d,h.year.to_numpy(),source,response,taxon)/np.outer(sc,sc);Xn=X/sc[:,None];yn=y/sc
  physical=sorted(set().union(*[set(r) for r in refs]));Z=np.array([[r.get(site,0) for site in physical] for r in refs]);V=np.eye(len(physical))
  for rr in geom:
   if rr['a'] in physical and rr['b'] in physical:
    ia=physical.index(rr['a']);ib=physical.index(rr['b']);V[ia,ib]=V[ib,ia]=rr['correlation']
  if np.linalg.eigvalsh(V).min()<-1e-8:raise ValueError('Invalid shared polygon covariance kernel')
  shared=Z@V@Z.T;shared=shared/np.sqrt(np.outer(np.diag(shared),np.diag(shared)));ids=h.series_id.to_numpy();time=h.year.to_numpy(int)
  projection=np.eye(N)-Xn@np.linalg.pinv(Xn)
  cache={}
  def at_rho(rho):
   rho=float(rho)
   if rho not in cache:
    ar=rho**np.abs(time[:,None]-time[None,:]);K_ind=(ids[:,None]==ids[None,:])*ar/np.outer(sc,sc);K_shared=shared*ar/np.outer(sc,sc)
    projected=np.column_stack([(projection@k@projection).ravel() for k in [K_ind,K_shared]])
    norms=np.linalg.norm(projected,axis=0);shared_identifiable=bool(np.all(norms>1e-10) and np.linalg.matrix_rank(projected/norms,tol=1e-7)==2)
    cache[rho]=fit(yn,Xn,S,[K_ind,K_shared] if shared_identifiable else [K_ind]);cache[rho]['shared_identifiable']=shared_identifiable
   return cache[rho]
  grid=[-.8,-.5,0.,.3,.6,.8,.95]
  objectives=[at_rho(rho)['obj'] for rho in grid];candidates=list(grid)
  intervals=[]
  if max(objectives)-min(objectives)>1e-8*(1+abs(min(objectives))):
   if objectives[0]<=objectives[1]:intervals.append((grid[0],grid[1]))
   if objectives[-1]<=objectives[-2]:intervals.append((grid[-2],grid[-1]))
   for k in range(1,len(grid)-1):
    if objectives[k]<=objectives[k-1] and objectives[k]<=objectives[k+1]:intervals.append((grid[k-1],grid[k+1]))
  for bounds in intervals:
   opt=minimize_scalar(lambda rho:at_rho(rho)['obj'],bounds=bounds,method='bounded',options={'xatol':.005})
   if not opt.success:raise RuntimeError('AR1 profile optimization failed')
   candidates.append(float(opt.x))
  rhohat=min(candidates,key=lambda rho:at_rho(rho)['obj'])
  if np.sum(at_rho(rhohat)['q'])<1e-8:rhohat=0. # No process variance: rho cannot be identified.
  contrast_rows=[];contrast_terms=[]
  for series in dict.fromkeys(t['series_id'] for t in terms):
   local=[(j,t) for j,t in enumerate(terms) if t['series_id']==series]
   if model!='post_trend':
    vec=np.zeros(ncoef)
    for j,t in local:
     if t['term']=='level_change':vec[j]=1.
     if t['term'] in ['recovery_slope','slope_change']:vec[j]=t['mean_post_time']
    meta=dict(local[0][1]);meta['term']='average_observed_post_change';contrast_rows.append(vec);contrast_terms.append(meta)
  T=np.vstack([np.eye(ncoef),*contrast_rows]) if contrast_rows else np.eye(ncoef);terms_final=terms+contrast_terms
  for label,rho in [('estimated_AR1',rhohat)]+[('fixed_AR1',v) for v in [0.,.3,.6]]:
   f=at_rho(rho);shared_identifiable=f['shared_identifiable'];df=N-ncoef;crit=student.ppf(.975,df)
   beta=T@f['beta'];bc=T@f['bc']@T.T
   for j,t in enumerate(terms_final):
    se=np.sqrt(max(bc[j,j],0));out.append(dict(**t,analysis=label,rho=rho,estimate=beta[j],se=se,low=beta[j]-crit*se,high=beta[j]+crit*se,df=df))
   residual=yn-Xn@f['beta'];white=np.linalg.solve(np.linalg.cholesky(f['C']),residual);normal_p=float(shapiro(white).pvalue) if len(white)>=3 else np.nan;variance_pattern=float(spearmanr(np.abs(white),Xn@f['beta']).statistic)
   diags.append(dict(Source=source,Taxa=taxon,Resp=response,model=model,analysis=label,rho=rho,ar_boundary=rho in [-.8,.95],ar_identifiable=bool(np.sum(f['q'])>=1e-8),comparisons=len(parts),years=N,parameters=ncoef,df=df,process_variance=f['q'][0],shared_process_variance=f['q'][1] if shared_identifiable else np.nan,shared_component_identifiable=shared_identifiable,boundary=bool(np.any(f['q']<1e-8)),upper_boundary=f['upper_boundary'],optimizer_attempts=f['optimizer_attempts'],whitened_normality_p=normal_p,absolute_residual_fitted_correlation=variance_pattern,max_abs_whitened_residual=float(np.max(np.abs(white))),condition_X=np.linalg.cond(Xn),covariance_min_eigenvalue=float(np.linalg.eigvalsh(f['bc']).min())))
   if label=='estimated_AR1':indexed.extend(terms_final);covblocks.append(bc)
pd.DataFrame(out).to_csv(OUT/'joint_model_results.csv',index=False);pd.DataFrame(diags).to_csv(OUT/'joint_model_diagnostics.csv',index=False);pd.DataFrame(indexed).to_csv(OUT/'joint_coefficient_index.csv',index=False);C=block_diag(*covblocks);np.savez_compressed(PRIVATE/'joint_coefficient_covariance.npz',covariance=C)
# Conditional source-specific summaries and across-program meta-analysis. The latter
# adds MPA-level heterogeneity; it is distinct from within-comparison recovery.
r=pd.DataFrame(out);r=r[r.analysis.eq('estimated_AR1')].reset_index(drop=True);idx=pd.DataFrame(indexed);assert len(r)==len(idx)
pooled=[]
for (taxon,response,model,term),g in r[r.term.isin(['recovery_slope','slope_change','level_change','average_observed_post_change'])].groupby(['Taxa','Resp','model','term']):
 for source,z in [('All programs',g)]+list(g.groupby('Source')):
  j=z.index.to_numpy();known=C[np.ix_(j,j)];values=z.estimate.to_numpy();n=len(z)
  if n<3:continue
  # One MPA random effect connects program estimates from the same reserve.
  K=(z.MPA.to_numpy()[:,None]==z.MPA.to_numpy()[None,:]).astype(float)
  Xmeta=pd.get_dummies(z.Source,dtype=float).to_numpy() if source=='All programs' else np.ones((n,1))
  if n<=Xmeta.shape[1]:continue
  f=fit(values,Xmeta,known,[K]);contrast=np.ones(Xmeta.shape[1])/Xmeta.shape[1];estimate=float(contrast@f['beta']);se=float(np.sqrt(contrast@f['bc']@contrast));clusters=z.MPA.nunique();crit=student.ppf(.975,max(1,clusters-1))
  pooled.append(dict(Taxa=taxon,Resp=response,model=model,term=term,Source=source,estimate=estimate,se=se,low=estimate-crit*se,high=estimate+crit*se,comparisons=n,MPAs=clusters,between_MPA_variance=f['q'][0],boundary=f['q'][0]<1e-8,upper_boundary=f['upper_boundary'],program_weighting='equal program means; source fixed effects' if source=='All programs' else 'within program'))
pd.DataFrame(pooled).to_csv(OUT/'pooled_results.csv',index=False);write_inputs();print('Joint coefficients',len(out),'pooled summaries',len(pooled))
