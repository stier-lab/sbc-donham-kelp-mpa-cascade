"""Conditional joint size/count observation model and temporal-process diagnostic.
Not a primary model or a claim of total uncertainty. No new supported cells.
"""
from pathlib import Path
import argparse, json, hashlib
from datetime import datetime, timezone
from scipy.stats import skew
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize, minimize_scalar
from scipy.stats import norm

def stable_id(parts):
    return hashlib.sha256('|'.join(map(str,parts)).encode()).hexdigest()[:20]

def design(d,center):
    return np.column_stack([np.ones(len(d)),d.zone.eq('MID').astype(float),(d.year.to_numpy()-center)/10])

def kernels(a,b):
    si=a.site.to_numpy()[:,None]==b.site.to_numpy()[None,:]
    sy=si & (a.year.to_numpy()[:,None]==b.year.to_numpy()[None,:])
    event=sy & (a.zone.to_numpy()[:,None]==b.zone.to_numpy()[None,:]) & (a.date.to_numpy()[:,None]==b.date.to_numpy()[None,:])
    return [si.astype(float),sy.astype(float),event.astype(float)]

def fit_reml(y,X,known,K,start=None):
    n,p=X.shape
    if np.linalg.matrix_rank(X)!=p or n<=p: raise ValueError('Unidentified size-model design or no residual degrees of freedom')
    def calc(theta,ret=False):
        cov=np.diag(known)+sum(t*k for t,k in zip(theta,K))+np.eye(n)*1e-9
        cf=cho_factor(cov,lower=True,check_finite=False);vx=cho_solve(cf,X,check_finite=False);vy=cho_solve(cf,y,check_finite=False)
        info=X.T@vx;bc=np.linalg.inv(info);beta=bc@(X.T@vy);res=y-X@beta
        loss=2*np.log(np.diag(cf[0])).sum()+np.linalg.slogdet(info)[1]+res@cho_solve(cf,res,check_finite=False)
        return (float(loss),beta,bc,cov,cf) if ret else float(loss)
    scale=max(float(np.var(y)),.01);upper=max(10.,10*scale)
    x0=np.full(len(K),scale/3) if start is None else np.clip(start,0,upper)
    bounds=[(0,upper) if np.any(k) else (0,0) for k in K]
    opt=minimize(calc,x0,method='L-BFGS-B',bounds=bounds,options={'maxiter':180,'ftol':1e-9})
    loss,beta,bc,cov,cf=calc(opt.x,True)
    return dict(theta=opt.x,beta=beta,beta_cov=bc,cov=cov,cf=cf,success=bool(opt.success),optimizer_message=str(opt.message),loss=loss,boundary=int((opt.x<1e-7).sum()))

def prediction(f,y,X,xt,Kcross,Ktt):
    if len(xt)==0: return np.empty(0),np.empty((0,0))
    cross=sum(t*k for t,k in zip(f['theta'],Kcross));gt=sum(t*k for t,k in zip(f['theta'],Ktt))
    cv=cho_solve(f['cf'],cross.T,check_finite=False).T
    delta=xt-cv@X;mu=xt@f['beta']+cv@(y-X@f['beta'])
    cov=gt-cv@cross.T+delta@f['beta_cov']@delta.T
    cov=(cov+cov.T)/2
    w,Q=np.linalg.eigh(cov)
    if w.min() < -1e-6: raise ValueError('non-PSD latent covariance')
    return mu,Q@np.diag(np.sqrt(np.maximum(w,0)))

def fit_trend(y,years,S,rho):
    if len(y)<4 or len(set(years))<4: raise ValueError('At least four observed years required for this process diagnostic')
    X=np.column_stack([np.ones(len(y)),years-years.mean()]);K=rho**np.abs(years[:,None]-years[None,:])
    def core(q,ret=False):
        C=S+q*K+np.eye(len(y))*1e-9;cf=cho_factor(C,lower=True,check_finite=False);ix=cho_solve(cf,X,check_finite=False);bc=np.linalg.inv(X.T@ix);b=bc@(X.T@cho_solve(cf,y,check_finite=False));res=y-X@b
        score=2*np.log(np.diag(cf[0])).sum()+np.linalg.slogdet(X.T@ix)[1]+res@cho_solve(cf,res,check_finite=False)
        return (b,bc,C,res,cf) if ret else score
    top=max(1.,10*np.var(y));opt=minimize_scalar(core,bounds=(0,top),method='bounded',options={'xatol':1e-8});q=float(opt.x)
    if core(0)<=opt.fun:q=0.
    b,bc,C,res,cf=core(q,True)
    whiten=np.linalg.solve(np.linalg.cholesky(C),res)
    return dict(estimate=b[1],se=np.sqrt(bc[1,1]),process_variance=q,beta=b,cov=C,residual=res,whitened=whiten,success=bool(opt.success),boundary=q<1e-7)

def run(args):
    if args.pool_df<=0 or args.draws<20 or not 1<=args.parameter_bootstraps<=args.draws or args.trend_bootstraps<10: raise ValueError('Invalid pooling/draw/bootstrap arguments')
    source_hash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    O=Path(args.output);O.mkdir(parents=True,exist_ok=True);(O/'private').mkdir(exist_ok=True)
    started=datetime.now(timezone.utc).isoformat()
    (O/'run_status.json').write_text(json.dumps(dict(status='running',started=started,source_sha256=source_hash)))
    for old in ['annual_observation_covariance.csv']:
        if (O/old).exists():
            (O/'archive_previous_labels').mkdir(exist_ok=True)
            (O/old).replace(O/'archive_previous_labels'/old)
    rng=np.random.default_rng(args.seed);N=args.draws;used=[]
    def read(path,**kw):
        p=Path(path);used.append(p);return pd.read_csv(p,**kw)
    visits=read(args.visits);visits=visits[visits.tolerance_days.eq(args.tolerance)].copy() if 'tolerance_days' in visits else visits.copy()
    assert len(visits)>0 and np.isfinite(visits.density).all() and visits.density.ge(0).all()
    visits=visits.reset_index(drop=True);visits['use_id']=np.arange(len(visits));visits['date']=visits.date.astype(str).str[:10]
    canonical=read('outputs/pisco_biomass_borrowing_20261102/visit_biomass_support.csv')
    visit_keys=['series_id','site','status','year','zone','date']
    assert not visits.duplicated(visit_keys).any() and set(map(tuple,visits[visit_keys].values))==set(map(tuple,canonical[visit_keys].values)), 'Incomplete or changed visit roster'
    assert visits.loc[visits.method.eq('unsupported'),'mass'].isna().all()
    assert visits.loc[visits.density.gt(0)&visits.method.ne('unsupported'),'mass'].gt(0).all()
    assert np.isfinite(visits.loc[visits.method.ne('unsupported'),'mass']).all()
    assert visits.loc[visits.density.eq(0),'biomass_point'].eq(0).all()
    samples=read('outputs/pisco_timing_donors_20261102/size_sample_inventory.csv');samples['date']=samples.date.astype(str).str[:10]
    raw=read('/Users/adrianstier/sbc-donham-kelp-mpa-data/data/ALL_sizefreq_2024.csv',low_memory=False)
    raw=raw[raw.campus.isin(['VRG','UCSB']) & raw.classcode.isin(['MESFRA','STRPUR','PANINT'])].copy();raw['Taxa']=raw.classcode.map({'MESFRA':'M. franciscanus','STRPUR':'S. purpuratus','PANINT':'P. interruptus'});raw=raw[raw['count'].gt(0)&raw['size'].gt(0)&(raw.Taxa.eq('P. interruptus')|raw['size'].ge(25))];raw['date']=pd.to_datetime(raw[['year','month','day']],errors='coerce').dt.strftime('%Y-%m-%d')
    raw['mass']=np.select([raw.Taxa.eq('M. franciscanus'),raw.Taxa.eq('S. purpuratus')],[.00059*raw['size']**2.917,.00059*raw['size']**2.870],default=.001352821*raw['size']**2.913963113)
    keys=['site','year','zone','date','Taxa'];records=[]
    for key,g in raw.merge(samples[keys].drop_duplicates(),on=keys,validate='many_to_one').groupby(keys):
        n=g['count'].sum();mean=np.average(g.mass,weights=g['count']);ss=np.sum(g['count']*(g.mass-mean)**2)
        records.append(dict(zip(keys,key),n=n,mean_mass_raw=mean,relative_ss=ss/mean**2,df=max(n-1,0),sample_id_physical=stable_id(key)))
    inv=pd.DataFrame(records);samples=samples.merge(inv,on=keys,validate='one_to_one');assert np.allclose(samples.mean_mass,samples.mean_mass_raw)
    samples['pooled_cv2']=samples.groupby(['series_id','status']).relative_ss.transform('sum')/samples.groupby(['series_id','status']).df.transform('sum')
    # Explicit partial-pooling assumption, never treating singleton variance as measured.
    samples['sampling_cv2']=(samples.relative_ss+args.pool_df*samples.pooled_cv2)/(samples.df+args.pool_df)
    samples['sampling_log_mean_variance']=np.log1p(samples.sampling_cv2/samples.n)
    assert np.isfinite(samples.sampling_log_mean_variance).all() and samples.sampling_log_mean_variance.gt(0).all()
    samples['variance_basis']=np.where(samples.n.eq(1),'pooled species/MPA/role dispersion; no within-sample estimate','sample dispersion shrunk to species/MPA/role pool')
    samples.to_csv(O/'size_observation_registry.csv',index=False)
    # Held-out site prediction checks: pool dispersion only from training sites.
    # Groups with fewer than three measured sites are explicitly unvalidated.
    cvrows=[]
    for (sid,role),sg in samples.groupby(['series_id','status']):
        if sg.site.nunique()<3:
            cvrows.append(dict(series_id=sid,status=role,available=False,reason='fewer than three measured sites; no independent protected-site validation'))
            continue
        for held in sorted(sg.site.unique()):
            tr=sg[sg.site.ne(held)].copy();te=sg[sg.site.eq(held)].copy();pool=tr.relative_ss.sum()/tr.df.sum()
            vt=np.log1p((tr.relative_ss+args.pool_df*pool)/(tr.df+args.pool_df)/tr.n)
            vte=np.log1p((te.relative_ss+args.pool_df*pool)/(te.df+args.pool_df)/te.n)
            center=float(tr.year.mean());xx=design(tr,center);xt=design(te,center);yy=np.log(tr.mean_mass.to_numpy())+.5*vt.to_numpy()
            ff=fit_reml(yy,xx,vt.to_numpy(),kernels(tr,tr));assert ff['success'],ff['optimizer_message'];mm,ll=prediction(ff,yy,xx,xt,kernels(te,tr),kernels(te,te));sd=np.sqrt(np.sum(ll**2,axis=1)+vte.to_numpy());mm=mm-.5*vte.to_numpy()
            for i,t in enumerate(te.itertuples()):
                cvrows.append(dict(series_id=sid,status=role,available=True,held_site=held,year=t.year,zone=t.zone,date=t.date,observed=t.mean_mass,predicted_median=np.exp(mm[i]),low=np.exp(mm[i]-1.96*sd[i]),high=np.exp(mm[i]+1.96*sd[i]),covered=bool(abs(np.log(t.mean_mass)-mm[i])<=1.96*sd[i]),converged=ff['success'],reason='Conditional plug-in heldout observed-mean interval; dependent samples within sites'))
    pd.DataFrame(cvrows).to_csv(O/'size_model_site_validation.csv',index=False)

    masses=np.full((len(visits),N),np.nan);model_rows=[];prediction_rows=[];fitparams=[]
    for (sid,role),s in samples.groupby(['series_id','status']):
        s=s.reset_index(drop=True);target=visits[(visits.series_id==sid)&(visits.status==role)&visits.density.gt(0)&visits.mass.notna()].copy();center=float(s.year.mean())
        X=design(s,center);xt=design(target,center);known=s.sampling_log_mean_variance.to_numpy();y=np.log(s.mean_mass.to_numpy())+.5*known;K=kernels(s,s);cross=kernels(target,s);tt=kernels(target,target)
        # One station cannot identify station intercept variance separately from fixed intercept.
        if s.site.nunique()==1:K[0]=np.zeros_like(K[0]);cross[0]=np.zeros_like(cross[0]);tt[0]=np.zeros_like(tt[0])
        f=fit_reml(y,X,known,K);assert f['success'],f['optimizer_message'];mu,L=prediction(f,y,X,xt,cross,tt)
        chunks=np.array_split(np.arange(N),args.parameter_bootstraps);boot_ok=0;bootstrap_failures=0
        observation_chol=np.linalg.cholesky(f['cov'])
        for bi,cols in enumerate(chunks):
            sim=X@f['beta']+observation_chol@rng.normal(size=len(s));fb=fit_reml(sim,X,known,K,f['theta'])
            attempt_success=fb['success'];attempt_message=fb['optimizer_message']
            if not attempt_success:
                bootstrap_failures+=1
                fb=fit_reml(sim,X,known,K)
            if not fb['success']: raise RuntimeError('Bootstrap fit failed after retry: '+fb['optimizer_message'])
            boot_ok+=1
            # Recondition original observations using bootstrap covariance parameters.
            # Fixed effects re-estimated by GLS; no second draw of donor sampling residuals.
            C=np.diag(known)+sum(t*k for t,k in zip(fb['theta'],K))+np.eye(len(s))*1e-9;cf=cho_factor(C,lower=True);vx=cho_solve(cf,X);bc=np.linalg.inv(X.T@vx);bf=bc@(X.T@cho_solve(cf,y));ff=dict(theta=fb['theta'],cf=cf,beta=bf,beta_cov=bc)
            mb,lb=prediction(ff,y,X,xt,cross,tt);masses[target.use_id.to_numpy()[:,None],cols]=np.exp(mb[:,None]+lb@rng.normal(size=(len(target),len(cols))))
            fitparams.append(dict(series_id=sid,status=role,bootstrap=bi,attempt_converged=attempt_success,attempt_message=attempt_message,retried=not attempt_success,converged=fb['success'],used_fallback=False,site_variance=fb['theta'][0],siteyear_variance=fb['theta'][1],event_variance=fb['theta'][2]))
        resid=y-X@f['beta'];wh=np.linalg.solve(np.linalg.cholesky(f['cov']),resid)
        model_rows.append(dict(series_id=sid,status=role,n_samples=len(s),sites=s.site.nunique(),years=s.year.nunique(),singletons=int(s.n.eq(1).sum()),converged=f['success'],boundary_components=f['boundary'],site_variance=f['theta'][0],siteyear_variance=f['theta'][1],event_variance=f['theta'][2],bootstrap_success=boot_ok,bootstrap_retries=bootstrap_failures,standardized_residual_max_abs=float(abs(wh).max()),sampling_variance_pool_df=args.pool_df,scope='Conditional latent-size diagnostic; borrowing eligibility unchanged'))
        for i,t in enumerate(target.itertuples()):prediction_rows.append(dict(use_id=t.use_id,series_id=sid,site=t.site,year=t.year,zone=t.zone,date=t.date,status=role,method=t.method,empirical_mass=t.mass,model_median_mass=float(np.exp(mu[i])),model_draw_mean=float(masses[t.use_id].mean()),model_draw_sd=float(masses[t.use_id].std(ddof=1))))
    masses[visits.density.eq(0)]=0
    # Count-unit weights shared across taxa and every use of a physical transect/event.
    counts=read(args.count_units);counts['date']=counts.date.astype(str).str[:10]
    required=['series_id','site','year','zone','date','method','transect','density']
    assert set(required)<=set(counts),counts.columns.tolist()
    counts['series_id']=counts.series_id.str.replace(r'\|Den$','|Bio',regex=True)
    phys=['site','year','zone','date','method','transect'];units=counts[phys].drop_duplicates().copy();units['unit_id']=units.apply(lambda r:stable_id(r.tolist()),axis=1)
    assert units.unit_id.is_unique;units['idx']=np.arange(len(units));counts=counts.merge(units,on=phys,validate='many_to_one')
    w=rng.exponential(size=(len(units),N));densitydraw=np.full((len(visits),N),np.nan);checks=[]
    for t in visits.itertuples():
        c=counts[(counts.series_id==t.series_id)&(counts.site==t.site)&(counts.year==t.year)&(counts.zone==t.zone)&(counts.date==t.date)]
        assert len(c)>0
        # Match original method means if more than one method appears in a dated zone.
        per=[];means=[]
        for _,cg in c.groupby('method'):
            ww=w[cg.idx];per.append(cg.density.to_numpy()@ww/ww.sum(axis=0));means.append(cg.density.mean())
        check=float(np.mean(means));assert np.isclose(check,t.density), (t.use_id,check,t.density)
        densitydraw[t.use_id]=np.mean(per,axis=0);checks.append(dict(use_id=t.use_id,expected=t.density,reconstructed=check,physical_units=len(c)))
    units.drop(columns='idx').to_csv(O/'count_observation_registry.csv',index=False);pd.DataFrame(checks).to_csv(O/'count_reconstruction.csv',index=False)
    draws=densitydraw*masses;pd.DataFrame(model_rows).to_csv(O/'size_model_diagnostics.csv',index=False);pd.DataFrame(fitparams).to_csv(O/'size_parameter_bootstraps.csv',index=False);pd.DataFrame(prediction_rows).to_csv(O/'mass_predictions.csv',index=False)
    off=read('outputs/all_program_station_refits_20261102/offsets.csv').set_index('series_id').offset
    base=read('outputs/pisco_biomass_borrowing_20261102/biomass_annual.csv');annual=[];drawrows=[]
    for sid,g in visits.groupby('series_id'):
        for year,z in g.groupby('year'):
            if not np.isfinite(draws[z.use_id]).all():continue
            role={}
            for status,rg in z.groupby('status'):
                sitevals=[]
                for site,sg in rg.groupby('site'):
                    zonevals=[draws[zz.use_id].mean(axis=0) for _,zz in sg.groupby('zone')];sitevals.append(np.mean(zonevals,axis=0))
                role[status]=np.mean(sitevals,axis=0)
            assert set(role)=={'mpa','reference'}
            y=np.log((role['mpa']+off[sid])/(role['reference']+off[sid]));t=base[(base.series_id==sid)&(base.year==year)].t.iloc[0]
            annual.append(dict(series_id=sid,year=year,t=t,estimate=y.mean(),latent_draw_sd=y.std(ddof=1),point_definition='mean conditional latent log-ratio draw'));drawrows.append(y)
    aa=pd.DataFrame(annual);Y=np.array(drawrows);aa.to_csv(O/'annual_observations.csv',index=False);np.savez_compressed(O/'private/annual_draws.npz',draws=Y)
    if len(aa)<2: raise ValueError('Fewer than two supported annual observations; covariance unavailable')
    C=np.cov(Y);covrows=[]
    for i,r in aa.iterrows():
        for j,sr in aa.iterrows():covrows.append(dict(series_id=r.series_id,year=r.year,other_series_id=sr.series_id,other_year=sr.year,latent_covariance=C[i,j]))
    pd.DataFrame(covrows).to_csv(O/'annual_latent_covariance.csv',index=False)
    # Slopes are finite-observed-year functionals of joint latent draws, not a new process likelihood.
    functions=[]
    for sid,g in aa[aa.t>0].groupby('series_id'):
        if g.year.nunique()<2: continue
        yy=Y[g.index];x=np.column_stack([np.ones(len(g)),g.year-g.year.mean()]);b=(np.linalg.pinv(x)@yy)[1]
        functions.append(dict(series_id=sid,n_years=len(g),observed_years=";".join(map(str,g.year)),mean_slope=b.mean(),draw_low=np.quantile(b,.025),draw_high=np.quantile(b,.975),scope='Conditional uncertainty in observed-year slope; no extra temporal process CI'))
    pd.DataFrame(functions).to_csv(O/'observed_year_slope_draws.csv',index=False)
    # Separate working measurement likelihood, only years with all direct size samples.
    # These are unsmoothed empirical estimates; never recycle latent posterior covariance.
    direct_mass=np.full((len(visits),N),np.nan);sample_draw={}
    for r in samples.itertuples():
        key=(r.site,r.year,r.zone,r.date,r.Taxa)
        sample_draw[key]=r.mean_mass*np.exp(-.5*r.sampling_log_mean_variance+np.sqrt(r.sampling_log_mean_variance)*rng.normal(size=N))
    for t in visits.itertuples():
        if t.density==0:direct_mass[t.use_id]=0
        elif t.method=='measured_same_visit':
            key=(t.site,t.year,t.zone,t.date,t.series_id.split('|')[1]);direct_mass[t.use_id]=sample_draw[key]
            ss=samples[(samples.site==t.site)&(samples.year==t.year)&(samples.zone==t.zone)&(samples.date==t.date)&(samples.Taxa==t.series_id.split('|')[1])]
            assert len(ss)==1 and np.isclose(t.biomass_point,t.density*ss.mean_mass.iloc[0])
    direct_draw=densitydraw*direct_mass;results=[];diag=[];bootrows=[];direct_annual=[];admission=[];direct_covrows=[]
    for sid,g in visits.groupby('series_id'):
        years=[];points=[];replicates=[]
        for year,z in g.groupby('year'):
            tt=base[(base.series_id==sid)&(base.year==year)].t.iloc[0]
            if tt<=0 or not np.isfinite(direct_draw[z.use_id]).all():continue
            roles={};rp={}
            for status,rg in z.groupby('status'):
                vals=[];pts=[]
                for _,sg in rg.groupby('site'):
                    vals.append(np.mean([direct_draw[zz.use_id].mean(axis=0) for _,zz in sg.groupby('zone')],axis=0))
                    pts.append(np.mean([zz.biomass_point.mean() for _,zz in sg.groupby('zone')]))
                roles[status]=np.mean(vals,axis=0);rp[status]=np.mean(pts)
            years.append(year);points.append(np.log((rp['mpa']+off[sid])/(rp['reference']+off[sid])));replicates.append(np.log((roles['mpa']+off[sid])/(roles['reference']+off[sid])))
        admission.append(dict(series_id=sid,direct_post_years=len(years),years=';'.join(map(str,years)),process_diagnostic_admitted=len(years)>=4,reason='At least four direct-measurement years' if len(years)>=4 else 'Too few direct-measurement years; retain count-led and conditional borrowing evidence'))
        for yr,pt,dr in zip(years,points,replicates):direct_annual.append(dict(series_id=sid,year=yr,observed_log_ratio=pt,bootstrap_sd=np.std(dr,ddof=1),sampling_log_bias=float(np.mean(dr)-pt),sampling_skew=float(skew(dr))))
        if len(years)<4:continue
        years=np.array(years);yy=np.array(points);S=np.cov(np.array(replicates))
        for i,yr in enumerate(years):
            for j,yr2 in enumerate(years):direct_covrows.append(dict(series_id=sid,year=yr,other_year=yr2,sampling_covariance=S[i,j]))
        for rho in [0.,.3,.6]:
            f=fit_trend(yy,years,S,rho);assert f['success'];bs=[]
            for j in range(args.trend_bootstraps):
                sim=np.column_stack([np.ones(len(yy)),years-years.mean()])@f['beta']+np.linalg.cholesky(f['cov'])@rng.normal(size=len(yy))
                ff=fit_trend(sim,years,S,rho);assert ff['success'];bs.append(ff['estimate']);bootrows.append(dict(series_id=sid,rho=rho,draw=j,slope=ff['estimate']))
            dif=np.array(bs)-f['estimate'];lo=f['estimate']-np.quantile(dif,.975);hi=f['estimate']-np.quantile(dif,.025)
            loo=[]
            if len(yy)>4:
                for j in range(len(yy)):
                    keep=np.arange(len(yy))!=j;ff=fit_trend(yy[keep],years[keep],S[np.ix_(keep,keep)],rho);assert ff['success'];loo.append(ff['estimate'])
            results.append(dict(series_id=sid,rho=rho,n_years=len(yy),estimate=f['estimate'],conditional_gls_se=f['se'],bootstrap_low=lo,bootstrap_high=hi,process_variance=f['process_variance'],process_boundary=f['boundary'],covariance_condition=np.linalg.cond(f['cov']),process_upper_bound_warning=f['process_variance']>.99*max(1.,10*np.var(yy)),leave_one_year_min=min(loo) if loo else np.nan,leave_one_year_max=max(loo) if loo else np.nan,interval_scope='Direct-measurement-only conditional working Gaussian likelihood; estimated sampling covariance held fixed; not total uncertainty',trend_bootstraps=args.trend_bootstraps))
            for k,year in enumerate(years):diag.append(dict(series_id=sid,rho=rho,year=year,residual=f['residual'][k],whitened_residual=f['whitened'][k]))
    pd.DataFrame(direct_covrows).to_csv(O/'direct_sampling_covariance.csv',index=False)
    pd.DataFrame(admission).to_csv(O/'direct_process_admission.csv',index=False);pd.DataFrame(direct_annual).to_csv(O/'direct_annual_observations.csv',index=False)
    pd.DataFrame(results).to_csv(O/'direct_trend_scenarios.csv',index=False);pd.DataFrame(diag).to_csv(O/'direct_trend_residuals.csv',index=False);pd.DataFrame(bootrows).to_csv(O/'private/direct_trend_bootstraps.csv',index=False)
    used.append(Path(__file__))
    (O/'input_hashes.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in used},indent=2))
    policy=dict(draws=N,seed=args.seed,tolerance_days=args.tolerance,within_sample_pool_df=args.pool_df,size_model='log sample mean mass + half sampling-log-variance ~ intercept + MID zone + calendar trend; random site, site-year and event; known partially pooled sampling variance',size_parameter_bootstraps=args.parameter_bootstraps,trend_bootstraps=args.trend_bootstraps,process_correlation_scenarios=[0,.3,.6],limits=['Exploratory hierarchical assumptions, not validated protected borrowing','No new supported biomass cells','Pooled within-sample variance treated as known','Allometry/offsets fixed','Independent counts and sizes conditional on measurements','Physical count weights shared across species; size model species/MPA/role-specific','Cross-year transect persistence not established; count weights event-specific','No process likelihood fitted to latent posterior summaries; direct-only working likelihood separate','Total uncertainty not claimed'])
    assert hashlib.sha256(Path(__file__).read_bytes()).hexdigest()==source_hash, 'Source changed during execution'
    (O/'run_status.json').write_text(json.dumps(dict(status='complete',started=started,finished=datetime.now(timezone.utc).isoformat(),source_sha256=source_hash),indent=2))
    (O/'policy.json').write_text(json.dumps(policy,indent=2));print(json.dumps(dict(size_models=len(model_rows),annual_observations=len(aa),post_comparison_years=int((aa.t>0).sum()),trend_scenarios=len(results),count_checks=len(checks),singleton_size_samples=int(samples.n.eq(1).sum())),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--visits',required=True);p.add_argument('--count-units',required=True);p.add_argument('--output',required=True);p.add_argument('--tolerance',type=int,default=7);p.add_argument('--draws',type=int,default=400);p.add_argument('--parameter-bootstraps',type=int,default=40);p.add_argument('--trend-bootstraps',type=int,default=300);p.add_argument('--pool-df',type=float,default=2.);p.add_argument('--seed',type=int,default=20261102);run(p.parse_args())
