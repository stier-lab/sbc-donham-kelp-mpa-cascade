"""Joint working Gaussian likelihood for raw size means and annual count contributions.
First-order log-biomass aggregation; conditional sampling variances and allometry.
No latent posterior summaries are observations. Review only, fixed support.
"""
from pathlib import Path
import argparse,json,hashlib,os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor,cho_solve,block_diag
from scipy.optimize import minimize
from biomass_observation_model import design,kernels
BASE=Path('outputs/biomass_next_steps_20261102')


def fit(y,X,S,K,start=None):
    n,p=X.shape
    if np.linalg.matrix_rank(X)!=p or n<=p:raise ValueError('Rank or residual df failure')
    Ks=np.array(K); upper=30.
    def core(q,ret=False):
        C=S+np.einsum('k,kij->ij',q,Ks)+np.eye(n)*1e-9
        cf=cho_factor(C,lower=True,check_finite=False)
        inv=cho_solve(cf,np.eye(n),check_finite=False);ix=inv@X;bc=np.linalg.inv(X.T@ix);b=bc@(ix.T@y);P=inv-ix@bc@ix.T;py=P@y
        obj=2*np.log(np.diag(cf[0])).sum()+2*np.log(np.diag(np.linalg.cholesky(X.T@ix))).sum()+y@py
        grad=np.einsum('ij,kji->k',P,Ks)-np.einsum('i,kij,j->k',py,Ks,py)
        if ret:return dict(beta=b,bc=bc,C=C,P=P,inv=inv,obj=float(obj),q=q,gradient=grad)
        return obj,grad
    bounds=[(0,upper) if np.any(k) else (0,0) for k in K]
    starts=[np.full(len(K),.1) if start is None else start]
    opts=[]
    for x in starts:
        opts.append(minimize(core,x,jac=True,method='L-BFGS-B',bounds=bounds,options={'maxiter':400,'ftol':1e-10,'gtol':1e-6,'maxls':50}))
    opt=min(opts,key=lambda x:x.fun)
    attempts=1
    if not opt.success:
        attempts+=1
        opt=minimize(core,np.full(len(K),.2),jac=True,method='L-BFGS-B',bounds=bounds,options={'maxiter':800,'ftol':1e-10,'maxls':50})
    if not opt.success:
        attempts+=1
        opt=minimize(core,np.clip(opt.x,0,upper),jac=True,method='SLSQP',bounds=bounds,options={'maxiter':1200,'ftol':1e-9})
    if not opt.success:raise RuntimeError(str(opt.message))
    f=core(opt.x,True);f.update(success=True,optimizer_attempts=attempts,iterations=int(opt.nit),upper_boundary=bool(np.any(opt.x>29.99)))
    return f


def aggregation(v,years,mass,offset):
    """Exact equal-date, equal-zone, equal-site arithmetic biomass and Jacobian."""
    A=np.zeros((len(years),2,len(v)))
    for j,yr in enumerate(years):
        for ri,role in enumerate(['mpa','reference']):
            g=v[(v.year==yr)&(v.status==role)]
            assert len(g)>0
            for _,sg in g.groupby('site'):
                for _,zg in sg.groupby('zone'):
                    A[j,ri,zg.index]=1/(g.site.nunique()*sg.zone.nunique()*len(zg))
    bm=v.density.to_numpy()*mass
    total=np.einsum('jri,i->jr',A,bm)
    h=np.log((total[:,0]+offset)/(total[:,1]+offset))
    W=A[:,0,:]*bm/(total[:,0,None]+offset)-A[:,1,:]*bm/(total[:,1,None]+offset)
    return h,W,A


def count_cov(v,A,mass,offset,counts,rng,draws=1500):
    density=np.empty((len(v),draws));units=sorted(counts.physical_unit_id.unique());uid={u:i for i,u in enumerate(units)}
    weights=rng.exponential(size=(len(units),draws))
    for i,t in v.iterrows():
        cc=counts[(counts.site==t.site)&(counts.year==t.year)&(counts.zone==t.zone)&(counts.date==t.date)]
        assert len(cc)>0
        gs=[];pts=[]
        for _,g in cc.groupby('method'):
            ww=weights[[uid[u] for u in g.physical_unit_id]]
            gs.append(g.density.to_numpy()@ww/ww.sum(axis=0));pts.append(g.density.mean())
        assert np.isclose(np.mean(pts),t.density)
        density[i]=np.mean(gs,axis=0)
    total=np.einsum('jri,id->jrd',A,density*mass[:,None])
    rr=np.log((total[:,0,:]+offset)/(total[:,1,:]+offset))
    return np.cov(rr),rr


def build(s,v,years,offset,counts,rng,rho=0.,expansion=None,pool_df=2):
    if not 0<=rho<1:raise ValueError('rho must be in [0,1)')
    if len(set(years))!=len(years) or len(years)<4:raise ValueError('At least four unique years required')
    keys=['site','status','year','zone','date']
    if s.duplicated(keys).any():raise ValueError('Duplicate physical size sample')
    if v.duplicated(keys).any():raise ValueError('Duplicate visit use')
    if 'method' in v and v.method.eq('unsupported').any():raise ValueError('Unsupported biomass is not eligible')
    v=v.reset_index(drop=True);s=s.reset_index(drop=True);ns=len(s);nt=len(v);ny=len(years)
    mass=v.mass.to_numpy().copy() if expansion is None else expansion.copy()
    mass[v.density.eq(0)]=0
    assert np.isfinite(mass).all() and np.all(mass[v.density.gt(0)]>0)
    logm=np.log(np.where(mass>0,mass,1))
    h,W,A=aggregation(v,years,mass,offset)
    S_count,rep=count_cov(v,A,mass,offset,counts,rng)
    # Sampling variance for mean size from original grouped frequencies, explicitly pooled.
    pool=s.groupby('status').relative_ss.transform('sum')/s.groupby('status').df.transform('sum')
    known=np.log1p((s.relative_ss+pool_df*pool)/(s.df+pool_df)/s.n).to_numpy()
    ys=np.log(s.mean_mass.to_numpy())+.5*known
    # All physical size samples enter once. Targets sharing samples inherit the same latent process.
    Xs=np.zeros((ns,6));Xt=np.zeros((nt,6));Kobs=[];Kcross=[];Ktt=[];labels=[]
    for ri,role in enumerate(['mpa','reference']):
        a=np.flatnonzero(s.status.eq(role));b=np.flatnonzero(v.status.eq(role));ss=s.iloc[a];vv=v.iloc[b];center=float(ss.year.mean())
        Xs[np.ix_(a,range(ri*3,ri*3+3))]=design(ss,center);Xt[np.ix_(b,range(ri*3,ri*3+3))]=design(vv,center)
        for kk,(ko,kc,kt) in enumerate(zip(kernels(ss,ss),kernels(vv,ss),kernels(vv,vv))):
            oo=np.zeros((ns,ns));cc=np.zeros((nt,ns));tt=np.zeros((nt,nt))
            if not(kk==0 and ss.site.nunique()==1):oo[np.ix_(a,a)]=ko;cc[np.ix_(b,a)]=kc;tt[np.ix_(b,b)]=kt
            Kobs.append(oo);Kcross.append(cc);Ktt.append(tt);labels.append(role+'_'+['site','siteyear','event'][kk])
    Z=np.column_stack([np.ones(ny),(years-years.mean())])
    X=np.block([[Xs,np.zeros((ns,2))],[-W@Xt,Z]])
    # c_obs = X_beta beta - W z + process + count sampling error.
    y=np.r_[ys,h-W@logm]
    S=block_diag(np.diag(known),S_count)
    K=[np.block([[ko,-kc.T@W.T],[-W@kc,W@kt@W.T]]) for ko,kc,kt in zip(Kobs,Kcross,Ktt)]
    process=np.zeros((ns+ny,ns+ny));process[ns:,ns:]=rho**np.abs(years[:,None]-years[None,:]);K.append(process);labels.append('annual_process')
    return dict(y=y,X=X,S=S,K=K,labels=labels,ns=ns,W=W,Xt=Xt,Xs=Xs,Kcross=Kcross,Ktt=Ktt,h=h,logm=logm,mass=mass,A=A,known=known,count_cov=S_count,count_replicates=rep)


def latent(d,f):
    G=sum(q*k for q,k in zip(f['q'][:6],d['Ktt']))
    cross=np.column_stack([sum(q*k for q,k in zip(f['q'][:6],d['Kcross'])), -G@d['W'].T])
    T=np.column_stack([d['Xt'],np.zeros((len(d['Xt']),2))])
    residual=d['y']-d['X']@f['beta'];mu=T@f['beta']+cross@f['inv']@residual
    delta=T-cross@f['inv']@d['X'];C=G-cross@f['inv']@cross.T+delta@f['bc']@delta.T
    return mu,(C+C.T)/2


def run(args):
    if args.bootstraps<19 or args.pool_df<=0:raise ValueError('Invalid bootstrap or pooling settings')
    if args.tolerance not in [0,3,7,14]:raise ValueError('Unsupported tolerance scenario')
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True);(out/'fits').mkdir(exist_ok=True)
    marker=out/'completed.json'
    if marker.exists():marker.unlink()
    paths=[BASE/'uncertainty/size_observation_registry.csv',BASE/'borrowing/visit_support.csv',BASE/'density/biomass_roster_count_units.csv',Path('outputs/all_program_station_refits_20261102/offsets.csv'),BASE/'borrowing/annual.csv',Path(__file__)]
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    visits=pd.read_csv(paths[1]);visits=visits[visits.tolerance_days.eq(args.tolerance)].copy();samples=pd.read_csv(paths[0]);counts=pd.read_csv(paths[2]);counts['series_id']=counts.biomass_series_id
    offsets=pd.read_csv(paths[3]).set_index('series_id').offset;ann=pd.read_csv(paths[4]);ann=ann[(ann.tolerance_days==args.tolerance)&(ann.t>0)&ann.supported]
    results=[];params=[];bootrows=[];annual=[];curvature=[];influence=[];fitindex=[];resids=[];support=[]
    for si,(sid,g) in enumerate(visits.groupby('series_id')):
        if args.series and args.series not in sid:continue
        years=ann.loc[ann.series_id.eq(sid),'year'].to_numpy(dtype=int);years=np.sort(years)
        assert len(years)>=4
        v=g[g.year.isin(years)].reset_index(drop=True)
        assert v.mass.notna().all() and not v.method.eq('unsupported').any()
        s=samples[samples.series_id.eq(sid)].reset_index(drop=True);cc=counts[counts.series_id.eq(sid)].copy();off=offsets[sid]
        support.append(dict(series_id=sid,n_years=len(years),years=';'.join(map(str,years)),size_samples=len(s),visit_uses=len(v),borrowed_uses=int((~v.method.isin(['measured_same_visit','zero_count'])).sum()),methods=';'.join(sorted(v.method.unique()))))
        for rho in args.rhos:
            rng=np.random.default_rng(args.seed+si*1000+int(rho*100));d=build(s,v,years,off,cc,np.random.default_rng(args.seed+si*1000),rho,pool_df=args.pool_df)
            f=fit(d['y'],d['X'],d['S'],d['K']);key=f's{si}_rho{rho:g}';ny=len(years);ns=len(s)
            np.savez_compressed(out/'fits'/f'{key}.npz',y=d['y'],X=d['X'],S=d['S'],K=np.array(d['K']),q=f['q'],beta=f['beta'],beta_cov=f['bc'],years=years)
            fitindex.append(dict(key=key,series_id=sid,rho=rho,ns=ns,n_years=ny))
            chol=np.linalg.cholesky(f['C']);bs=[];ses=[]
            for bi in range(args.bootstraps):
                yy=d['X']@f['beta']+chol@rng.normal(size=len(d['y']));fb=fit(yy,d['X'],d['S'],d['K'],start=f['q'])
                bs.append(fb['beta'][-1]);ses.append(np.sqrt(fb['bc'][-1,-1]));bootrows.append(dict(series_id=sid,rho=rho,replicate=bi,slope=bs[-1],se=ses[-1],converged=fb['success'],optimizer_attempts=fb['optimizer_attempts'],upper_boundary=fb['upper_boundary']))
            se=np.sqrt(f['bc'][-1,-1]);pivot=(np.array(bs)-f['beta'][-1])/np.array(ses)
            lo=f['beta'][-1]-np.quantile(pivot,.975)*se;hi=f['beta'][-1]-np.quantile(pivot,.025)*se
            # Keep normal and basic bootstrap alternatives visible for sparse histories.
            dif=np.array(bs)-f['beta'][-1]
            results.append(dict(series_id=sid,rho=rho,n_years=ny,estimate=f['beta'][-1],se=se,low=lo,high=hi,basic_low=f['beta'][-1]-np.quantile(dif,.975),basic_high=f['beta'][-1]-np.quantile(dif,.025),process_variance=f['q'][-1],boundary_components=int((f['q']<1e-7).sum()),upper_boundary=f['upper_boundary'],converged=f['success'],bootstraps=args.bootstraps,condition=np.linalg.cond(f['C']),interval_scope='Studentized parametric bootstrap; joint first-order working likelihood; sampling variances/allometry/expansion fixed'))
            for label,q in zip(d['labels'],f['q']):params.append(dict(series_id=sid,rho=rho,component=label,variance=q))
            mu,C=latent(d,f);e,Q=np.linalg.eigh(C);assert min(e)>-1e-6
            draws=mu[:,None]+(Q*np.sqrt(np.maximum(e,0)))@rng.normal(size=(len(mu),300))
            total=np.einsum('jri,id->jrd',d['A'],v.density.to_numpy()[:,None]*np.exp(draws));exact=np.log((total[:,0]+off)/(total[:,1]+off));linear=(d['h']-d['W']@d['logm'])[:,None]+d['W']@draws
            for j,yr in enumerate(years):
                curvature.append(dict(series_id=sid,rho=rho,year=yr,mean_exact_minus_linear=np.mean(exact[j]-linear[j]),q95_abs_error=np.quantile(abs(exact[j]-linear[j]),.95)))
                annual.append(dict(series_id=sid,rho=rho,year=yr,empirical_ratio=d['h'][j],joint_linear_latent_ratio=(d['h']-d['W']@d['logm']+d['W']@mu)[j],fitted_process=f['beta'][-2]+(yr-years.mean())*f['beta'][-1]))
            white=np.linalg.solve(chol,d['y']-d['X']@f['beta'])
            for j,val in enumerate(white):resids.append(dict(series_id=sid,rho=rho,row=j,kind='size' if j<ns else 'annual count contribution',whitened_residual=val))
            if rho==0:
                # Recenter the first-order expansion to fitted sizes, on identical raw observations.
                expmass=np.exp(mu);expmass[v.density.eq(0)]=0
                dd=build(s,v,years,off,cc,np.random.default_rng(args.seed+si*1000),rho,expmass,args.pool_df);ff=fit(dd['y'],dd['X'],dd['S'],dd['K'])
                results[-1]['recentered_estimate']=ff['beta'][-1]
                for j,yr in enumerate(years):
                    keep=np.arange(ns+ny)!=ns+j;fj=fit(d['y'][keep],d['X'][keep],d['S'][np.ix_(keep,keep)],[k[np.ix_(keep,keep)] for k in d['K']],start=f['q'])
                    influence.append(dict(series_id=sid,omitted_count_year=yr,estimate=fj['beta'][-1],difference=fj['beta'][-1]-f['beta'][-1],scope='Remove annual count equation only; original size observations retained'))
            print(sid,rho,results[-1]['estimate'],lo,hi,flush=True)
            # Save progress without marking completion.
            pd.DataFrame(results).to_csv(out/'joint_trends.csv',index=False)
    for name,rows in [('support',support),('variance_components',params),('bootstrap_fits',bootrows),('annual_comparison',annual),('linearization_checks',curvature),('leave_one_count_year',influence),('fit_index',fitindex),('residuals',resids)]:pd.DataFrame(rows).to_csv(out/f'{name}.csv',index=False)
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items())
    (out/'input_hashes.json').write_text(json.dumps(hashes,indent=2));(out/'policy.json').write_text(json.dumps(vars(args),indent=2))
    (out/'completed.json').write_text(json.dumps(dict(review_only=True,model='joint first-order Gaussian likelihood',fits=len(results),joint_process_with_borrowing_implemented=True,total_uncertainty_complete=False),indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',default='outputs/joint_borrowed_biomass_20261102');p.add_argument('--tolerance',type=int,default=7);p.add_argument('--seed',type=int,default=20261104);p.add_argument('--pool-df',type=float,default=2);p.add_argument('--bootstraps',type=int,default=999);p.add_argument('--rhos',nargs='+',type=float,default=[0,.3,.6]);p.add_argument('--series');run(p.parse_args())
