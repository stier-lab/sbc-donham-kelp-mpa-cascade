"""Paired conditional nonlinear calibration of the joint working likelihood.
Exact size aggregation and sample-mean/count resampling; fixed support and nuisance
sampling inputs. This is not a fitted full nonlinear biological likelihood.
"""
from pathlib import Path
import argparse,json,hashlib,copy
import numpy as np
import pandas as pd
from scipy.linalg import block_diag
from scipy.stats import norm
from joint_borrowed_biomass import fit,build,latent,aggregation,BASE
J=Path('outputs/joint_borrowed_biomass_20261102')


def root(C):
    e,Q=np.linalg.eigh((C+C.T)/2)
    if e.min() < -1e-6:raise ValueError('Non-PSD covariance')
    return Q*np.sqrt(np.maximum(e,0))


def exact_ratio(A,density,z,offset):
    """Column-wise exact arithmetic aggregation; zeros retained, no log(0) mass."""
    total=np.einsum('yri,id->yrd',A,density*np.exp(z))
    return np.log((total[:,0]+offset)/(total[:,1]+offset))


def density_replicates(v,counts,rng,n):
    units=sorted(counts.physical_unit_id.unique());uid={u:i for i,u in enumerate(units)};weights=rng.exponential(size=(len(units),n));dr=np.empty((len(v),n))
    for i,t in v.iterrows():
        cc=counts[(counts.site==t.site)&(counts.year==t.year)&(counts.zone==t.zone)&(counts.date==t.date)]
        vals=[]
        for _,g in cc.groupby('method'):
            w=weights[[uid[x] for x in g.physical_unit_id]];vals.append(g.density.to_numpy()@w/w.sum(axis=0))
        if not vals:raise ValueError('Missing count units')
        dr[i]=np.mean(vals,axis=0)
    return dr


def simulated_vectors(d,f,v,s,counts,offset,rng,n):
    """Common latent effects, process and Gaussian errors across paired branches."""
    ns=d['ns'];nt=len(v);ny=len(d['h']);W=d['W'];q=f['q'];alpha=f['beta'][:6];Z=d['X'][ns:,-2:]
    # Recover size-observation random covariance from the SAME kernel as joint fit.
    Gss=sum(a*k[:ns,:ns] for a,k in zip(q[:6],d['K'][:6]));Gts=sum(a*k for a,k in zip(q[:6],d['Kcross']));Gtt=sum(a*k for a,k in zip(q[:6],d['Ktt']))
    G=np.block([[Gss,Gts.T],[Gts,Gtt]])
    u=root(G)@rng.normal(size=(ns+nt,n));zs=d['Xs']@alpha[:,None]+u[:ns];zt=d['Xt']@alpha[:,None]+u[ns:]
    es=np.sqrt(d['known'])[:,None]*rng.normal(size=(ns,n));ec=root(d['count_cov'])@rng.normal(size=(ny,n));ep=root(q[-1]*d['K'][-1][ns:,ns:])@rng.normal(size=(ny,n))
    process=Z@f['beta'][-2:,None]+ep
    linear_h=d['h'][:,None]+W@(zt-d['logm'][:,None]);h=exact_ratio(d['A'],v.density.to_numpy()[:,None],zt,offset)
    c0=d['h']-W@d['logm']
    ylin=np.vstack([zs+es,process-W@zt+ec])
    ymass=np.vstack([zs+es,process-h+c0[:,None]+ec])
    # Sample individuals with the fitted arithmetic mean and pooled individual CV.
    # Reuse one sample mean per physical measurement; target borrowing adds no observations.
    ys_exact=np.empty((ns,n))
    for i,t in s.iterrows():
        n_ind=int(t['n']);assert n_ind==t['n'] and n_ind>0
        v_ind=np.log1p(t.sampling_cv2)
        relative=np.exp(-.5*v_ind+np.sqrt(v_ind)*rng.normal(size=(n_ind,n))).mean(axis=0)
        ys_exact[i]=zs[i]+np.log(relative)+.5*d['known'][i]
    dr=density_replicates(v,counts,rng,n)
    count_exact=exact_ratio(d['A'],dr,zt,offset)-h
    yfull=np.vstack([ys_exact,process-h+c0[:,None]+count_exact])
    return dict(gaussian=ylin,exact_mass=ymass,exact_mass_samples_counts=yfull),dict(curvature=h-linear_h,count_exact=count_exact,size_error=ys_exact-zs,latent_targets=zt,latent_samples=zs)


def reexpress(d,y,v,years,offset,zstar):
    """Same observations expressed around another mass vector; S remains fixed."""
    mass=np.exp(zstar);mass[v.density.eq(0)]=0
    hh,WW,_=aggregation(v,years,mass,offset);cstar=hh-WW@zstar;c0=d['h']-d['W']@d['logm'];ns=d['ns']
    dd=copy.copy(d);dd['y']=y.copy();dd['y'][ns:]+=cstar-c0;dd['W']=WW;dd['h']=hh;dd['logm']=zstar;dd['mass']=mass
    dd['X']=d['X'].copy();dd['X'][ns:,:6]=-WW@d['Xt']
    dd['K']=[np.block([[ko[:ns,:ns],-kc.T@WW.T],[-WW@kc,WW@kt@WW.T]]) for ko,kc,kt in zip(d['K'][:6],d['Kcross'],d['Ktt'])]+[d['K'][-1]]
    return dd


def iterate(d,v,years,offset,start,maxiter=20,tol=1e-5):
    current=d;f=start;trace=[];converged=False
    for j in range(maxiter):
        mu,_=latent(current,f)
        # Damping guards against large expansion jumps. This remains an approximation.
        zstar=.5*current['logm']+.5*mu
        nxt=reexpress(current,current['y'],v,years,offset,zstar);ff=fit(nxt['y'],nxt['X'],nxt['S'],nxt['K'],f['q'])
        step=abs(ff['beta'][-1]-f['beta'][-1]);massstep=float(np.max(abs(zstar-current['logm'])[v.density.gt(0)]))
        trace.append(dict(iteration=j+1,estimate=ff['beta'][-1],slope_step=step,max_logmass_step=massstep))
        current=nxt;f=ff
        if step<tol and massstep<1e-3:converged=True;break
    return f,trace,converged


def run(args):
    if args.replicates<100:raise ValueError('At least100 paired replicates required')
    O=Path(args.output);O.mkdir(parents=True,exist_ok=True)
    for file in ['completed.json','report_verified.json']:
        if (O/file).exists():(O/file).unlink()
    sources=[BASE/'uncertainty/size_observation_registry.csv',BASE/'borrowing/visit_support.csv',BASE/'density/biomass_roster_count_units.csv',Path('outputs/all_program_station_refits_20261102/offsets.csv'),J/'fit_index.csv',J/'joint_trends.csv',J/'bootstrap_fits.csv',Path(__file__),Path('code/tools/joint_borrowed_biomass.py')]
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    samples=pd.read_csv(sources[0]);vis=pd.read_csv(sources[1]);vis=vis[vis.tolerance_days==7];counts=pd.read_csv(sources[2]);counts['series_id']=counts.biomass_series_id
    offsets=pd.read_csv(sources[3]).set_index('series_id').offset;idx=pd.read_csv(sources[4]);idx=idx[idx.rho.isin(args.rhos)];results=[];trace=[];obs=[];decomp=[];support=[]
    for row in idx.itertuples():
        p=J/'fits'/f'{row.key}.npz';hashes[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();a=np.load(p);years=a['years'];s=samples[samples.series_id==row.series_id].reset_index(drop=True);v=vis[(vis.series_id==row.series_id)&vis.year.isin(years)].reset_index(drop=True);cc=counts[counts.series_id==row.series_id];off=offsets[row.series_id]
        si=int(row.key.split('_')[0][1:]);d=build(s,v,years,off,cc,np.random.default_rng(20261104+si*1000),row.rho);f=fit(d['y'],d['X'],d['S'],d['K'],a['q']);assert np.allclose(f['beta'],a['beta'],atol=1e-5)
        support.append(dict(series_id=row.series_id,rho=row.rho,years=';'.join(map(str,years)),n_years=len(years),size_samples=len(s),visits=len(v)))
        ff,tr,conv=iterate(d,v,years,off,f)
        obs.append(dict(series_id=row.series_id,rho=row.rho,original=f['beta'][-1],iterated=ff['beta'][-1],change=ff['beta'][-1]-f['beta'][-1],converged=conv,iterations=len(tr),interval_claim='No new interval; iterated expansion diagnostic only'))
        for rr in tr:trace.append(dict(series_id=row.series_id,rho=row.rho,**rr))
        rng=np.random.default_rng(args.seed+si*1000+int(row.rho*100));vectors,detail=simulated_vectors(d,f,v,s,cc,off,rng,args.replicates)
        for j in range(args.replicates):
            for kind,Y in vectors.items():
                fitj=fit(Y[:,j],d['X'],d['S'],d['K'],f['q']);se=float(np.sqrt(fitj['bc'][-1,-1]));estimate=fitj['beta'][-1]
                results.append(dict(series_id=row.series_id,rho=row.rho,replicate=j,generator=kind,true_slope=f['beta'][-1],estimate=estimate,se=se,error=estimate-f['beta'][-1],pivot=(estimate-f['beta'][-1])/se,optimizer_attempts=fitj['optimizer_attempts'],upper_boundary=fitj['upper_boundary']))
        x=years-years.mean();slope_error=x@detail['curvature']/(x@x)
        for j,yr in enumerate(years):decomp.append(dict(series_id=row.series_id,rho=row.rho,year=yr,mean_mass_curvature=detail['curvature'][j].mean(),q95_abs_mass_curvature=np.quantile(abs(detail['curvature'][j]),.95),mean_exact_count_error=detail['count_exact'][j].mean(),sd_exact_count_error=detail['count_exact'][j].std(ddof=1)))
        print(row.series_id,row.rho,'iterated change',obs[-1]['change'],'paired fits',3*args.replicates,flush=True)
        pd.DataFrame(results).to_csv(O/'paired_fits.csv',index=False)
    for name,rows in [('observed_iteration',obs),('iteration_trace',trace),('nonlinear_components',decomp),('support',support)]:pd.DataFrame(rows).to_csv(O/f'{name}.csv',index=False)
    assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items())
    (O/'input_hashes.json').write_text(json.dumps(hashes,indent=2));(O/'policy.json').write_text(json.dumps(vars(args),indent=2));(O/'completed.json').write_text(json.dumps(dict(review_only=True,paired_replicates=args.replicates,model_combinations=len(idx),fit_count=len(results),full_nonlinear_likelihood_fitted=False,total_uncertainty_complete=False),indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',default='outputs/nonlinear_biomass_calibration_20261102');p.add_argument('--replicates',type=int,default=999);p.add_argument('--rhos',nargs='+',type=float,default=[0,.6]);p.add_argument('--seed',type=int,default=20261106);run(p.parse_args())
