"""Identification, influence and conditional simulation checks for joint review fits."""
from pathlib import Path
import json,argparse
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import skew
from joint_borrowed_biomass import fit,build,BASE
O=Path('outputs/joint_borrowed_biomass_20261102')
def run(args):
    idx=pd.read_csv(O/'fit_index.csv');ident=[];profiles=[];cal=[];multistart=[];siteinfluence=[];countchecks=[];centering=[]
    samples=pd.read_csv(BASE/'uncertainty/size_observation_registry.csv');vis=pd.read_csv(BASE/'borrowing/visit_support.csv');vis=vis[vis.tolerance_days.eq(7)];counts=pd.read_csv(BASE/'density/biomass_roster_count_units.csv');counts['series_id']=counts.biomass_series_id
    offsets=pd.read_csv('outputs/all_program_station_refits_20261102/offsets.csv').set_index('series_id').offset
    for row in idx.itertuples():
        a=np.load(O/'fits'/f'{row.key}.npz');y,X,S,K,q,b,years=[a[k] for k in ['y','X','S','K','q','beta','years']];f=fit(y,X,S,K,q)
        active=np.array([np.any(k) for k in K]);ka=K[active];P=f['P'];info=np.array([[.5*np.trace(P@i@P@j) for j in ka] for i in ka]);diag=np.sqrt(np.diag(info));corr=info/diag[:,None]/diag[None,:]
        ident.append(dict(series_id=row.series_id,rho=row.rho,active_components=active.sum(),projected_kernel_rank=np.linalg.matrix_rank(info),information_condition=np.linalg.cond(corr),max_information_correlation=np.max(abs(corr-np.eye(len(corr)))),process_vs_sizeyear_correlation=max(abs(corr[-1,:-1])),fixed_design_rank=np.linalg.matrix_rank(X),fixed_design_columns=X.shape[1],min_covariance_eigenvalue=np.linalg.eigvalsh(f['C']).min()))
        if row.rho!=0:continue
        for start in [np.zeros(len(q)),np.full(len(q),1.)]:
            ff=fit(y,X,S,K,start);multistart.append(dict(series_id=row.series_id,start=float(start[0]),objective_difference=ff['obj']-f['obj'],slope_difference=ff['beta'][-1]-f['beta'][-1]))
        # Refit nuisance variances for a coarse profile; slope intervals remain bootstrap-based.
        for qproc in [0.,.25*q[-1],q[-1],4*q[-1]]:
            Sfix=S+qproc*K[-1];ff=fit(y,X,Sfix,K[:-1],q[:-1]);profiles.append(dict(series_id=row.series_id,process_variance_fixed=qproc,restricted_objective_difference=ff['obj']-f['obj'],slope=ff['beta'][-1],scope='Coarse sensitivity profile, not a profile confidence interval'))
        s=samples[samples.series_id.eq(row.series_id)].reset_index(drop=True);v=vis[(vis.series_id==row.series_id)&vis.year.isin(years)].reset_index(drop=True);cc=counts[counts.series_id==row.series_id]
        d=build(s,v,years,offsets[row.series_id],cc,np.random.default_rng(901),rho=0)
        bias=d['count_replicates'].mean(axis=1)-d['h'];yc=y.copy();yc[-len(years):]-=bias;fc=fit(yc,X,S,K,q)
        centering.append(dict(series_id=row.series_id,estimate=fc['beta'][-1],change=fc['beta'][-1]-f['beta'][-1],scope='Resampling-centering sensitivity, not a validated bias correction'))
        for j,yr in enumerate(years):countchecks.append(dict(series_id=row.series_id,year=yr,log_count_sampling_bias=d['count_replicates'][j].mean()-d['h'][j],sampling_sd=d['count_replicates'][j].std(ddof=1),sampling_skew=skew(d['count_replicates'][j])))
        for site in sorted(s[s.status=='reference'].site.unique()):
            # Remove size observations only; preserve count histories and donor-target random structure.
            keep=np.r_[s.site.ne(site).to_numpy(),np.ones(len(years),bool)]
            if np.linalg.matrix_rank(X[keep])!=X.shape[1]:continue
            ff=fit(y[keep],X[keep],S[np.ix_(keep,keep)],[k[np.ix_(keep,keep)] for k in K],q)
            siteinfluence.append(dict(series_id=row.series_id,omitted_size_site=site,removed_samples=int((~s.site.ne(site)).sum()),slope=ff['beta'][-1],change=ff['beta'][-1]-f['beta'][-1],scope='Remove this site size measurements only; retain count roster and structural size model'))
        # Finite simulation calibration at the fitted working model and at zero process variance.
        # This is a smoke test with wide Monte Carlo uncertainty, not nominal-coverage certification.
        for scenario,qtrue in [('fitted',q),('zero_process',np.r_[q[:-1],0.])]:
            C=S+np.einsum('k,kij->ij',qtrue,K)+np.eye(len(y))*1e-9;L=np.linalg.cholesky(C);rng=np.random.default_rng(552+int(row.Index)*101+(scenario=='zero_process'))
            for sim in range(args.simulations):
                rng=np.random.default_rng(552+int(row.Index)*100001+(scenario=='zero_process')*10001+sim)
                yy=X@b+L@rng.normal(size=len(y));ff=fit(yy,X,S,K,qtrue);se=np.sqrt(ff['bc'][-1,-1]);piv=[];LL=np.linalg.cholesky(ff['C'])
                for j in range(args.bootstraps):
                    yb=X@ff['beta']+LL@rng.normal(size=len(y));fb=fit(yb,X,S,K,ff['q']);piv.append((fb['beta'][-1]-ff['beta'][-1])/np.sqrt(fb['bc'][-1,-1]))
                lo=ff['beta'][-1]-np.quantile(piv,.975)*se;hi=ff['beta'][-1]-np.quantile(piv,.025)*se
                cal.append(dict(series_id=row.series_id,scenario=scenario,simulation=sim,true_slope=b[-1],estimate=ff['beta'][-1],low=lo,high=hi,covered=lo<=b[-1]<=hi,upper_boundary=ff['upper_boundary'],bootstraps=args.bootstraps))
        print(row.series_id,'diagnostics complete',flush=True)
    for name,rows in [('identification',ident),('process_profiles',profiles),('simulation_calibration',cal),('multiple_starts',multistart),('leave_one_size_site',siteinfluence),('count_sampling_checks',countchecks),('count_centering_sensitivity',centering)]:pd.DataFrame(rows).to_csv(O/f'{name}.csv',index=False)
    summary=pd.DataFrame(cal).groupby(['series_id','scenario']).agg(simulations=('covered','size'),covered=('covered','sum'),coverage=('covered','mean'),bias=('estimate',lambda x:0)).reset_index()
    # Wilson uncertainty intervals for estimated coverage (finite simulation count).
    n=summary.simulations;p=summary.coverage;z=1.96;mid=(p+z*z/(2*n))/(1+z*z/n);half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
    summary['coverage_mc_low']=mid-half;summary['coverage_mc_high']=mid+half
    bias=pd.DataFrame(cal).assign(error=lambda d:d.estimate-d.true_slope).groupby(['series_id','scenario']).error.mean()
    summary['bias']=[bias.loc[(r.series_id,r.scenario)] for r in summary.itertuples()]
    summary.to_csv(O/'simulation_summary.csv',index=False)
    (O/'diagnostic_checks.json').write_text(json.dumps(dict(simulations=args.simulations,inner_bootstraps=args.bootstraps,scope='Conditional first-order Gaussian model checks, not nonlinear coverage certification',execution_complete=True,nominal_coverage_certified=False),indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--simulations',type=int,default=20);p.add_argument('--bootstraps',type=int,default=199);run(p.parse_args())
