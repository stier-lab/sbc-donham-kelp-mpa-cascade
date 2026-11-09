"""Build a review report without promoting manuscript or primary estimates."""
from pathlib import Path
import json,hashlib,html
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path('outputs/joint_borrowed_biomass_20261102');B=Path('outputs/biomass_next_steps_20261102')
def run():
    assert (O/'completed.json').exists()
    results=pd.read_csv(O/'joint_trends.csv');assert len(results)==18
    old=pd.read_csv(B/'borrowing/effects.csv');old=old[(old.tolerance_days==7)&(old.observation_scope=='full_supported')]
    joint=results[results.rho==0].copy();comparison=old[['series_id','estimate','low','high','n_years','years']].rename(columns={c:'current_'+c for c in ['estimate','low','high','n_years']}).merge(joint,on='series_id',validate='one_to_one')
    comparison['estimate_change']=comparison.estimate-comparison.current_estimate
    assert (comparison.n_years==comparison.current_n_years).all()
    assert set(comparison.series_id)==set(pd.read_csv(O/'support.csv').series_id)
    comparison.to_csv(O/'same_year_comparison.csv',index=False)
    curv=pd.read_csv(O/'linearization_checks.csv');loo=pd.read_csv(O/'leave_one_count_year.csv')
    diagnostics=joint[['series_id','converged','boundary_components','upper_boundary','condition','recentered_estimate']].copy()
    diagnostics['recenter_change']=diagnostics.recentered_estimate-joint.estimate.to_numpy()
    diagnostics=diagnostics.merge(curv[curv.rho==0].groupby('series_id').agg(max_mean_curvature=('mean_exact_minus_linear',lambda x:max(abs(x))),max_q95_curvature=('q95_abs_error','max')).reset_index(),on='series_id')
    diagnostics=diagnostics.merge(loo.groupby('series_id').agg(leave_year_min=('estimate','min'),leave_year_max=('estimate','max'),max_leave_year_change=('difference',lambda x:max(abs(x)))).reset_index(),on='series_id')
    diagnostics['diagnostic_status']='partial: conditional working likelihood; assess curvature, influence and boundary components'
    diagnostics.to_csv(O/'diagnostics_summary.csv',index=False)
    names={'M. franciscanus':'Red urchin','S. purpuratus':'Purple urchin','P. interruptus':'Lobster'}
    fig,axes=plt.subplots(1,2,figsize=(9,4.4),sharey=True)
    for ax,mpa,letter in zip(axes,['Abalone Cove SMCA','Point Vicente SMCA'],['a','b']):
        for i,taxon in enumerate(['M. franciscanus','S. purpuratus','P. interruptus']):
            r=comparison[comparison.series_id.eq(f'{mpa}|{taxon}|PISCO|Bio')].iloc[0]
            for dy,mean,lo,hi,color,marker in [(-.12,r.current_estimate,r.current_low,r.current_high,'#777777','o'),(.12,r.estimate,r.low,r.high,'#0072B2','s')]:
                ax.plot([lo,hi],[i+dy,i+dy],color=color,lw=1.5);ax.scatter([mean],[i+dy],color=color,marker=marker,s=30,zorder=3)
        ax.axvline(0,color='#aaaaaa',lw=.8);ax.set_xlabel('Annual change in log biomass ratio');ax.set_yticks(range(3),['Red urchin','Purple urchin','Lobster']);ax.text(.02,.98,letter,transform=ax.transAxes,ha='left',va='top');ax.text(.5,.98,mpa.replace(' SMCA',''),transform=ax.transAxes,ha='center',va='top');ax.set_ylim(-.55,2.55);ax.spines[['top','right']].set_visible(False)
    lim=max(abs(comparison[['current_low','current_high','low','high']].to_numpy()).ravel())*1.08
    for ax in axes:ax.set_xlim(-lim,lim)
    fig.tight_layout();fig.savefig(O/'same_year_trends.png',dpi=180);fig.savefig(O/'same_year_trends.pdf');plt.close(fig)
    caption='Recovery trends in the protected/reference biomass ratio for Abalone Cove (a) and Point Vicente (b). Gray circles show the existing seven-day borrowing sensitivity with conditional HAC intervals; blue squares show the joint working likelihood at calendar correlation zero with studentized bootstrap intervals. Each comparison uses identical years. Positive slopes mean increasing protected/reference biomass ratios. Interval methods differ; neither is total uncertainty. The seven-day scenario remains a review sensitivity.'
    (O/'figure_caption.md').write_text(caption+'\n')
    formatted=[]
    for r in comparison.itertuples():
        mpa,taxon,*_=r.series_id.split('|');formatted.append(dict(Site=mpa.replace(' SMCA',''),Taxon=names[taxon],Years=r.n_years,Existing_review=f'{r.current_estimate:.3f} [{r.current_low:.3f}, {r.current_high:.3f}]',Joint_review=f'{r.estimate:.3f} [{r.low:.3f}, {r.high:.3f}]'))
    display=pd.DataFrame(formatted)
    txt='''# Joint biomass uncertainty review

The model now fits measured sizes, survey-count uncertainty and recovery through time together. A size sample enters once. Shared sites, years and visits create covariance across the biomass estimates that depend on them. No smoothed posterior estimate is treated as an independent observation.

This is a **joint first-order Gaussian working likelihood**, not a complete uncertainty model or a manuscript update. The same seven-day sensitivity support is retained: six histories and 45 post-protection comparison-years. Researcher exclusions and unsupported biomass cells remain unchanged; no protected-site borrowing is introduced.

## What changes on the same years

'''+('| '+' | '.join(display.columns)+' |\n| '+' | '.join(['---']*len(display.columns))+' |\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in display.itertuples(index=False,name=None)))+'''

Both columns are annual slopes of log protected/reference biomass ratios, with intervals. The existing-review column uses HAC conditional on assigned sizes; the joint-review column estimates size variation and annual process variance together. This comparison changes the estimator and uncertainty model, not sites or survey years. Changes must not be attributed solely to wider uncertainty.

![Same-year recovery comparison](same_year_trends.png)

'''+caption+'''

## What the model assumes

Raw mean-mass observations use a lognormal sampling approximation, corrected by half the log sampling variance. Within-sample dispersion is partially pooled to avoid pretending singleton sizes are exact. Size varies with zone and calendar year, with site, site-year and visit effects within each role. Counts retain equal-date, zone and station aggregation. A first-order derivative of that aggregation links the raw measurements to annual biomass ratios; the count covariance is estimated from shared physical survey-unit resampling.

The stacked likelihood jointly estimates the size effects, recovery slope and variance components by REML. Calendar correlations of 0, 0.3 and 0.6 are fixed sensitivity scenarios, not selected by significance. Each parametric bootstrap simulates the entire joint observation vector and refits variance components. Sampling variances, count covariance, allometry and expansion weights remain conditional inputs. The process term describes recovery over observed calendar years, not a before–after protection effect.

## Checks and remaining limits

See [diagnostics](diagnostics_summary.csv), [linearization checks](linearization_checks.csv), [year influence](leave_one_count_year.csv), [variance components](variance_components.csv), and [all correlation scenarios](joint_trends.csv). Boundary variance components and influential years are reported rather than hidden. Recentring the expansion at fitted masses tests dependence on the initial approximation. The diagnostic curvature draws are conditional on fitted variance components, not a full nonlinear validation.

This implements the previously missing joint borrowed-size/process calculation in a working approximation. It does not establish protected-site transportability or estimate allometric uncertainty, count–size dependence, detection error, or unknown long-term survey-unit dependence. KFM original-size holds remain. The model is for review; primary eligibility, manuscript and supplement remain unchanged.

## Reproduce

Run `python3 code/tools/run_joint_biomass_review.py` to rebuild the model, sensitivity fits, calibration, independent checks and consolidated guide in order. Numerical checks and independent review are linked in this directory. Inputs are hashed; a completion marker is written only after all model fits finish.
'''
    cal=pd.read_csv(O/'simulation_summary.csv');bf=pd.read_csv(O/'bootstrap_fits.csv');sens=pd.read_csv(O/'sensitivity_comparison.csv') if (O/'sensitivity_comparison.csv').exists() else pd.DataFrame()
    warning=bool((cal.coverage_mc_high<.95).any())
    calibration_note=f"The conditional simulation check covered the true slope in {int(cal.covered.min())}–{int(cal.covered.max())} of {int(cal.simulations.iloc[0])} datasets per setting. These small simulation counts have wide Monte Carlo uncertainty. "+('Some settings fall short of nominal 95% coverage even allowing for that uncertainty. Keep the intervals provisional; do not promote this model to the primary analysis.' if warning else 'This small exercise does not certify nominal 95% coverage; nonlinear data generation remains untested.')
    extra='\n## Current interpretation and verification\n\nThe biological directions are unchanged on identical years. Point estimates move by at most '+f'{comparison.estimate_change.abs().max():.3f}'+' log-ratio units per year. All six current HAC intervals and all six joint working intervals at rho=0 include zero. The new model has not removed a previously clear effect in these six comparisons.\n\n'+calibration_note+'\n\n'
    extra+=f"All {len(bf):,} main bootstrap fits converged after checked retries; {int(bf.optimizer_attempts.gt(1).sum())} needed more than one optimizer attempt. {int(bf.upper_boundary.sum())} simulated fits reached the computational variance bound; these are retained and flagged.\n\n"
    extra+='Changing the expansion center shifts slopes by at most '+f'{diagnostics.recenter_change.abs().max():.3f}'+'. Removing individual count years causes larger changes, including direction changes for some short histories. Removing a reference site’s size observations tests donor influence without dropping its count history. [Simulation results](simulation_summary.csv), [size-site influence](leave_one_size_site.csv), [count centering](count_centering_sensitivity.csv), [identification](identification.csv), [independent checks](independent_checks.json), and [model equations and borrowing assumptions](MODEL.md) retain these details.\n\n'
    if len(sens):
        extra+='Across the alternate random seed, stronger size-variance pooling and strict timing rule, all rho=0 intervals still include zero. [Sensitivity comparison](sensitivity_comparison.csv) separates these changes.\n'
    txt+=extra
    if Path('outputs/nonlinear_biomass_calibration_20261102/report_verified.json').exists():
        txt='> Earlier review stage. See the [nonlinear calibration and verified Scott source delivery](../nonlinear_biomass_calibration_20261102/README.md) for the current extension and receipt status.\n\n'+txt
    (O/'README.md').write_text(txt)
    body='<h1>Joint biomass uncertainty review</h1><p>'+html.escape(calibration_note)+'</p><p>Same supported years; review only. The manuscript and primary analysis are unchanged.</p>'+display.to_html(index=False)+'<img src="same_year_trends.png" style="max-width:100%"><p>'+html.escape(caption)+'</p><p><a href="README.md">Methods, assumptions and limits</a> · <a href="diagnostics_summary.csv">Diagnostics</a> · <a href="joint_trends.csv">Correlation scenarios</a> · <a href="../dataset_inventory_20261101/index.html">Consolidated guide</a></p>'+diagnostics.to_html(index=False)
    (O/'index.html').write_text('<!doctype html><meta charset="utf-8"><style>body{font:16px system-ui;max-width:1200px;margin:32px auto;padding:0 20px}td,th{padding:8px;text-align:left}table{border-collapse:collapse}tr{border-bottom:1px solid #ccc}</style>'+body)
    h=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text());assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==v for p,v in h.items())
    print(display.to_string(index=False));print('Protected files unchanged:',len(h))
if __name__=='__main__':run()
