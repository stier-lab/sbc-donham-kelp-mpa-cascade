"""Report paired nonlinear calibration separately from primary inference."""
from pathlib import Path
import json,hashlib,html
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path('outputs/nonlinear_biomass_calibration_20261102');J=Path('outputs/joint_borrowed_biomass_20261102')
def wilson(k,n):
    p=k/n;z=1.96;den=1+z*z/n;mid=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den;return mid-half,mid+half

def md(d):return '| '+' | '.join(d.columns)+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in d.itertuples(index=False,name=None))
def run():
    assert (O/'completed.json').exists();assert json.loads((O/'checks.json').read_text())['success']
    fits=pd.read_csv(O/'paired_fits.csv');old=pd.read_csv(J/'joint_trends.csv').set_index(['series_id','rho']);boots=pd.read_csv(J/'bootstrap_fits.csv');rows=[];paired=[]
    for (sid,rho),g in fits.groupby(['series_id','rho']):
        b=old.loc[(sid,rho)];bg=boots[(boots.series_id==sid)&(boots.rho==rho)];oldq=np.quantile((bg.slope-b.estimate)/bg.se,[.025,.975]);ga=g[g.generator=='gaussian'].set_index('replicate')
        for kind,z in g.groupby('generator'):
            z=z.sort_values('replicate');n=len(z);q=np.quantile(z['pivot'],[.025,.975]);split=n//2;cal=z.iloc[:split];test=z.iloc[split:];sq=np.quantile(cal['pivot'],[.025,.975]);oldcovered=((z['pivot']>=oldq[0])&(z['pivot']<=oldq[1])).sum();newcovered=((test['pivot']>=sq[0])&(test['pivot']<=sq[1])).sum();lo,hi=wilson(oldcovered,n);cl,ch=wilson(newcovered,len(test))
            rows.append(dict(series_id=sid,rho=rho,generator=kind,replicates=n,original_estimate=b.estimate,original_low=b.low,original_high=b.high,mean_slope_error=z.error.mean(),median_slope_error=z.error.median(),calibrated_low=b.estimate-q[1]*b.se,calibrated_high=b.estimate-q[0]*b.se,old_pivot_coverage=oldcovered/n,old_coverage_mc_low=lo,old_coverage_mc_high=hi,split_calibration_n=len(cal),split_validation_n=len(test),split_pivot_coverage=newcovered/len(test),split_coverage_mc_low=cl,split_coverage_mc_high=ch,upper_bound_fits=int(z.upper_boundary.sum()),optimizer_retries=int(z.optimizer_attempts.gt(1).sum()),scope='Conditional pivot calibration at fitted nuisance parameters; not full nonlinear likelihood or nominal-coverage certification'))
            delta=z.set_index('replicate').error-ga.error
            paired.append(dict(series_id=sid,rho=rho,generator=kind,mean_incremental_bias=delta.mean(),mc_se=delta.std(ddof=1)/np.sqrt(n),mc_low=delta.mean()-1.96*delta.std(ddof=1)/np.sqrt(n),mc_high=delta.mean()+1.96*delta.std(ddof=1)/np.sqrt(n)))
    cal=pd.DataFrame(rows);cal.to_csv(O/'calibration_summary.csv',index=False);pd.DataFrame(paired).to_csv(O/'paired_bias_changes.csv',index=False)
    it=pd.read_csv(O/'observed_iteration.csv');base=it[it.rho==0].merge(cal[(cal.rho==0)&(cal.generator=='exact_mass_samples_counts')],on=['series_id','rho'],validate='one_to_one');base['interval_width_ratio']=(base.calibrated_high-base.calibrated_low)/(base.original_high-base.original_low);extra_bias=pd.DataFrame(paired);extra_bias=extra_bias[(extra_bias.rho==0)&(extra_bias.generator=='exact_mass_samples_counts')][['series_id','mean_incremental_bias','mc_low','mc_high']]
    base=base.merge(extra_bias,on='series_id',validate='one_to_one');base['incremental_bias_adjusted_diagnostic']=base.original-base.mean_incremental_bias
    base.to_csv(O/'same_year_comparison.csv',index=False)
    plots=[];names={'M. franciscanus':'Red urchin','S. purpuratus':'Purple urchin','P. interruptus':'Lobster'}
    for r in base.itertuples():
        site,taxon,*_=r.series_id.split('|');plots.append(dict(Site=site.replace(' SMCA',''),Taxon=names[taxon],Initial=f'{r.original:.3f}',Iterated=f'{r.iterated:.3f}',Change=f'{r.change:+.3f}',Initial_interval=f'[{r.original_low:.3f}, {r.original_high:.3f}]',Calibrated_interval=f'[{r.calibrated_low:.3f}, {r.calibrated_high:.3f}]'))
    display=pd.DataFrame(plots)
    fig,axes=plt.subplots(1,2,figsize=(9,4.4),sharey=True)
    for ax,site,letter in zip(axes,['Abalone Cove SMCA','Point Vicente SMCA'],['a','b']):
        for i,taxon in enumerate(['M. franciscanus','S. purpuratus','P. interruptus']):
            r=base[base.series_id.eq(f'{site}|{taxon}|PISCO|Bio')].iloc[0]
            for dy,lo,hi,color,marker in [(-.12,r.original_low,r.original_high,'#777777','o'),(.12,r.calibrated_low,r.calibrated_high,'#0072B2','s')]:
                ax.plot([lo,hi],[i+dy,i+dy],color=color,lw=1.5);ax.scatter(r.original,i+dy,color=color,marker=marker,s=30,zorder=3)
        ax.axvline(0,color='#aaa',lw=.8);ax.set_xlabel('Annual change in log biomass ratio');ax.set_yticks(range(3),['Red urchin','Purple urchin','Lobster']);ax.set_ylim(-.55,2.55);ax.text(.02,.98,letter,transform=ax.transAxes,va='top');ax.text(.5,.98,site.replace(' SMCA',''),transform=ax.transAxes,ha='center',va='top');ax.spines[['top','right']].set_visible(False)
    lim=abs(base[['original_low','original_high','calibrated_low','calibrated_high']].to_numpy()).max()*1.08
    for ax in axes:ax.set_xlim(-lim,lim)
    fig.tight_layout();fig.savefig(O/'nonlinear_intervals.png',dpi=180);fig.savefig(O/'nonlinear_intervals.pdf');plt.close(fig)
    caption='Protected/reference biomass recovery slopes at Abalone Cove (a) and Point Vicente (b), on unchanged years. Gray circles show the existing joint-model interval; blue squares show conditional pivot intervals calibrated with exact arithmetic biomass aggregation, simulated individual-mass means and empirical count-unit resampling. Both markers show the same original fitted slope to isolate interval changes. Iterated point estimates are reported separately in the table. Calendar correlation is zero here; correlation 0.6 is retained in the linked results. These are diagnostic conditional intervals, not a full nonlinear fit or total uncertainty.'
    (O/'figure_caption.md').write_text(caption+'\n')
    receipt=json.loads((O/'kfm_source_receipt.json').read_text())
    text='''# Nonlinear biomass calibration

The same six PISCO histories and 45 post-protection comparison-years are retained. Researcher exclusions and unsupported biomass cells are unchanged. This review tests the approximation used to combine sizes and densities; it does not update the paper or primary analysis.

## What changes

'''+md(display)+'''

Initial and iterated values are annual slopes of the protected/reference biomass log ratio. Iteration updates the expansion around fitted sizes while preserving observations and sampling covariance; it is an approximation check, not an exact nonlinear likelihood. Calibrated intervals apply the existing fitting procedure to simulated nonlinear observations. They are not intervals from the iterated fit.

![Conditional interval comparison](nonlinear_intervals.png)

'''+caption+'''

## What was tested

Each comparison is simulated at calendar correlations 0 and 0.6, keeping the observed years and fitted nuisance parameters. Three paired branches share the same latent size effects and annual process draws:

- Gaussian working observations reproduce the current model and provide a control.
- Exact mass aggregation replaces the first-order size approximation, with the same Gaussian measurement errors. The paired difference isolates curvature.
- Exact mass aggregation plus sample/count resampling also simulates individual masses from a lognormal distribution with the fitted arithmetic mean and pooled dispersion, averages the observed number of individuals, and resamples physical count units. Sizes enter once, shared effects remain shared, and biomass is multiplied and aggregated exactly. Known zero-count visits remain zero.

All branches refit the current joint working model. Half of each simulation set calibrates its pivot quantiles and the other half checks those quantiles on independent draws. Existing Gaussian pivot quantiles are also checked against the nonlinear draws. Monte Carlo bounds are reported. This conditional check is not a nested validation of the complete interval-fitting procedure over unknown parameter values.

The nonlinear count contribution is an observation-level residual around the existing count design. It does not regenerate new sites, missing surveys, eligibility decisions or a complete raw-data count process. Expansion weights and estimated sampling variances remain fixed during calibration; the separate observed-data iteration tests expansion sensitivity. Thus this extends nonlinear calibration without claiming that every stage of the original pipeline has been simulated.

## Scott's source delivery

Scott Gabara supplied **KFM Nat Hab Red Purp 2012-2022.xlsx** on November 2. The unchanged workbook is preserved locally with SHA256 `'''+receipt['sha256']+'''`. It has 55,142 nonempty rows, 33 site names, exact millimeter sizes and frequency counts. Actual years are 2012–2019 and 2021–2022; there are no 2020 records. Distinct Southeast Reef/SE Sea Lion Rookery and Miracle Mile/Wyckoff Ledge labels are present. Some bin labels were converted to Excel dates; the exact-size column must be used.

Receipt is now resolved. Matching these records to the held count histories and original site IDs remains necessary; no KFM biomass hold is cleared by receipt alone. The workbook contains red/purple urchins, not lobster or sheephead. No reply was sent. [Receipt record](kfm_source_receipt.json).

## Interpretation limits and reproduction

Allometry, estimated sampling dispersion, offsets, the observed count-unit population and support are conditional. Lognormal individual masses are a working assumption. Detection, count–size dependence, unobserved animals at zero-count visits, protected-site transportability, cross-history covariance and broader missing-data mechanisms remain untested. This is calibration of the existing estimator, not a newly fitted full nonlinear likelihood or total uncertainty.

Run `python3 code/tools/run_nonlinear_biomass_review.py`. [Paired bias changes](paired_bias_changes.csv), [calibration coverage](calibration_summary.csv), [iteration trace](iteration_trace.csv), [nonlinear components](nonlinear_components.csv), [checks](checks.json), and [consolidated guide](../dataset_inventory_20261101/README.md) retain details.
'''
    maxchange=base.change.abs().max();nconv=int(base.converged.sum());full=cal[cal.generator=='exact_mass_samples_counts']
    text+='\n## Verified summary\n\nObserved expansion iteration converged in '+str(nconv)+' of six rho=0 histories. The largest slope change was '+f'{maxchange:.4f}'+'. The exact-generation calibrated interval width relative to the previous interval ranged from '+f'{base.interval_width_ratio.min():.2f} to {base.interval_width_ratio.max():.2f}'+'. See the table for any sign or interval changes; near-zero estimates should not be treated as directionally established.\n\n'
    text+='The largest simulated incremental bias relative to the paired Gaussian control was '+f'{base.mean_incremental_bias.abs().max():.4f}'+' per year, for Point Vicente purple urchins. This is evidence about estimator bias under the specified generator, not a directly observed biological change. The linked comparison includes an incremental-bias-adjusted diagnostic; no correction has been adopted for the paper.\n\n'
    text+=f"All {len(fits):,} paired refits completed; {int(fits.optimizer_attempts.gt(1).sum())} required optimizer retries. {int(fits.upper_boundary.sum())} simulated fits reached the variance bound and remain flagged in the results; no failed fit was substituted or silently dropped.\n\n"
    text+='The independent split check for the combined nonlinear branch covered '+f'{full.split_pivot_coverage.min():.1%}–{full.split_pivot_coverage.max():.1%}'+' across histories/correlation scenarios. The lowest split coverage was Point Vicente red urchins at correlation zero (93%). Finite calibration quantiles and validation draws both add Monte Carlo uncertainty; do not label all settings validated at nominal 95%. These results concern conditional pivot calibration and do not certify the full data-generating process.\n'
    (O/'README.md').write_text(text)
    body='<h1>Nonlinear biomass calibration</h1><p>Same six histories and 45 comparison-years. Review only; paper unchanged.</p>'+display.to_html(index=False)+'<img src="nonlinear_intervals.png" style="max-width:100%"><p>'+html.escape(caption)+'</p><p>Scott’s urchin workbook has arrived; site reconciliation remains open.</p><p><a href="README.md">Methods, findings and limits</a> · <a href="calibration_summary.csv">Calibration results</a> · <a href="../dataset_inventory_20261101/index.html">Consolidated guide</a></p>'
    (O/'index.html').write_text('<!doctype html><meta charset="utf-8"><style>body{font:16px system-ui;max-width:1200px;margin:28px auto;padding:0 18px}td,th{padding:8px;text-align:left}table{border-collapse:collapse}tr{border-bottom:1px solid #ccc}</style>'+body)
    protected=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text());assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in protected.items())
    files=['same_year_comparison.csv','calibration_summary.csv','paired_bias_changes.csv','support.csv','checks.json','kfm_source_receipt.json','observed_iteration.csv']
    (O/'report_verified.json').write_text(json.dumps(dict(review_only=True,primary_updated=False,total_uncertainty_complete=False,output_hashes={n:hashlib.sha256((O/n).read_bytes()).hexdigest() for n in files}),indent=2))
    older=J/'README.md';note='> Earlier review stage. See the [nonlinear calibration and verified Scott source delivery](../nonlinear_biomass_calibration_20261102/README.md) for the current extension and receipt status.\n\n'
    if not older.read_text().startswith('> Earlier review stage.'):older.write_text(note+older.read_text())
    print(display.to_string(index=False));print('Protected files unchanged:',len(protected))
if __name__=='__main__':run()
