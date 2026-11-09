from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path('outputs/pisco_biomass_borrowing_20261102');e=pd.read_csv(O/'biomass_effects.csv');c=pd.read_csv(O/'biomass_coverage.csv');v=pd.read_csv(O/'visit_biomass_support.csv');cv=pd.read_csv(O/'heldout_interval_validation.csv');cv=cv[cv.predicted.notna()].copy();cv['interval_contains_observed']=(cv.observed>=cv.low)&(cv.observed<=cv.high)
validation=cv.groupby(['Taxa','status','kind'],as_index=False).agg(targets=('site','size'),heldout_sites=('site','nunique'),observed_coverage=('interval_contains_observed','mean'),median_calibration_siteyears=('cal_clusters','median'));validation.to_csv(O/'interval_validation_summary.csv',index=False)
summary=c.merge(e[['series_id','estimate','low','high','prediction_draw_low','prediction_draw_high']],on='series_id',how='left');summary.to_csv(O/'summary.csv',index=False)
fig,axs=plt.subplots(1,2,figsize=(10,4.3),sharex=True,sharey=True)
for ax,lo,hi,letter in zip(axs,['low','prediction_draw_low'],['high','prediction_draw_high'],'AB'):
 y=np.arange(len(e));ax.hlines(y,e[lo],e[hi],color='#177c80',lw=2);ax.scatter(e.estimate,y,color='#222',s=24);ax.axvline(0,color='.6',lw=.8);ax.set_yticks(y,[sid.split('|')[0].replace(' SMCA','')+' · '+{'M. franciscanus':'red urchin','S. purpuratus':'purple urchin','P. interruptus':'lobster'}[sid.split('|')[1]] for sid in e.series_id]);ax.text(.03,.98,letter,transform=ax.transAxes,va='top');ax.set_xlabel('Change in log biomass ratio / year');ax.spines[['top','right']].set_visible(False)
axs[0].set_ylim(len(e)-.5,-.8);fig.tight_layout();fig.savefig(O/'biomass_uncertainty_components.png',dpi=180);fig.savefig(O/'biomass_uncertainty_components.pdf');plt.close(fig)
text='''# Count-led review and borrowing uncertainty

The review pipeline now retains all-month density histories and provides a separate biomass sensitivity with empirical prediction uncertainty. Scope: the six affected PISCO animal comparisons at Abalone Cove and Point Vicente. Existing station/year membership, researcher exclusions and shared INNER/MID zones remain unchanged. The manuscript, canonical primary outputs and other programs have not been updated.

## What is implemented

`selected_density_annual.csv` carries the all-month density series forward unchanged:63post-protection comparison-years. `review_analysis_policy.json` records the response selection and uncertainty limits. Biomass is explicitly a conditional review sensitivity.

Each positive count visit receives one of three statuses: sizes measured at that site/zone/date; borrowed sizes with an empirical prediction-error distribution; or unsupported biomass. Counts of zero remain zero. Missing measurements never become zero, and an unsized positive visit cannot disappear silently from an annual reference average.

The preferred borrowing source is another visit in the same site/zone/year if its temporal-transfer uncertainty can be estimated from other sites. Otherwise, use the nearest dated measurements from other retained sites with the same protection role, zone, species and year. Donor sites receive equal weights; animal multiplicity is preserved within each site sample. No cross-protection, cross-zone or cross-year borrowing is introduced. A predicted date gap beyond the calibration data's observed support is left unsupported. At least two calibration site-year groups are needed to define a variable empirical distribution; this is a minimum computational guard, not a claim of adequate ecological replication.

## Error calibration without using the held-out station

For each target station, remove that station from the calibration data and from every calibration donor pool. Predict measurements at the remaining stations using other eligible measurements. Their log observed/predicted ratios form the empirical error distribution. This nested exclusion prevents the target station from helping predict the samples used to calibrate its own uncertainty. Calibration site-years receive equal weight, rather than letting large numbers of measured animals become independent validation replicates. The point conversion uses the empirical median correction.

For measured samples,1,000frequency-aware Bayesian size draws retain original measurement multiplicity and shared sample identities. For borrowed samples,1,000prediction draws apply the empirical transfer errors. One error multiplier is shared by borrowed visits of the same species/site/year and borrowing type, so those visits are not treated as independent errors. This is an explicit within-year dependence assumption; it does not capture all persistence across years or dependence among different target stations using a common donor.

We do not add a second donor-size bootstrap to borrowed predictions because the empirical residuals already contain donor and target sampling variation. Residuals predict measured sample means, not an independently observed true population mean. Allometric coefficients, counts, site membership and offsets remain fixed.

## Coverage and results

There are624positive visits with same-date measurements,131cross-site borrowed visits,8borrowed visits from another date at the same site,38unsupported positive visits and285zero-count visits. Borrowed values are labeled in `visit_biomass_support.csv`, including donors, date gaps and calibration sample counts.

Complete paired biomass coverage rises from25years under the preceding local-site/zone analysis to44of63comparison-years. The earlier25included some other-date local samples without this new transfer calibration, so the support rules are not identical. New coverage is: Abalone redurchin7/11,purpleurchin7/11,lobster10/11; Point Vicente redurchin5/10,purpleurchin6/10,lobster9/10. Remaining years are withheld from biomass only. The density series retains every one of the63years.

All six biomass slopes can now be estimated, but their conditional temporal95%intervals span zero. This does not mean that density patterns vanished or that borrowing itself caused the change: the biomass estimates use different supported years, and matched-year original estimates are provided. The plot separates temporal-model uncertainty (A) from local-size/borrowing draw variation (B); dots are conditional point estimates. PanelBsegments are central95%draw ranges, not confidence intervals for protection effects. The two components are not mechanically added, because they can overlap. Neither is total biomass uncertainty, and no causal before–after claim is made from these recovery slopes.

## Held-out interval check

For cross-site borrowing,420reference-site targets support honest outer-site interval checks. Nominal95%prediction ranges contain93.4%of redurchin sample means,93.6%of lobster sample means and93.6%of purpleurchin sample means. These are descriptive coverage estimates with dependent site/year samples, not independent-binomial guarantees. Abalone purpleurchin reference coverage is lower, about88%, so species-wide averages do not certify every location.

Only eight same-site other-date validations support the temporal-transfer intervals; all are reference samples and calibration is sparse. No protected-site borrowing is certified by these checks. Protected-site donor/calibration gaps remain unsupported, including Abalone's single protected station. The empirical ranges may be inaccurate for unsampled conditions, and individual measured means can be noisy when few animals were sized. Singleton local samples have no estimable within-sample variance under the bootstrap; they are not evidence of zero biological uncertainty.

## Remaining limits

This implements borrowing uncertainty for a review sensitivity, not a fully validated latent biomass model or total-error interval. Count sampling, coefficient uncertainty, seasonal abundance effects, persistent site biases and calibration-distribution uncertainty are not fully propagated. The earlier date-pairing sensitivities remain supporting diagnostics; no narrow matching window is imposed and no primary result is selected by direction or significance.

## Reproduce and verify

Run `python3 code/tools/run_pisco_biomass_review.py`. It builds the selected density table, biomass uncertainty outputs, independent R checks and report. Upstream source audits must already exist. Independent R checks verify3,208outer calibration predictions,428held-out intervals,6biomass models and18draw regressions. All487protected manuscript/result hashes remain unchanged. The analysis does not send messages, push changes or modify source data.
'''
(O/'README.md').write_text(text)
page='<h1>Density and biomass borrowing uncertainty</h1><p>Review pipeline:63density comparison-years retained;44biomass years supported conditionally. Paper unchanged.</p><p><a href="README.md">Implemented rules, validation and limitations</a></p><img src="biomass_uncertainty_components.png" style="max-width:100%"><p>A:conditional temporal95%intervals. B:central95%local-size/borrowing draw ranges. B is not a confidence interval; neither panel represents total uncertainty.</p>'
for name,df in [('Coverage and biomass results',summary),('Held-out interval validation',validation),('Visit support and borrowing',v)]:page+='<h2>'+name+'</h2><div style="overflow:auto">'+df.to_html(index=False,float_format=lambda x:f'{x:.5g}')+'</div>'
(O/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Biomass borrowing uncertainty</title><style>body{font:16px/1.5 system-ui;margin:30px}td,th{padding:6px;border-bottom:1px solid #ddd;text-align:left}table{font-size:13px;border-collapse:collapse}</style>'+page)
inputs=[Path('outputs/pisco_timing_donors_20261102/size_sample_inventory.csv'),Path('outputs/pisco_all_months_20261102/visit_values.csv'),Path('outputs/pisco_all_months_20261102/annual_values.csv'),Path('/Users/adrianstier/sbc-donham-kelp-mpa-data/data/ALL_sizefreq_2024.csv')];(O/'input_hashes.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},indent=2))
print('Review report and figure written')
