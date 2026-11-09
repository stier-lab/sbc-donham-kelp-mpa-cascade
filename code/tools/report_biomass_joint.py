from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,html
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path('outputs/biomass_joint_checks_20261102')
e=pd.read_csv(O/'conditional_joint_results.csv');m=pd.read_csv(O/'count_size_matching.csv');h=pd.read_csv(O/'count_size_history_summary.csv');info=json.loads((O/'verification.json').read_text())
# Companion artwork contains identifiers and axes only; explanation is in the report.
kelp=e[e.series_id.str.contains('M. pyrifera',regex=False)];sids=kelp.series_id.unique()
fig,ax=plt.subplots(1,3,figsize=(11,4),sharey=True,sharex=True)
for aa,scheme in zip(ax,['count_only','size_only','joint']):
 g=kelp[kelp.scheme==scheme].set_index('series_id').loc[sids];y=np.arange(len(g));aa.hlines(y,g.measurement_draw_low,g.measurement_draw_high,color='#177c80',lw=2);aa.scatter(g.estimate,y,c='#222222',s=24,zorder=3);aa.axvline(0,color='.7',lw=.8);aa.set_yticks(y,[s.split('|')[0].replace(' SMR 2003','').replace(' SMR','') for s in sids]);aa.set_xlabel('Change in log biomass ratio / year');aa.text(.03,.98,{'count_only':'A','size_only':'B','joint':'C'}[scheme],transform=aa.transAxes,va='top');aa.spines[['top','right']].set_visible(False);aa.grid(axis='x',alpha=.15)
ax[0].set_ylim(len(sids)-.6,-.75);fig.tight_layout();fig.savefig(O/'kelp_sampling.png',dpi=180);fig.savefig(O/'kelp_sampling.pdf');plt.close(fig)
# Every review result is tabulated, including held histories and existing direct measurement checks.
direct=pd.read_csv('outputs/biomass_support_20261102/direct_support_summary.csv');direct['current_scope']='Existing direct measurement screen; new joint audit covers KFM kelp only. Direct fish and LTER kelp require their own nested sampling model.';direct.to_csv(O/'direct_measurement_scope.csv',index=False)
m[m.fraction_matching_recorded_zone<m.fraction_zone_evaluable-1e-12].to_csv(O/'size_zone_mismatches.csv',index=False)
text='''# Biomass matching and joint sampling checks

Update: [the subsequent source check](../pisco_size_zone_resolution_20261102/README.md) explains all23PISCO zone flags. Raw counts exist on the size dates; seasonal filtering removed those counts but retained the sizes. The original flags below describe filtered coverage, not missing raw surveys.

This is a review-only audit. The manuscript, primary estimates, researcher exclusions and annual coverage are unchanged.

## What is now checked

- **Count reconstruction:** all 4,236 station-year density checks reproduce the retained data. Count zeros remain in the survey roster. Reused stations receive the same resampling weights across comparisons.
- **Count–size matching:** 1,143 donor-bearing protected/reference year groups in 37 pooled-animal histories are listed individually. Six histories fail the source-identity guard. Matching labels elsewhere does not certify ecological representativeness.
- **KFM kelp:** 65,215 delivered rows contain 18,413 exact duplicates. Removing them leaves 46,802 unique site/date/plant-marker records, with no conflicting records under those keys. This validates the existing exact-row cleanup, not deduplication by equal stipe counts. Monthly conversion factors stay paired with each plant.
- **Joint sampling:** 1,000 conditional draws for 34 histories (29 animal and five KFM kelp). Six alias-affected histories are withheld; two Point Dume histories cannot support a post-protection slope with fewer than three observed years. These holds affect this diagnostic, not density retention or the frozen paper.
- **Regression guard:** known excluded source IDs and the two misleading combined urchin labels cannot enter this new review resampling. The check is source-table specific: original kelp plants at Wyckoff Ledge are not automatically rejected because an unrelated combined urchin label is ambiguous. Production rebuilding remains on hold; this is not yet a production-pipeline patch.

## What the matching tells us

PISCO has 380 of 529 donor-bearing role-years with at least some borrowed or unresolved-site sizes. In 23 role-years, evaluable same-site size samples include a zone absent from the count roster. Same-site nearest survey-date gaps reach 290 days. These are flags for assessment, not automatically excluded observations. A nearest survey date does not show that every sampled count animal was measured.

LTER has 18 of 128 donor-bearing role-years with some borrowed sizes. Evaluable same-site donors match a count-survey date, but that does not validate the borrowed portion. Only 18 PISCO role-years have any evaluable same-site numeric depth pairs; KFM/LTER count deliveries lack corresponding depth fields. Numeric depth agreement is reported only where both sources provide depths; zone names across programs are not treated as equivalent. PISCO's current 25 mm size rule is checked separately from unfiltered KFM/LTER sizes; changing only measured sizes would not standardize the counted populations.

## What the sampling checks show

The five KFM kelp histories retain their current slope directions in the conditional combined draw ranges. That is evidence about the available count and plant samples, not a new significance test. Count resampling can matter substantially: Matlahuayl lobster and red-urchin slopes have broad combined draw ranges that cross zero. All comparisons use their existing observed post-protection years and fixed offsets; no baseline or site is selected by its result.

Figure panels A, B and C show count sampling, plant sampling and their combination, respectively. Dots are current estimates; segments are the central 95% of conditional measurement draws. They are **not confidence intervals for protection effects**. Existing temporal-model intervals are reported separately in the result table.

## Assumptions and remaining limits

Counts are resampled at the retained spatial sampling units: KFM plot IDs and LTER transects have shared weights through time; PISCO transects are resampled within site/year/zone, preserving equal zone weighting. This is a diagnostic assumption about spatial sampling, not proof that every plot label follows an unchanged physical plot. Each shared animal-size source row receives the same weight wherever it is used, with recorded animal multiplicity preserved. Kelp plants are resampled within site/year, with their survey-month conversion intact. Count and size weights are independent conditional on the retained samples; unsupplied covariance between abundance and sampled sizes is not estimated.

The analysis does not resolve wrong-site donors, unmeasured sizes, allometric coefficient error, survey detection, season/depth bias, or dependence among plants measured during the same visit. Singleton samples cannot reveal size variance. Direct LTER lobster, fish and LTER/PISCO kelp already have separate measurement screens; their full joint uncertainty is not completed by this pooled-size/KFM-kelp analysis. Do not add the measurement variance to the existing temporal variance mechanically: that can count part of the same variability twice. A coordinated observation/model uncertainty analysis remains necessary before calling an interval total uncertainty.

The exact expected KFM kelp conversion (count density times mean measured plant mass) is also checked against the cached simulated point estimate. Differences reach about 3% in individual positive station-years. The expectation avoids Monte Carlo noise; adoption requires the coordinated rebuild, not an isolated manuscript edit.

## What can be resolved now versus what needs evidence

The donor alias guard, same-observation count/size audit, KFM plant identity check and conditional joint draws are implemented. Keep supported density. Resolve the six alias histories from Scott's original site-specific sizes before final biomass estimates. Investigate the listed PISCO zone mismatches and borrowed LTER pools before adopting restrictions. No new sample-size, season or depth cutoff has been chosen.

## Verification and reproduction

Run `Rscript code/tools/export_biomass_joint_sources.R`, `python3 code/tools/audit_biomass_matching.py`, `python3 code/tools/audit_biomass_joint.py`, `Rscript code/tools/verify_biomass_joint.R`, and `python3 code/tools/report_biomass_joint.py`. Three guard tests pass. Independent R checks reproduce 34 baseline slopes and temporal standard errors, 306 sampled-draw slopes, and 296 KFM station-year expected masses. All 487 protected paper/result hashes remain unchanged. Record-level extracts remain private.
'''
(O/'README.md').write_text(text)
tables={'History decisions':pd.read_csv(O/'resampling_admission.csv'),'Count–size matching':h,'Conditional sampling results':e,'PISCO zone flags':pd.read_csv(O/'size_zone_mismatches.csv'),'KFM plant support':pd.read_csv(O/'kelp_plant_support.csv')}
body='<h1>Biomass matching and sampling checks</h1><p>Review only. Paper unchanged. Draw ranges are not total confidence intervals.</p><p><a href="README.md">Full findings, assumptions and limits</a></p><img src="kelp_sampling.png" style="max-width:100%"><p>Dots: current recovery slopes. Lines: central 95% of conditional measurement draws. Panels A, B and C show counts, plant sizes and combined sampling, respectively, using identical observed years.</p>'
for name,df in tables.items():body+='<h2>'+html.escape(name)+'</h2><div style="overflow:auto">'+df.to_html(index=False,float_format=lambda x:f'{x:.5g}')+'</div>'
(O/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>Biomass sampling audit</title><style>body{font:16px system-ui;margin:30px;line-height:1.5}td,th{padding:6px;text-align:left;border-bottom:1px solid #ddd}table{border-collapse:collapse;font-size:13px}</style>'+body)
inputs=[Path('outputs/biomass_resolution_20261102/review_annual_input.csv'),Path('outputs/biomass_support_20261102/private/pool_members.csv'),Path('outputs/all_program_station_refits_20261102/stage/reconstruction.RData'),Path('/Users/adrianstier/sbc-donham-kelp-mpa-data/data/MBON/KFM_Macrocystis_RawData_1984-2023.csv'),Path('data/external/sbc1275/Macrocystis_Biomass_Relationships_20261016.csv')]
(O/'input_hashes.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},indent=2))
print('Wrote report and figure')
