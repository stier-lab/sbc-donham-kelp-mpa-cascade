from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path('outputs/pisco_all_months_20261102');e=pd.read_csv(O/'effects.csv');a=pd.read_csv(O/'annual_values.csv');cov=pd.read_csv(O/'coverage.csv');info=json.loads((O/'verification.json').read_text());su=pd.read_csv(O/'size_matching.csv');sv=pd.read_csv(O/'station_years.csv')
sc='allmonths_sharedzones_visit_sizes';summary=[]
for sid in a.series_id.drop_duplicates():
 q=e[(e.series_id==sid)&(e.scenario==sc)].set_index('comparison');c=cov[(cov.series_id==sid)&(cov.scenario==sc)].iloc[0];old=q.loc['current_all_years'];new=q.loc['revised_matched_years'] if 'revised_matched_years' in q.index else None;matched=q.loc['current_matched_years'] if 'current_matched_years' in q.index else None
 summary.append(dict(series_id=sid,current_slope=old.estimate,current_matched_slope=matched.estimate if matched is not None else np.nan,revised_slope=new.estimate if new is not None else np.nan,revised_low=new.low if new is not None else np.nan,revised_high=new.high if new is not None else np.nan,original_post_years=c.original_years,supported_post_years=c.retained_years,lost_years=c.lost_years,status='Review estimate; conditional temporal uncertainty' if new is not None else 'Insufficient supported years for a slope'))
s=pd.DataFrame(summary);s.to_csv(O/'summary.csv',index=False)
density=s[s.series_id.str.endswith('|Den')];sids=list(density.series_id);fig,axs=plt.subplots(1,3,figsize=(12,4.5),sharex=True,sharey=True)
for ax,stage,letter in zip(axs,['original','allmonths_oldpool','allmonths_sharedzones_visit_sizes'],'ABC'):
 rows=[]
 for sid in sids:
  z=e[(e.series_id==sid)&(e.scenario==('allmonths_oldpool' if stage=='original' else stage))&(e.comparison==('current_all_years' if stage=='original' else 'revised_matched_years'))].iloc[0];rows.append(z)
 z=pd.DataFrame(rows);y=np.arange(len(z));ax.hlines(y,z.low,z.high,color='#177c80',lw=2);ax.scatter(z.estimate,y,c='#222',s=24);ax.axvline(0,c='.6',lw=.8);ax.set_yticks(y,[sid.split('|')[0].replace(' SMCA','')+' · '+{'M. franciscanus':'red urchin','S. purpuratus':'purple urchin','P. interruptus':'lobster'}[sid.split('|')[1]] for sid in sids]);ax.set_xlabel('Change in log density ratio / year');ax.text(.03,.98,letter,transform=ax.transAxes,va='top');ax.spines[['top','right']].set_visible(False)
axs[0].set_ylim(5.5,-.8);fig.tight_layout();fig.savefig(O/'density_all_months.png',dpi=180);fig.savefig(O/'density_all_months.pdf');plt.close(fig)
body='''# All-month PISCO survey rerun

Adrian approved retaining valid visits from every calendar month, matching sizes to counts where possible, averaging repeat visits consistently and checking timing differences. This review implements that rule for the six animal biomass comparisons implicated in the23zone flags (Abalone Cove and Point Vicente: redurchins, purpleurchins and lobsters), alongside their density responses. It keeps the existing retained sites, station membership and observed year windows. This is not a full-paper or all-program rerun. Researcher exclusions and manuscript outputs remain unchanged.

## Rules actually used

- No calendar-month exclusion. Existing program/method filters remain. The survey roster includes true zeros, reconstructed from visits on which other taxa were observed. Missing measurements are not zeros.
- Average transects within a dated zone visit, visits within zone/year, zones within station/year, and stations within protected/reference role. Thus repeat visits do not increase a station's weight. Keep method and calendar date in the survey identity.
- INNER and MID are present at every retained station-year in these six histories. The shared-zone version holds that set constant across protected/reference stations and years. Source zone categories do not establish identical measured depths; no claim of exact habitat equivalence is made.
- For biomass, search the original size file for same-site/zone/species/date donors, preserving frequency multiplicity and the existing25mm urchin threshold. If none exist, use the nearest dated sample within the same site/zone/year, record the gap and compare against exact-date-only support. Ties use all tied donors. No other site, zone or year is silently substituted. This matching is an explicit review assumption, not proof that a long date gap is biologically unimportant.
- A positive count visit without eligible local sizes leaves station biomass unknown. Withhold that paired biomass year rather than silently omit a reference station or average around its missing visit. Density remains available.

## What changes

All63existing post-protection comparison-years remain available for density. The locally matched biomass version supports25of63comparison-years: Abalone Cove red3/11,purple2/11,lobster9/11; Point Vicente red0/10,purple4/10,lobster7/10. Two histories therefore have insufficient supported years for a slope. This is a diagnostic of stricter local size support, not an instruction to remove the original density histories or a claim that pooled-size biomass is inherently invalid.

Among801positive species/zone/visit uses in the shared zones,624have same-date sizes,24use another visit in that site/zone/year (up to294days away), and153lack eligible local sizes. These are survey uses across comparisons/species, not801independent surveys. The exact-date-only scenario shows the influence and coverage cost of using other visits. Dates and donor sample sizes are listed individually.

Searching the raw size file was essential: it recovered measurements for four Point Vicente2013species/zone visits that were absent from the old seasonal pool (lobster and purpleurchin in INNER/MID). An audit based only on previously retained donors would falsely report those measurements as missing. Remaining gaps are kept explicit; small animals below the retained size threshold do not supply eligible adult biomass measurements.

The density rerun shows declining protected/reference urchin ratios and an increasing lobster ratio at Abalone Cove. Point Vicente shows a declining purpleurchin ratio, with redurchin and lobster trend intervals spanning zero. These are recovery trends in protected/reference ratios, not causal before–after effects of protection. The tables report point estimates and conditional temporal intervals; total sampling/allometric uncertainty is not included.

## Changes separated

Five versions are saved: May–October with equal repeat-visit weighting and existing size pools; all months with those same pools; all months and constant shared zones with those pools; shared zones with local same/nearest-visit sizes; and shared zones with exact-date sizes only. Old-pool versions are sensitivity comparisons, not certified local size support. Restored positive counts without an old donor pool remain unknown in those versions. Every scenario includes the original fit on its exact supported years, distinguishing a changed estimate from a changed time window. Deterministic mean donor mass replaces cached Monte Carlo point estimation in these diagnostics, so the first reconstruction step also includes that numerical difference.

Figure: A=current density inputs; B=all-month density with equal repeat-visit averaging; C=all-month density in constant shared INNER/MID zones. All panels use identical years for each history. Lines are95%conditional temporal-model intervals.

## Timing still matters

Including November or April surveys avoids an arbitrary calendar cutoff but does not guarantee contemporaneous protected/reference sampling. Across all retained years, Abalone protected stations were surveyed roughly41–50days earlier on average than references, with some annual gaps near127days; Point Vicente mean differences were about3–7days, with annual gaps up to76days. The role timing trends also differ. These descriptive mean day-of-year summaries can hide within-year spread and assume these calendar-year survey schedules; they are not a seasonal correction. Inspect the dated visits before a final paired-survey timing model. Do not label this issue resolved merely because every month is now eligible.

## Verification and holds

Independent R checks reproduce172slopes/temporal standard errors and2625station-scenario-year aggregations from dated visit means. The comparison tables preserve all researcher exclusions and all487protected manuscript/result hashes. No new source file, manuscript, primary pooled model or pathway model was changed. The six KFM alias histories are unrelated and remain provisional.

Reproduce with `python3 code/tools/refit_pisco_all_months.py`, `Rscript code/tools/verify_pisco_all_months.R`, then `python3 code/tools/report_pisco_all_months.py`. Raw record identities remain in their private source files. These results support dropping the hard month cutoff for this review, while keeping local-size limitations and protected/reference timing differences explicit before primary integration.
'''
(O/'README.md').write_text(body)
page='<h1>All-month PISCO rerun</h1><p>Six affected comparisons. All existing density years retained; local biomass size support varies. Paper unchanged.</p><p><a href="README.md">Rules, findings and remaining limits</a></p><img src="density_all_months.png" style="max-width:100%"><p>A:current density; B:all months and equal repeat-visit weights; C:all months and shared INNER/MID zones. Lines: conditional temporal95%intervals.</p>'
for name,df in [('Results and supported years',s),('Size match records',su),('Timing differences',pd.read_csv(O/'timing_summary.csv')),('All scenarios',e)]:page+='<h2>'+name+'</h2><div style="overflow:auto">'+df.to_html(index=False,float_format=lambda x:f'{x:.5g}')+'</div>'
(O/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>All-month PISCO rerun</title><style>body{font:16px/1.5 system-ui;margin:30px}td,th{padding:6px;border-bottom:1px solid #ddd;text-align:left}table{font-size:13px;border-collapse:collapse}</style>'+page)
inputs=[D for D in [Path('/Users/adrianstier/sbc-donham-kelp-mpa-data/data/PISCO/MLPA_kelpforest_swath_2024.csv'),Path('/Users/adrianstier/sbc-donham-kelp-mpa-data/data/ALL_sizefreq_2024.csv'),Path('outputs/biomass_support_20261102/private/pool_members.csv'),Path('outputs/all_program_station_refits_20261102/verified_station_values.csv'),Path('outputs/biomass_resolution_20261102/review_annual_input.csv')]]
(O/'input_hashes.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},indent=2));print(s.to_string(index=False))
