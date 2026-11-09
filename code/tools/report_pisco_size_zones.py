from pathlib import Path
import pandas as pd,numpy as np,json,html
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path('outputs/pisco_size_zone_resolution_20261102');e=pd.read_csv(O/'effects.csv');d=pd.read_csv(O/'source_zone_details.csv');l=pd.read_csv(O/'coverage_changes.csv')
assert d.all_size_dates_have_raw_same_zone_count.all();assert all(set(map(int,str(v).split(';'))).isdisjoint(range(5,11)) for v in d.raw_zone_months)
sids=list(e.series_id.drop_duplicates());summary=[]
for sid in sids:
 q=e[e.series_id.eq(sid)].set_index('scenario');a=q.loc['current_all_years'];b=q.loc['same_season__current_matched_years'];c=q.loc['same_season__revised_matched_years'];ll=l[(l.series_id==sid)&(l.scenario=='same_season')].iloc[0]
 summary.append(dict(series_id=sid,current_all_years=a.estimate,current_same_years=b.estimate,revised_same_years=c.estimate,coverage_change=b.estimate-a.estimate,size_change=c.estimate-b.estimate,original_years=ll.original_years,retained_years=ll.retained_years,lost_years=ll.lost_years,conditional_low=c.low,conditional_high=c.high,source_resolution='All flagged zone surveys present in raw counts; inconsistent seasonal filtering'))
s=pd.DataFrame(summary);s.to_csv(O/'summary.csv',index=False)
fig,axs=plt.subplots(1,3,figsize=(12,4.5),sharex=True,sharey=True)
for ax,sc,letter in zip(axs,['current_all_years','same_season__current_matched_years','same_season__revised_matched_years'],'ABC'):
 z=e[e.scenario==sc].set_index('series_id').loc[sids];y=np.arange(len(z));ax.hlines(y,z.low,z.high,color='#177c80',lw=2);ax.scatter(z.estimate,y,color='#222',s=23);ax.axvline(0,color='.6',lw=.8);ax.set_yticks(y,[a.split('|')[0].replace(' SMCA','')+' · '+{'M. franciscanus':'red urchin','S. purpuratus':'purple urchin','P. interruptus':'lobster'}[a.split('|')[1]] for a in sids]);ax.set_xlabel('Change in log biomass ratio / year');ax.text(.03,.98,letter,transform=ax.transAxes,va='top');ax.spines[['top','right']].set_visible(False)
axs[0].set_ylim(len(sids)-.5,-.8);fig.tight_layout();fig.savefig(O/'seasonal_size_comparison.png',dpi=180);fig.savefig(O/'seasonal_size_comparison.pdf');plt.close(fig)
body='''# PISCO size-zone flags resolved to seasonal filtering

All 23 flagged species/year/role groups arise because the counts were filtered to May–October but their size pools were not. All 43 underlying station/year/species/zone groups have raw count surveys in the same zone on the size-measurement dates. Those dates are in March, April, November or December. No flag requires a zone relabeling or a missing count record to explain it. This supersedes the earlier unresolved explanation of the 23 flags in A61; it does not establish that every count record in the delivery is complete.

For example, Abalone Cove's 2014 MID-zone sizes and counts were both measured on December 22; only the count records were removed. Point Vicente West's 2013 shallow-zone sizes and counts were both measured on April 24, leaving only June DEEP-zone counts after filtering. Three Palms East's 2021 shallow-zone surveys were on November 2, twelve days after the retained October DEEP survey.

## Same-observation checks

I tested five changes separately: removing only the flagged zone measurements; requiring same-site sizes; requiring same-site/zone sizes; applying the existing May–October rule to the size pools; and requiring site, zone and season together. These are diagnostic restrictions, not adopted exclusions. The count data, offsets, researcher exclusions and models remain fixed. Each revised fit has a current-input fit on exactly the same retained years. Annual mean masses are reweighted from source frequencies and used to scale current biomass, so these remain conditional sensitivities rather than a complete production rebuild.

Applying May–October to the size pools and requiring site/zone/season gives identical recovery estimates for these six histories. All six retain matching-season support in most years, but three species-specific paired biomass years would be lost: Abalone Cove lobster2014, Abalone Cove purpleurchin2018 and Point Vicente redurchin2013. A whole paired biomass year is withheld if a positive role lacks eligible sizes; density remains supported. Zero counts remain zero.

The largest changes in these comparisons come from losing years, rather than from revising average animal size within the same years. Abalone Cove lobster changes from0.232 to0.190 simply by comparing the ten supported years, then to0.184 when sizes are restricted. Its temporal interval already includes zero after the year is removed, before changing sizes. Point Vicente redurchin changes from0.016 to−0.061 when2013 is removed, then to−0.063 with the size restriction; intervals include zero throughout. On matched years, the largest absolute slope adjustment from seasonal size matching is0.0091 log-ratio units/year. The other four histories show small changes in slope, with the same signs.

## Recommendation

Use a consistent survey-selection rule for counts and their size measurements. These findings do not independently justify May–October as the right biological season: some excluded visits fall only days outside its boundary. Before freezing the primary analysis, choose the count-survey timing policy and apply it to the corresponding size support. Do not exclude these sites because of the earlier zone flags, or borrow out-of-season sizes silently to preserve a biomass year. Retain valid density when biomass lacks support. The existing user hold on selecting a new seasonal cutoff remains in force.

The figure separates the current full history (A), current inputs on the years retained by the seasonal size check (B), and season-matched sizes on those same years (C). Lines are conditional temporal-model95% intervals, not total biomass uncertainty. These are post-protection recovery slopes using calendar years, not estimates of before–after protection effects. No pooled or downstream pathway model has been changed.

## Verification

All23original flags and43underlying station/species/zone groups were traced to raw files. Every flagged size date has a matching raw swath zone/date; all those count dates are outside May–October. Independent R checks reproduce66slopes and temporal standard errors. All487protected paper/result files are unchanged. Source hashes, pool counts, exact survey dates, coverage changes and matched-year fits are saved here. Code: `code/tools/audit_pisco_size_zones.py`, `code/tools/verify_pisco_size_zones.R`, `code/tools/report_pisco_size_zones.py`.
'''
(O/'README.md').write_text(body)
page='<h1>PISCO size-zone check</h1><p>All 23 flags trace to inconsistent seasonal filtering: matching count surveys are present in the raw files.</p><p>Paper unchanged. No new seasonal rule adopted.</p><p><a href="README.md">Findings, methods and recommendation</a></p><img style="max-width:100%" src="seasonal_size_comparison.png"><p>A: current full history. B: current values on supported years. C: season-matched sizes on those same years. Lines show conditional temporal uncertainty, not total biomass uncertainty.</p>'
for name,df in [('Summary',s),('Source records',d),('Coverage',l),('All fits',e)]:page+='<h2>'+name+'</h2><div style="overflow:auto">'+df.to_html(index=False,float_format=lambda x:f'{x:.5g}')+'</div>'
(O/'index.html').write_text('<!doctype html><meta charset="utf-8"><title>PISCO size-zone resolution</title><style>body{font:16px/1.5 system-ui;margin:30px}td,th{padding:6px;border-bottom:1px solid #ddd;text-align:left}table{font-size:13px;border-collapse:collapse}</style>'+page)
print(s.to_string(index=False))
