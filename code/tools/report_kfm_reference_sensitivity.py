from pathlib import Path
import pandas as pd,numpy as np,json,html
R=Path(__file__).resolve().parents[2];B=R/'outputs/integrated_review_20261103';O=R/'outputs/kfm_reference_sensitivity_20261103'
names=['main','without_webster','without_miracle','without_wyckoff','matched_main_without_wyckoff','without_flagged','matched_main_without_flagged']
def folder(n):return B if n=='main' else B/'sensitivities'/('kfm_reference_'+n)
coverage=[];pooled=[];local=[];annual=[]
base=pd.read_csv(B/'annual_values.csv');basekeys=set(zip(base.loc[base.response_supported,'series_id'],base.loc[base.response_supported,'year']))
for n in names:
 p=folder(n);a=pd.read_csv(p/'annual_values.csv');m=pd.read_csv(p/'model_eligibility.csv');keys=set(zip(a.loc[a.response_supported,'series_id'],a.loc[a.response_supported,'year']))
 coverage.append(dict(scenario=n,supported_years=len(keys),lost_supported_years=len(basekeys-keys),new_supported_years=len(keys-basekeys),fitted_histories=int(m.eligible_model.ne('insufficient_post_years').sum()),before_after=int(m.eligible_model.eq('before_after_and_post_trend').sum()),post_only=int(m.eligible_model.eq('post_trend_only').sum())))
 pp=pd.read_csv(p/'pooled_results.csv');pp['scenario']=n;pooled.append(pp)
 j=pd.read_csv(p/'joint_model_results.csv');j=j[j.analysis.eq('estimated_AR1') & j.Source.eq('KFM') & j.MPA.isin(['Harris Point SMR','Santa Barbara Island SMR'])];j['scenario']=n;local.append(j)
 a=a[a.Source.eq('KFM') & a.MPA.isin(['Harris Point SMR','Santa Barbara Island SMR'])];a['scenario']=n;annual.append(a)
coverage=pd.DataFrame(coverage);coverage.to_csv(O/'coverage_changes.csv',index=False)
pooled=pd.concat(pooled,ignore_index=True);pooled.to_csv(O/'pooled_results.csv',index=False)
local=pd.concat(local,ignore_index=True);local.to_csv(O/'local_results.csv',index=False)
pd.concat(annual,ignore_index=True).to_csv(O/'affected_annual_values.csv',index=False)
keys=['Source','Taxa','Resp','model','term'];decomp=[]
for variant in ['without_wyckoff','without_flagged']:
 frames={n:pooled[pooled.scenario.eq(n)].set_index(keys) for n in ['main','matched_main_'+variant,variant]}
 for k in set.intersection(*[set(v.index) for v in frames.values()]):
  b,matched,new=[frames[n].loc[k] for n in ['main','matched_main_'+variant,variant]]
  decomp.append(dict(zip(keys,k),scenario=variant,full_main=b.estimate,matched_main=matched.estimate,omission=new.estimate,coverage_change=matched.estimate-b.estimate,reference_change=new.estimate-matched.estimate,total_change=new.estimate-b.estimate))
pd.DataFrame(decomp).to_csv(O/'pooled_decomposition.csv',index=False)
# Individual coefficients elsewhere may change through reestimated shared covariance.
comparison=local[local.Taxa.eq('S. pulcher')&local.Resp.eq('Den')&local.model.eq('post_trend')&local.term.eq('recovery_slope')].copy()
for col in ['estimate','low','high']:comparison[col+'_pct']=100*np.expm1(comparison[col])
comparison.to_csv(O/'sheephead_trends.csv',index=False)
# Exact reference names and yearly role values for the email/figure.
sites=pd.read_csv(R.parent/'sbc-donham-kelp-mpa-data/config/kfm_physical_sites.csv').drop_duplicates('site').set_index('site').source_name.to_dict()
st=pd.read_csv(B/'station_annual_values.csv');st=st[st.series_id.str.contains('Harris Point|Santa Barbara Island')&st.series_id.str.contains('S. pulcher')&st.series_id.str.endswith('|KFM|Den')];st['site_name']=st.site.map(sites);st.to_csv(O/'sheephead_station_values.csv',index=False)
print(coverage.to_string(index=False));print(comparison[['scenario','MPA','n_years','estimate_pct','low_pct','high_pct']].to_string(index=False))
# Narrative below is the reviewed Nov3 interpretation; fail rather than print stale values.
expected_local={('main','Harris Point SMR'):-20.8,('matched_main_without_wyckoff','Harris Point SMR'):-16.5,('without_wyckoff','Harris Point SMR'):-13.8,('without_miracle','Harris Point SMR'):-20.4,('main','Santa Barbara Island SMR'):3.4,('without_webster','Santa Barbara Island SMR'):-0.4,('without_flagged','Santa Barbara Island SMR'):3.2}
for (scenario,mpa),expected in expected_local.items():
 value=comparison.loc[comparison.scenario.eq(scenario)&comparison.MPA.eq(mpa),'estimate_pct'].iloc[0]
 assert round(value,1)==expected, 'Update narrative for changed local result'
assert coverage.set_index('scenario').loc['without_flagged','lost_supported_years']==217
for scenario,expected in [('main',229.3),('without_flagged',424.2)]:
 value=pooled.loc[pooled.scenario.eq(scenario)&pooled.Source.eq('All programs')&pooled.Taxa.eq('P. interruptus')&pooled.Resp.eq('Den')&pooled.model.eq('flat_before_after')&pooled.term.eq('average_observed_post_change'),'estimate'].iloc[0]
 assert round(100*np.expm1(value),1)==expected, 'Update narrative for changed before-after result'
labels={'main':'Current review','without_webster':'Omit Webster’s Arch','without_miracle':'Omit Miracle Mile','without_wyckoff':'Omit Wyckoff Ledge','without_flagged':'Omit all flagged references','matched_main_without_wyckoff':'Current references, Wyckoff-test years','matched_main_without_flagged':'Current references, combined-test cohort'}
def fmt(v):return f'{v:+.1f}' if abs(v)<10000 else f'{v:+.2g}'
def pctcell(row):return f"{fmt(100*np.expm1(row.estimate))}% ({fmt(100*np.expm1(row.low))}, {fmt(100*np.expm1(row.high))})"
short=coverage[coverage.scenario.isin(['main','without_webster','without_miracle','without_wyckoff','without_flagged'])].copy();short.scenario=short.scenario.map(labels)
rows=[]
for (tx,resp),g in pooled[pooled.Source.eq('All programs')&pooled.model.eq('post_trend')&pooled.term.eq('recovery_slope')].groupby(['Taxa','Resp']):
 rows.append({'Variable':tx+' '+resp,**{labels[n]:pctcell(g[g.scenario.eq(n)].iloc[0]) for n in ['main','matched_main_without_flagged','without_flagged']}})
pooled_table=pd.DataFrame(rows);pooled_table.to_csv(O/'pooled_trends_readable.csv',index=False)
local_table=comparison[['scenario','MPA','n_years','estimate_pct','low_pct','high_pct']].copy();local_table.scenario=local_table.scenario.map(labels)
text='''# KFM reference-site sensitivity

November 3, 2026. These tests answer Scott Gabara’s concerns about the current integrated review. They do not change the manuscript, primary roster or researcher exclusions.

## Findings

- **Harris Point:** Hare Rock is the protected station. Wyckoff Ledge is the only current reference for six responses; sheephead also uses Miracle Mile. Omitting Miracle Mile changes the sheephead post-protection ratio trend from −20.8% to −20.4% per year. Omitting Wyckoff leaves no reference for the six other responses. Sheephead retains 15 post-protection years, with no before baseline. On those same 15 years, the current-reference trend is −16.5% (95% interval −25.6 to −6.4), compared with −13.8% (−24.2 to −2.0) using Miracle Mile alone. The larger full-record change from −20.8% partly reflects lost years.
- **Santa Barbara Island:** Webster’s Arch affects sheephead density only; its benthic exclusion already applies. Removing it retains Arch Point and Cat Canyon as references and changes the local sheephead trend from +3.4% (−3.8 to +11.2) to −0.4% (−7.9 to +7.6) per year. Protected stations are SE Sea Lion Rookery, Graveyard Canyon and Southeast Reef, as available each year. No years are lost.
- **Rodes Reef:** no current measurement uses; omitting it adds no change. Gull KFM remains excluded under the current researcher-scope rules; Gull PISCO is separate and unchanged.
- **Combined omission:** removing Wyckoff, Miracle Mile and Webster’s Arch removes all seven Harris Point histories, loses 217 supported comparison-years, and reduces fitted histories from 162 to 155. Before–after eligible histories decline from 62 to 57. These are response histories and uses, not independent sites or unique calendar years.
- **Pooled post-protection trends:** kelp changes from +16.9% to +16.5% per year; purple-urchin density from −9.7% to −10.3%; lobster density from +7.1% to +10.4%; sheephead density from −3.1% to −1.6%. Predator trend intervals still include increases and decreases. Red-urchin density changes only modestly (−9.6% to −10.2%), although its interval endpoint crosses zero; do not turn that boundary change into a different biological conclusion.
- **Before–after lobster:** the flat-baseline average observed post contrast changes materially, from +229% (−40 to +1695) to +424% (+77 to +1451). The matched-cohort calculation attributes this to removing Harris Point, not to changing the remaining reference averages. Allowing the pre-protection trend to continue gives +158% (−97 to +19050) after omission. The evidence remains sensitive to the counterfactual model.
- **Before–after kelp:** the flat-baseline contrast changes from +280% (+56 to +827) to +292% (+63 to +843). The pre-trend alternative remains uncertain: about +10% before omission and −9% afterward, with very wide intervals in both cases.

## What the tests isolate

Four omission cases reuse the integrated estimator and joint uncertainty model. Dates, retained within-station observations, zone support, frozen ratio offsets, nonlinear size draws and within-station aggregation stay fixed. Equal weights among remaining stations are renormalized. A year without both roles cannot produce a paired estimate. No alternative station is introduced and no outcome determines the omission.

For Wyckoff and the combined case, a matched-main refit uses the original reference definition on exactly the alternative’s supported series/year keys. Full main → matched main measures the consequence of changed coverage. Matched main → omission measures the reference-definition change on that coverage, including reestimated shared process covariance. These differences add on the log-estimate scale; they are not separately identified causal mechanisms. Pooled weights and nuisance covariance can change after histories are removed.

Annual values and measurement draws for unaffected histories are unchanged. Their fitted coefficients can still move through reestimated covariance within a program/species/response group. This explains why the SBI estimate after combined omission (+3.2%, interval −19.1 to +31.8) differs from omitting Webster alone. It is not a second change to SBI’s observations. In the combined case, the KFM sheephead post-trend AR estimate reaches its upper bound of 0.95; fixed-AR alternatives also span zero at SBI. Treat this local uncertainty as dependence-sensitive.

Do not interpret the local trend as density change at the protected station alone: it is change relative to the reference. These are not habitat-validation tests. Miracle Mile is also on the southern side that Scott questioned, so its survival as a numerical comparison does not make it a preferred ecological replacement.

## Recommendation and remaining question

Keep these as explicit sensitivities while asking Scott whether he recommends retaining or omitting the Hare Rock versus Wyckoff/Miracle Mile comparison, and whether the Webster’s Arch habitat concern also warrants excluding its fish records. Do not choose by which result better supports a cascade. Rodes and the existing benthic exclusions need no further numerical repair. No new reference omission has been adopted here.

## Figure captions

**Local sheephead trends.** Points and lines show estimates and conditional 95% intervals for annual change in the protected/reference sheephead density ratio, using the unchanged zero-handling offsets. Panel A separates full-record, common-year and reference changes at Harris Point. Panel B shows the Webster-only and combined omissions at SBI. Current means the integrated review, not the unchanged manuscript. The intervals are for individual estimates, not tests of the difference between scenarios. No before–after result is plotted.

**Pooled trends.** Gray circles: current full cohort. Ochre squares: current references after removing the histories/years absent from the combined test. Green triangles: combined reference omission. Intervals are conditional 95% intervals. Overlapping ochre/green estimates show where the change is explained by loss of Harris Point rather than changes to remaining references. Positive values mean an increasing protected/reference ratio. No before–after result is plotted.

## Verification and reproducibility

All four omissions and two matched-main joint refits completed. Checks confirm exact retained measurement records, unchanged unaffected annual values and draws, fixed offsets, matched support keys, absence of omitted references, original ecological exclusion reasons and unchanged primary numerical files. All coefficient covariance matrices are symmetric and nonnegative within numerical tolerance; none of the 60 model groups per scenario reaches an upper variance bound. Depending on the scenario, 11–13 groups reach the AR search boundary and 36–38 have a zero variance component. These remain conditional fits with the same diagnostic qualifications as the current review; they are not estimates of total uncertainty.

The initial numerical run was followed by metadata-only ledger cleanup: existing exclusion reasons were preserved and dropped-row station metadata labeled as original candidate coverage. Neither changes model values or draws. The rerun script includes that cleanup from the start.

Run `python3 code/tools/run_kfm_reference_sensitivity.py`, then `check_kfm_reference_sensitivity.py`, `report_kfm_reference_sensitivity.py` and `plot_kfm_reference_sensitivity.py` from the analysis root. Individual source/draw files stay in the existing private sensitivity folders. `verification.json`, `checks.csv`, `diagnostics.csv` and the input records accompany this report. Exact scenarios and matching definitions are in `design.json`.
'''
(O/'README.md').write_text(text)
body='<h1>KFM reference-site tests</h1><p>November 3, 2026 · Sensitivities only; primary analysis and paper unchanged.</p><p><a href="README.md">Findings, methods and limitations</a> · <a href="SCOTT_DRAFT.md">Draft note to Scott (not sent)</a> · <a href="../dataset_inventory_20261101/index.html">Current decisions</a></p>'
body+='<h2>Coverage</h2>'+short.to_html(index=False)+'<h2>Sheephead comparisons</h2><img src="sheephead_reference_trends.png" alt="Harris Point and Santa Barbara Island sheephead trends">'
body+='<p>Conditional estimates and 95% intervals. Current means the integrated review. Harris uses 19 or 15 post-protection years as labeled; SBI uses 20. The common-year comparisons separate coverage loss from reference changes. These are relative density trends, not before–after protection effects.</p>'
body+='<h2>Pooled post-protection trends</h2><img src="pooled_reference_trends.png" alt="Pooled trends under current, matched-cohort and omitted-reference cases"><p>Gray circles: current full cohort. Ochre squares: current references on the combined test’s cohort. Green triangles: combined omission. Conditional 95% intervals; no before–after effects plotted.</p>'+pooled_table.to_html(index=False)
body+='<h2>A material before–after change</h2><p>The flat-baseline lobster estimate rises from +229% (95% interval −40 to +1695) to +424% (+77 to +1451), entirely through removing Harris Point. The pre-trend sensitivity remains very uncertain: +158% (−97 to +19050). Do not treat this as robust causal evidence or choose exclusions from the result.</p>'
body+='<h2>Detailed tables</h2><ul>'+''.join(f'<li><a href="{f}.csv">{html.escape(f.replace("_"," "))}</a></li>' for f in ['coverage_changes','sheephead_trends','sheephead_station_values','local_results','pooled_results','pooled_decomposition','diagnostics','checks'])+'</ul><p>AR and zero-variance boundary flags remain. Read the methods and limitations before interpreting the intervals.</p>'
(O/'index.html').write_text('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>KFM reference tests</title><style>body{font:16px/1.5 system-ui;max-width:1200px;margin:32px auto;padding:20px;color:#243832}img{width:100%}table{border-collapse:collapse;display:block;overflow:auto;font-size:13px}th,td{padding:9px;border-bottom:1px solid #ccd8d0;text-align:left}a{color:#245d4c}</style></head><body>'+body+'</body></html>')
