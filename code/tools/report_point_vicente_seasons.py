from pathlib import Path
import pandas as pd,numpy as np,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path('outputs/point_vicente_seasons_20261101');a=pd.read_csv(O/'annual_input.csv');local=pd.read_csv(O/'local_annual.csv');s=pd.read_csv(O/'local_slopes.csv');p=pd.read_csv(O/'pooled.csv');dates=pd.read_csv(O/'survey_dates_depths.csv');grid=pd.read_csv(O/'coverage_grid.csv');pct=lambda x:100*np.expm1(x)
scs=['Previous_review','AprOct_INNER_MID_full11_anchored','AprOct_INNER_MID_gap90_full11'];colors=['#777777','#287E8D','#C56632'];fig,axes=plt.subplots(1,2,figsize=(12,4.8))
for sc,col in zip(scs,colors):
 q=a[(a.scenario==sc)&(a.MPA=='Point Vicente SMCA')&(a.Source=='PISCO')&(a.Taxa=='M. pyrifera')&(a.Resp=='Bio')&(a.t>0)].sort_values('year');axes[0].plot(q.year,np.exp(q.lnDiff),'-o',color=col,ms=4,lw=1.6)
axes[0].set_yscale('log');axes[0].axhline(1,color='#aaaaaa',lw=.8,ls=':');axes[0].set_xlabel('Survey year');axes[0].set_ylabel('Protected/reference kelp biomass ratio')
choose=[('Previous_review','Previous comparison'),('AprOct_INNER_MID_old10_anchored','Same 10 years'),('AprOct_INNER_MID_full11_anchored','All 11 years'),('AprOct_INNER_MID_late9_anchored','Later 9 years'),('AprOct_INNER_MID_gap90_full11','90-day limit')]
for j,(sc,label) in enumerate(choose):
 r=s[(s.scenario==sc)&(s.Taxa=='M. pyrifera')&(s.Resp=='Bio')].iloc[0];col=colors[0] if j==0 else colors[2] if j==4 else colors[1];axes[1].plot(pct([r.low,r.high]),[j,j],color=col,lw=1.8);axes[1].plot(pct(r.estimate),j,'o',color=col)
axes[1].set_yticks(range(5),[x[1] for x in choose]);axes[1].invert_yaxis();axes[1].axvline(0,color='#aaaaaa',lw=.8,ls=':');axes[1].set_xlabel('Annual change in protected/reference ratio (%)')
for ax,letter in zip(axes,'ab'):ax.text(.02,.98,letter,transform=ax.transAxes,va='top',fontweight='bold',fontsize=13);ax.spines[['top','right']].set_visible(False);ax.margins(y=.15)
fig.tight_layout();fig.savefig(O/'shallow_kelp_history.png',dpi=185);fig.savefig(O/'shallow_kelp_history.pdf');plt.close(fig)
# Dates for the four-station April-Oct shallow history: physical dates, not outcome-weighted means.
refs=['RESORT_POINT','ROCKY_POINT_N','ROCKY_POINT_S'];d=dates[(dates.site.isin(['POINT_VICENTE_W']+refs))&dates.zone.isin(['INNER','MID'])&dates.month.between(4,10)].copy();d['month_fraction']=d.month+(d.day-1)/31
fig,ax=plt.subplots(figsize=(8,4.5))
for sites,col,marker in [(refs,colors[0],'o'),(['ROCKY_POINT_N'],colors[2],'s'),(['POINT_VICENTE_W'],colors[1],'D')]:
 q=d[d.site.isin(sites)].drop_duplicates(['site','year','month','day']);ax.scatter(q.year,q.month_fraction,color=col,marker=marker,s=28,alpha=.8)
ax.set_yticks(range(4,12),['Apr','May','Jun','Jul','Aug','Sep','Oct','Nov']);ax.set_xlabel('Survey year');ax.set_ylabel('Survey date');ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.15);fig.tight_layout();fig.savefig(O/'survey_timing.png',dpi=185);fig.savefig(O/'survey_timing.pdf');plt.close(fig)
rows=['| Point Vicente kelp comparison | Years | Annual ratio change (95% interval) |','|---|---|---:|']
for sc,label in choose:
 r=s[(s.scenario==sc)&(s.Taxa=='M. pyrifera')&(s.Resp=='Bio')].iloc[0];rows.append(f'| {label} | {r.years.replace(";", ", ")} | {pct(r.estimate):+.1f}% ({pct(r.low):+.1f}, {pct(r.high):+.1f}) |')
biotable=['| Response at Point Vicente | Previous review | April–October INNER+MID, 2013–2023 |','|---|---:|---:|']
for tx,resp,label in [('M. pyrifera','Bio','Kelp biomass'),('M. franciscanus','Den','Red-urchin density'),('S. purpuratus','Den','Purple-urchin density'),('P. interruptus','Den','Lobster density')]:
 vals=[]
 for sc in ['Previous_review','AprOct_INNER_MID_full11_anchored']:
  r=s[(s.scenario==sc)&(s.Taxa==tx)&(s.Resp==resp)].iloc[0];vals.append(f'{pct(r.estimate):+.1f}% ({pct(r.low):+.1f}, {pct(r.high):+.1f})')
 biotable.append('| '+label+' | '+' | '.join(vals)+' |')
ptable=['| Pooled response, all 169 histories | Previous review | Recovered shallow history |','|---|---:|---:|']
for tx,resp,label in [('M. pyrifera','Bio','Kelp biomass'),('M. franciscanus','Den','Red-urchin density'),('S. purpuratus','Den','Purple-urchin density'),('P. interruptus','Den','Lobster density')]:
 vals=[]
 for sc in ['Previous_review','AprOct_INNER_MID_full11_anchored']:
  r=p[(p.scenario==sc)&(p.metric=='Post_slope')&(p.Taxa==tx)&(p.Resp==resp)].iloc[0];vals.append(f'{pct(r.estimate):+.2f}% ({pct(r.low):+.2f}, {pct(r.high):+.2f}); q={r.q:.3f}')
 ptable.append('| '+label+' | '+' | '.join(vals)+' |')
text='''# Point Vicente: recovered shallow surveys and seasonal comparability

November 1, 2026. Review only. The manuscript, supplement, canonical outputs and earlier audits remain unchanged.

## Answer

The April surveys recover a continuous shallow-water history for 2013–2023. Habitat coverage can be made consistent, but survey season cannot be fully matched throughout that history. The recovered record is useful sensitivity evidence; it is not yet a seasonally equivalent primary comparison or a causal estimate of protection.

Use POINT_VICENTE_W as the protected station and RESORT_POINT, ROCKY_POINT_N and ROCKY_POINT_S as references, with INNER+MID at every station in every year. Apply April–October to all years and both roles. Recorded depths overlap, approximately 3–12 m. An additional INNER+MID+OUTER sensitivity retains the same stations. The former DEEP-only intersection is not used as a recovery estimate.

## What the source records recover

- In 2013, INNER/MID/OUTER at Point Vicente West were surveyed in April; the June observations were DEEP. May–October had kept only the latter.
- In 2014, Point Vicente West was surveyed in April. That actual survey year was absent from the former paired history. April–October adds it without inventing observations.
- LONG_POINT_E cannot join a fixed INNER+MID history across all eleven years: its 2013 surveys were in March and its 2014 surveys covered only DEEP/OUTER in December. Even admitting all months does not supply the missing 2014 shallow zones.
- Monthly kelp mass-per-frond coefficients were applied to the recorded survey month. That improves biomass conversion; it does not remove seasonal changes in kelp abundance.
- The unknown 2023 UNDERWATER_ARCH frond measurement was not imputed. That station lacks continuous survey coverage independently and does not enter these fixed comparisons.

## Seasonal rules tested consistently

May–October cannot provide a continuous shallow protected history from 2013. April–October provides one protected station and three references. April–November, March–November and all-month rules each permit one protected station and six references. Coverage decisions use survey dates, zones and measurement completeness for all four focal responses, not their abundance or slopes.

The principal season contrasts hold the April–October stations and zones fixed. On those same stations, expanding to November or all months barely changes the kelp slope (about +7.5% rather than +7.6%/year). The late 2015–2023 observations are exactly identical under April–October and May–October. Thus the important changes are recovery of the early shallow surveys and which stations/depths are represented, rather than a new result produced by a cutoff difference during the later years.

These are exploratory seasonal rules motivated by documented filtering problems, not preregistered or established biological thresholds. Applying a new seasonal rule in a primary cross-MPA synthesis would require evaluating the same rule at other affected comparisons; this task tests Point Vicente, without changing their inputs.

## Remaining timing mismatch

For the fixed three-reference April–October comparison, the closest reference dates are 63 days after the protected survey in 2013 and 84 days after it in 2014. Other retained references are as much as 107 and 189 days apart, respectively. Identical month windows therefore do not mean simultaneous surveys.

A separate date audit required every within-zone protected/reference date combination to remain within a common maximum gap, for every year. No fixed full-history reference passed 30 or 60 days. Under a 90-day rule:

- April–October INNER+MID retains ROCKY_POINT_N alone; its largest annual gaps range from 7 to 84 days.
- April–November INNER+MID retains HAWTHORNE_REEF and ROCKY_POINT_N.
- April–November INNER+MID+OUTER retains HAWTHORNE_REEF alone.

Every feasible 90-day case was fit, not just a favorable result. These wider tolerances are diagnostics, not evidence that kelp is seasonally interchangeable over three months. Spring sampling occurs only in the early protected years, so time since protection and seasonal timing cannot be cleanly separated by adding a month term to this short record. The shallower protected surveys were not repeated in summer in those early years; June DEEP surveys cannot supply that calibration.

## Point Vicente kelp results

All estimates are annual percentage changes in the protected/reference ratio, with individual calendar-HAC3 intervals using residual degrees of freedom. They are neither absolute population growth nor causal protection effects.

'''+ '\n'.join(rows)+'''

“Same 10 years” omits 2014 but uses the identical shallow station/zone footprint as “All 11 years.” “Later 9 years” retains that footprint but uses 2015–2023. The original comparison changed its contributing habitat and stations; its +27.4% estimate is not an appropriate direct estimate of shallow recovery. The previously rejected deep-only −27.6% estimate is also not a shallow-water result.

The 90-day April–October pairing estimates about +9.1%/year, with an interval including zero. April–November alternatives estimate about +9.6% (INNER+MID) and +12.9% (including OUTER). Their reference membership differs, so those differences cannot be attributed to season alone. One interval excludes zero; that is not a reason to choose that rule. The early data influence the positive slope, and their seasonal mismatch remains unresolved.

## Other responses on exactly the same shallow footprint

'''+ '\n'.join(biotable)+'''

All four histories use the same sites, zones, years and seasonal rule. The corresponding new urchin/lobster size-based biomass histories were not reconstructed: shared size-pool linkage and April size coverage would require their own validation. The existing biomass histories remain unchanged and must not be described as sharing the newly reconstructed density footprint.

## Pooled consequences

Only these four Point Vicente response histories change. Full-history variants add observed 2014; same-ten-year variants do not. All other MPA histories remain identical. The roster stays at 169 post-protection histories; all 55 before–after effects and variances are unchanged. Recovered post-protection surveys do not create a baseline.

'''+ '\n'.join(ptable)+'''

Pooled conclusions change little. The late Point Dume influence identified previously remains; this local Point Vicente audit does not resolve it. BH q-values reflect the nine-response family and can change for an unchanged response when another test changes. The earlier 140-history coverage proposal is also rerun separately in pooled.csv and remains unadopted.

## Recommendation

Retain the recovered shallow record and its complete sampling metadata in the review pipeline. Treat it as a clearly labeled seasonally imperfect sensitivity, alongside the same-footprint later-year comparison. Do not return to DEEP-only support or cherry-pick April in one year. Before adopting a primary seasonal rule, decide how much temporal separation is biologically defensible, evaluate that rule consistently across the affected sites, and reconcile corresponding biomass support. No new rule is adopted here.

## Figures

shallow_kelp_history.png/pdf: (a) annual protected/reference kelp ratios, gray previous review, teal April–October INNER+MID with three references, orange the 90-day April–October pairing with ROCKY_POINT_N; (b) annual slope estimates with 95% individual calendar-HAC3 intervals. The same-ten-year, all-eleven-year and later-nine-year fits hold stations and zones fixed. Curves connect observed years without imputing gaps. Ratios retain original additive corrections.

survey_timing.png/pdf: survey dates for the fixed shallow station set. Teal diamonds: protected Point Vicente West; gray circles: all three references; orange squares: ROCKY_POINT_N, also included among the gray reference points. Depth zones are INNER+MID; multiple dates at a station/year remain visible. The graph illustrates temporal mismatches rather than claiming a seasonal correction.

## Reproduce and evidence

Run audit_point_vicente_seasons.py, check_point_vicente_date_pairs.py, fit_point_vicente_date_rules.py, Rscript refit_point_vicente_seasons.R, verify_point_vicente_seasons.py and report_point_vicente_seasons.py from code/tools. The raw source is read only. coverage_grid.csv and date_rule_grid.csv record all tested coverage rules; anchored_footprints.csv and date_matched_footprints.csv record exact membership. reconstruction.csv verifies the former four annual responses; source_hashes.json, verification.json and the manifest preserve the audit trail. The annual_input and local_annual files are model-response tables; obsolete normalized role columns are deliberately omitted. Actual role values are in anchored_role_values.csv and date_matched_role_values.csv. Numerical checks do not certify ecological or causal equivalence.
'''
(O/'README.md').write_text(text)
print('\n'.join(rows));print('\n'.join(ptable))
