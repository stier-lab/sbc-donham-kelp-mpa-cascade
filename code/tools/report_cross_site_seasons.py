from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path('outputs/cross_site_seasons_20261101');c=pd.read_csv(O/'selected_coverage.csv').fillna({'years':'','lost_original':'','added_years':''});p=pd.read_csv(O/'pooled.csv');f=pd.read_csv(O/'local_effects.csv');e=pd.read_csv(O/'effects.csv');dates=pd.read_csv(O/'survey_dates_depths.csv')
mpas=['Abalone Cove SMCA','Gull Island SMR','Point Dume SMCA','Point Vicente SMCA','South Point SMR'];short={m:m.replace(' SMCA','').replace(' SMR','') for m in mpas};pct=lambda x:100*np.expm1(x)
scs=['AprOct_IM_d0','AprOct_IM_d60','AprOct_IM_d90','AprNov_IM_d0','AprNov_IM_d60','AprNov_IM_d90','AllMonths_IM_d60'];labels=['Apr–Oct\nNo date limit','Apr–Oct\n≤60 days','Apr–Oct\n≤90 days','Apr–Nov\nNo date limit','Apr–Nov\n≤60 days','Apr–Nov\n≤90 days','All months\n≤60 days']
mat=np.array([[c[(c.MPA==m)&(c.scenario==s)].n.iloc[0] for s in scs] for m in mpas]);fig,ax=plt.subplots(figsize=(12,4.6));ax.imshow(mat,cmap='Blues',vmin=0,vmax=20,aspect='auto')
for i in range(5):
 for j in range(7):ax.text(j,i,str(mat[i,j]),ha='center',va='center',fontsize=14,color='white' if mat[i,j]>=12 else '#222222')
ax.set_yticks(range(5),[short[m] for m in mpas]);ax.set_xticks(range(7),labels);ax.tick_params(length=0);ax.spines[:].set_visible(False);fig.tight_layout();fig.savefig(O/'retained_years.png',dpi=190);fig.savefig(O/'retained_years.pdf');plt.close(fig)
# Fixed physical station pairs across every displayed timing-rule contrast.
show=['Previous_review','Anchor_AprOct_IM_d60','Anchor_AprOct_IM_d90','Anchor_AprNov_IM_d60','Anchor_AprNov_IM_d90'];labs=['Previous review','Apr–Oct, ≤60 days','Apr–Oct, ≤90 days','Apr–Nov, ≤60 days','Apr–Nov, ≤90 days'];responses=[('M. pyrifera','Bio','Kelp biomass'),('M. franciscanus','Den','Red-urchin density'),('S. purpuratus','Den','Purple-urchin density'),('P. interruptus','Den','Lobster density')]
fig,axs=plt.subplots(2,2,figsize=(12,7.5))
for ax,(tx,resp,label),letter in zip(axs.ravel(),responses,'abcd'):
 for j,sc in enumerate(show):
  q=p[(p.rule==sc)&(p.variant==('Previous_review' if sc=='Previous_review' else 'Recovered'))&(p.Taxa==tx)&(p.Resp==resp)].iloc[0];color='#777777' if j==0 else '#287E8D';ax.plot(pct([q.low,q.high]),[j,j],color=color,lw=1.5);ax.scatter(pct(q.estimate),j,color=color,s=30)
 ax.set_yticks(range(5),labs);ax.invert_yaxis();ax.axvline(0,color='#bbbbbb',ls=':',lw=1);ax.set_xlabel('Annual change in protected/reference ratio (%)');ax.text(.02,1.02,f'{letter}   {label}',transform=ax.transAxes,fontsize=11);ax.spines[['top','right']].set_visible(False);ax.margins(y=.17)
fig.tight_layout();fig.savefig(O/'pooled_changes.png',dpi=190);fig.savefig(O/'pooled_changes.pdf');plt.close(fig)
# Exact selected-pair depth overlap is a diagnostic, not a new outcome-driven gate.
rows=[]
for r in c[c.scenario=='AprNov_IM_d90'].itertuples():
 for y in map(int,r.years.split(';')):
  for z in ['INNER','MID']:
   q=dates[(dates.year==y)&dates.zone.eq(z)&dates.month.between(4,11)];a=q[q.site==r.protected];b=q[q.site==r.reference]
   rows.append(dict(MPA=r.MPA,protected=r.protected,reference=r.reference,year=y,zone=z,protected_min=a.depth_min.min(),protected_max=a.depth_max.max(),reference_min=b.depth_min.min(),reference_max=b.depth_max.max(),depth_known=bool(np.isfinite([a.depth_min.min(),a.depth_max.max(),b.depth_min.min(),b.depth_max.max()]).all()),overlap=(max(a.depth_min.min(),b.depth_min.min())<=min(a.depth_max.max(),b.depth_max.max()) if np.isfinite([a.depth_min.min(),a.depth_max.max(),b.depth_min.min(),b.depth_max.max()]).all() else np.nan)))
pd.DataFrame(rows).to_csv(O/'selected_depth_overlap.csv',index=False)
def cell(q):
 r=q.iloc[0];return f'{pct(r.estimate):+.1f}% ({pct(r.low):+.1f}, {pct(r.high):+.1f}); q={r.q:.3f}'
coverage=['| MPA | Previous years | Apr–Oct, 60 d | Apr–Oct, 90 d | Apr–Nov, 60 d | Apr–Nov, 90 d |','|---|---:|---:|---:|---:|---:|']
for m in mpas:
 old=len(c[c.MPA==m].iloc[0].original_years.split(';'));nn=[str(int(c[(c.MPA==m)&(c.scenario==s)].n.iloc[0])) for s in ['AprOct_IM_d60','AprOct_IM_d90','AprNov_IM_d60','AprNov_IM_d90']];coverage.append('| '+short[m]+f' | {old} | '+' | '.join(nn)+' |')
pt=['| Response | Previous review | Apr–Oct, 60 d | Apr–Nov, 60 d | Apr–Nov, 90 d |','|---|---|---|---|---|']
for tx,resp,label in responses:
 vals=[cell(p[(p.rule==sc)&(p.variant==('Previous_review' if sc=='Previous_review' else 'Recovered'))&(p.Taxa==tx)&(p.Resp==resp)]) for sc in ['Previous_review','Anchor_AprOct_IM_d60','Anchor_AprNov_IM_d60','Anchor_AprNov_IM_d90']];pt.append('| '+label+' | '+' | '.join(vals)+' |')
local=['| MPA | Previous kelp trend | Fixed pair, Apr–Nov ≤90 d | Years lost / added relative to previous |','|---|---|---|---|']
for m in mpas:
 b=e[(e.variant=='Previous_review')&(e.MPA==m)&(e.Source=='PISCO')&(e.Taxa=='M. pyrifera')].iloc[0];r=f[(f.scenario=='AprNov_IM_d90')&(f.variant=='Recovered')&(f.MPA==m)&(f.Taxa=='M. pyrifera')].iloc[0];cr=c[(c.MPA==m)&(c.scenario=='AprNov_IM_d90')].iloc[0];local.append(f'| {short[m]} | {pct(b.estimate):+.1f}% | {pct(r.estimate):+.1f}% ({pct(r.low):+.1f}, {pct(r.high):+.1f}) | Lost {cr.lost_original or "none"}; added {cr.added_years or "none"} |')
at=['| Step, fixed Apr–Nov ≤90 d pair | Kelp | Red-urchin density | Purple-urchin density | Lobster density |','|---|---|---|---|---|']
for v in ['Retention_only','Old_common_years','Corrected_common_years','Recovered']:
 vals=[cell(p[(p.rule=='Anchor_AprNov_IM_d90')&(p.variant==v)&(p.Taxa==tx)&(p.Resp==resp)]) for tx,resp,_ in responses];at.append('| '+v.replace('_',' ')+' | '+' | '.join(vals)+' |')
text='''# Cross-site seasonal and survey-date audit

November 1, 2026. Review only: manuscript, supplement, sources and canonical analyses unchanged.

## Findings

April–October is too restrictive to retain all shallow histories, but widening the season alone does not ensure matched survey dates. April–November restores useful years at Gull Island, South Point and Point Dume. A 60-day rule still costs many years at Abalone Cove and removes Point Vicente's 2014 survey; 90 days retains more data but does not establish seasonal equivalence.

The effects do not uniformly disappear. On the fixed-pair April–November/90-day diagnostic, pooled kelp is +10.5%/year, lobster density +8.7%, purple-urchin density −8.3%, and red-urchin density −1.3%. The corresponding adjusted q-values are .080, .122, .110 and 1.000. This supports sensitivity of the inferred trajectories to sampling support, not a newly established cascade. The controlled steps below show most of these shifts arise from rebuilding observations on fixed shallow support rather than recovering a few extra years.

## Scope and interpretation

The same rules were evaluated at all five PISCO MPA assignments with documented seasonal zone losses: Abalone Cove, Gull Island, Point Dume SMCA, Point Vicente and South Point. Researcher exclusions remain. The rebuild covers kelp biomass and red-urchin, purple-urchin and lobster densities. It does not establish matching seasonal coverage for fish, size-based invertebrate biomass, KFM, LTER or satellite records.

Rules use fixed physical stations and identical INNER+MID zones on both roles; INNER+MID+OUTER is also tested. Windows are May–October, April–October, April–November, March–November and all months. Maximum within-zone separation between every protected/reference survey date is unrestricted, 30, 60 or 90 days. These thresholds are diagnostics, not established biological equivalence limits. All retained visits count; we do not select convenient visits or silently impute missing surveys. This maximum-gap rule is conservative: admitting more months can introduce a distant extra visit and cause a year to fail the timing limit, so an all-month window need not retain more years. Monthly kelp mass conversion does not correct seasonal kelp abundance.

## What is kept or lost

Counts below are observed years from the best-covered fixed station pair for each rule. Pair choice uses number of years, span, earliest year and alphabetical IDs, in that order, before fitted effects. Thus the table shows recoverable coverage, not an adopted new reference assignment. Exact years, physical stations, candidate alternatives, gaps and reasons are in selected_coverage.csv, pair_coverage.csv and pair_year_gaps.csv. Two observations cannot supply a slope; three allow estimation but provide weak recovery evidence.

'''+ '\n'.join(coverage)+'''

The full-history target includes every post-designation year through 2023 with any raw swath observations in both roles, including formerly filtered years. It does not invent annual observations in calendar gaps. Requiring complete target coverage under April–November/90 days leaves only Point Vicente among these five MPAs. Other affected swath histories fail because of missing shallow coverage and/or unmatched dates. This is not a recommendation to discard their other data.

## What admitting winter surveys recovers

All-month sampling with the same 60-day date limit retains 3, 11, 6, 10 and 14 years at Abalone Cove, Gull Island, Point Dume, Point Vicente and South Point, respectively, with coverage-selected fixed pairs. Point Dume recovers 2013, 2019 and 2020, giving six actual years through 2023 rather than a three-year late record. This retains the original LECHUZA/NICHOLAS_CANYON_W pair and INNER+MID zones. The early observations are in winter and later observations in fall, so matching protected and reference dates within a year does not remove changes in sampling season through time. The large gap between 2013 and 2019 remains.

In that Point Dume six-year diagnostic, kelp changes +4.2%/year (95% interval −24.8 to +44.4), red-urchin density −6.3% (−17.3 to +6.3), purple-urchin density −39.9% (−49.0 to −29.2), and lobster density −14.7% (−45.8 to +34.1). Thus the very steep three-year kelp increase is not sustained across the longer observed period. Across the full pooled synthesis, this all-month/60-day alternative gives kelp +7.9%, red urchins −9.7%, purple urchins −17.2%, and lobster +11.4%, with q=.070, .070, .047 and .220. It also changes reference choice at some other MPAs and is not a pure month-window contrast. These values are not a reason to select this rule.

## Controlled pooled comparison

For the table below, each MPA uses its April–November/90-day coverage-selected pair under every displayed rule. This holds reference identity fixed: it is separate from the maximum-coverage table, where a different pair may retain more years. For example, Abalone Cove's anchored pair has only two years within 60 days, so its affected response slopes drop from those anchored pooled cases even though another pair provides four years.

All other histories retain the previous review estimates and variances. Failed affected histories are explicitly omitted from these sensitivities. Original additive offsets and shared-reference clusters remain fixed. Pooled intervals use the existing REML and cluster-robust procedure; q is BH over nine responses. These are annual changes in protected/reference ratios, not absolute population growth or causal protection effects. Other responses can change q through multiplicity without changing their own estimates.

'''+ '\n'.join(pt)+'''

## Changes within MPAs

The following uses the coverage-selected April–November/90-day pair. Point Vicente has one reference here, unlike the earlier three-reference April–October diagnostic. A numerical difference between those reviews is therefore not solely a seasonal-window effect.

'''+ '\n'.join(local)+'''

## Separating retention, years and reconstructed observations

The following controlled sequence uses the same April–November/90-day pairs. Retention-only pools the original effects for retained histories. Old-common-years refits original annual values on the intersection of old/new years. Corrected-common-years uses corrected observations on exactly those years. Recovered uses all eligible observed years. Each step re-estimates variances and pooled weights, so differences are sequential attribution, not an additive decomposition of pure biological components.

'''+ '\n'.join(at)+'''

All 55 qualifying before–after effects and variances remain unchanged: the 20 targeted swath histories have no qualifying before–after fit in the previous review. Recovered post-protection surveys do not create a baseline. Full-history-only scenarios and all tested windows/zones are retained in pooled.csv; do not read them as primary recommendations.

## Remaining limits

Fixed zone names do not certify habitat equivalence. In the April–November/90-day selected pairs, Gull's protected MID zone is typically about 1.8 m deeper than the reference, while South Point's protected INNER zone is about 2.3 m shallower. Some differences are small; Gull 2004 has no depth measurements, which is unknown rather than a failed depth match. Depth ranges alone do not establish a biological exclusion threshold. No depth-based post hoc exclusion was added. selected_depth_overlap.csv reports actual depth ranges for every selected April–November/90-day pair/year/zone. Missing counts/fronds stay unknown; taxa absent from an otherwise observed event are zeros. No April invertebrate mass is invented from an unvalidated size pool. Narrowing to a single fixed pair changes the spatial population represented; the chosen pair is a coverage diagnostic, not a collaborator-approved substitute for the existing reference group. Date thresholds can remove early years or recent endpoints and change the recovery period.

## Recommendation and diagnostic limits

Do not adopt a blanket seasonal exclusion solely because it yields a cleaner-looking sample or a favorable result. Retain the year-round source records and the complete coverage ledger. These tests support evaluating paired survey dates and measured habitat together, with seasonally matched subsets as sensitivity checks. A final rule still needs a biological justification for allowed date separation and an agreed spatial target; maximizing pair coverage cannot supply that justification. Point Dume should remain visibly separate as a sparse, season-changing recovery history even when six years are recoverable.

Numerical reconstruction, coverage selection, interval/BH arithmetic, unchanged histories and source checks are verified. Model diagnostics remain partial: max-Cook influence is recorded, calendar-HAC covariance is retained, and short/exact-fit failures are explicit, but linear recovery is an approximation and these checks do not resolve seasonal confounding, depth equivalence or all size-sampling uncertainty. No new primary policy or paper change is made.

## Figures and reproducibility

retained_years.png/pdf gives numbers of retained years in the best-covered fixed pair, all INNER+MID; darker blue means more years, with exact values printed. No inference is attached to the colors.

pooled_changes.png/pdf gives pooled annual protected/reference ratio trends and 95% intervals, using the same anchored pair within each MPA across timing rules. Gray is the previous review and teal the diagnostic alternatives. The scope is the four reconstructed responses at five MPAs; other components remain unchanged.

Run code/tools/audit_cross_site_seasons.py, refit_cross_site_seasons.R, verify_cross_site_seasons.py, then report_cross_site_seasons.py. PLAN.md records choices made before inspecting refitted outcomes. Raw source paths can be supplied with DONHAM_DATA_REPO_PISCO. Source hashes, exact role values, memberships, model rows, residual influence checks and numerical verification are saved here. The continuing findings register is docs/AUDIT_FINDINGS_REGISTER.md. No new primary rule is adopted.
'''
(O/'README.md').write_text(text)
print('\n'.join(pt));print('\n'.join(at));print('Depth overlaps:',pd.DataFrame(rows).groupby('MPA').overlap.agg(['sum','count']).to_string())
