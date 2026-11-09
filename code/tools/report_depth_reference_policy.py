from pathlib import Path
import pandas as pd,numpy as np,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path('outputs/depth_reference_policy_20261101');m=pd.read_csv(O/'mpa_decisions.csv').fillna({'lost_years':'','added_years':'','years':''});p=pd.read_csv(O/'pair_decisions.csv');s=pd.read_csv(O/'selected_pair_zone_years.csv');s=s[s.scenario=='Measured_shallow_d90_n5'].copy();s['difference']=s.depth_center_p-s.depth_center_r
names=['Abalone Cove SMCA','Gull Island SMR','Point Dume SMCA','Point Vicente SMCA','South Point SMR'];short=lambda x:x.replace(' SMCA','').replace(' SMR','')
rr=[]
for (mpa,ps,rs,z),q in s.groupby(['MPA','site_p','site_r','zone']):
 rr.append(dict(MPA=mpa,protected=ps,reference=rs,zone=z,n=len(q),median_difference=q.difference.median(),depth_difference_slope_m_per_year=np.polyfit(q.year-q.year.min(),q.difference,1)[0]))
dr=pd.DataFrame(rr);dr.to_csv(O/'selected_depth_drift.csv',index=False)
# Distinct pair summaries, not independent replicates or uncertainty intervals.
fig,axs=plt.subplots(1,2,figsize=(11,4.8),sharey=True)
for ax,z,letter in zip(axs,['INNER','MID'],'ab'):
 for i,name in enumerate(names):
  q=dr[(dr.MPA==name)&(dr.zone==z)];ax.scatter(q.median_difference,np.full(len(q),i),color='#287E8D',s=24,alpha=.65);ax.scatter(q.median_difference.median(),i,marker='D',color='#222222',s=40)
 ax.axvline(0,color='#aaaaaa',lw=1,ls=':');ax.set_xlabel('Protected minus reference depth (m)');ax.set_yticks(range(5),[short(x) for x in names]);ax.text(.02,1.02,f'{letter}   {z.title()}',transform=ax.transAxes);ax.spines[['top','right']].set_visible(False)
axs[0].invert_yaxis();fig.tight_layout();fig.savefig(O/'depth_differences.png',dpi=190);fig.savefig(O/'depth_differences.pdf');plt.close(fig)
rows=['| MPA | Previous paired years | Proposed paired years | Removed from this recovery set | Recovered |','|---|---:|---:|---|---|']
for name in names:
 r=m[(m.MPA==name)&(m.scenario=='Measured_shallow_d90_n5')].iloc[0];rows.append(f'| {short(name)} | {len(r.original_years.split(";"))} | {r.n_years} | {r.lost_years.replace(";", ", ") or "None"} | {r.added_years.replace(";", ", ") or "None"} |')
refs=['| MPA | Protected stations retained | Reference stations retained | Qualifying pairs (not independent samples) |','|---|---|---|---:|']
for name in names:
 r=m[(m.MPA==name)&(m.scenario=='Measured_shallow_d90_n5')].iloc[0];refs.append(f'| {short(name)} | {r.protected.replace(";", ", ")} | {r.reference.replace(";", ", ")} | {r.n_pairs} |')
sens=['| MPA | 30-day limit | 60-day limit | Proposed 90-day limit | Depth labels with unknown depths allowed, 90 days |','|---|---:|---:|---:|---:|']
for name in names:
 r=m[m.MPA==name].set_index('scenario');vals=[r.loc[k,'n_years'] for k in ['Measured_shallow_d30_n5','Measured_shallow_d60_n5','Measured_shallow_d90_n5','Zone_only_d90_n5']];sens.append('| '+short(name)+' | '+' | '.join(map(str,vals))+' |')
# A duration-only whole-review screen is illustrative, not a complete new primary cohort.
f=pd.read_csv('outputs/full_followup_recovery_20261101/effects.csv');f=f[(f.scenario=='Full')&(f.metric=='Post_slope')].copy();f['passes_duration_only']=(f.na>=5)&(f.last_after-f.first_after>=5);f[['series_id','MPA','Source','Taxa','Resp','na','first_after','last_after','passes_duration_only']].to_csv(O/'whole_review_duration_screen.csv',index=False)
text='''# Recommended primary recovery selection rule

November 1, 2026. This is an outcome-blind coverage recommendation, not an adopted analysis. No new biological effects were fitted, and the manuscript, supplement and canonical results remain unchanged. The scoped numerical audit covers the five PISCO assignments with documented seasonal-zone problems; the same policy would need coverage verification elsewhere before claiming a fully revised primary cohort.

## Recommendation

Keep the inherited reference groups and documented researcher exclusions. Retain all eligible protected/reference combinations within each group; do not choose a new single reference because it provides the longest series or a preferred trend.

For primary post-protection recovery, I recommend:

1. Use actual dated surveys in all months, paired within the same calendar year and provider survey campaign. Within each depth stratum, choose the closest available dates, at most 90 days apart. Both strata must come from the same campaign; where multiple campaigns occur in a calendar year, select the tightest joint date match. This avoids mixing January observations from an earlier campaign with fall observations from the next one. Calendar time remains the clock.
2. Use the consistently defined shallow INNER and MID survey strata, with at least two usable transects per retained visit. Require measured positive depths within the shallow target and overlapping observed depth support for each retained station pair/stratum. Never substitute DEEP for missing shallow surveys. Retain missing-depth records for a separately labeled sensitivity rather than inventing depths.
3. Require at least five paired post-protection years spanning at least five elapsed years for each contributing pair. Keep shorter histories as supporting evidence. This is a proposed recovery-support minimum, not Emily's baseline criterion. Do not impose an early-start requirement; state the actual observed interval and never extrapolate to protection year.
4. Fit at the station/stratum level and standardize over shared observed depth support, accounting for date/season and changes in contributing stations. Do not average the selected pairs as independent data: they reuse stations and surveys. Keep reference dependence at MPA/group level. Constant habitat differences and changing sampling depth are different problems.
5. Retain before–after change as a separate question, using the already selected baseline policy and fixed baseline support. A recovered post-only history does not become a before–after comparison.

The 90-day cap is a declared quarter-year operational tolerance, not a proven safe interval for kelp seasonality. Use 60 days as the stricter sensitivity, with 30 days shown too. The shallow ceiling 16.25 m is the midpoint between the published 12.5 m middle and 20 m deep targets; it is our target definition, not a provider requirement. Using 15 m or 20 m instead gives the same retained station pairs and MPA-year sets in this audit. Freeze these rules before the coordinated outcome refit and do not alter them in response to effect direction or significance.

## What the coauthor advice supports

Scott Hamilton's explanation distinguishes habitat/geography from data quality: some otherwise good reference records are inappropriate for a local paired comparison. Emily described preserving the same stations across before/after periods while treating after-only coverage separately. Jenn described expert and community-composition evidence for the Naples reference choice, but stressed that no close match exists. Those messages support respecting the established map and documenting its limitations; they do not validate our new 90-day or depth rules, nor every inherited reference individually.

Dan Pondella also identified sedimentation, turbidity and runoff as alternative influences at Abalone Cove, Point Vicente and Malibu. Matching dates and depth cannot remove those environmental differences. No new site exclusion is inferred from that comment. Detailed dated message evidence stays in the private manuscript correspondence record 2026-11-01_depth_reference/EVIDENCE.md. No emails were sent.

The [monitoring methods](https://opc.dataone.org/metacat/d1/mn/v2/object/doi:10.25494%2FP6%2FMLPA_kelpforest_methods.6) describe non-fixed benthic transects, approximate 5/12.5/20 m targets and generally June–October surveys. Actual sampling departs from that general schedule. A numeric transect identifier does not establish a fixed physical transect across years. Fish sampling differs and is not certified by this swath audit.

## Retained years

The comparison is against the latest audited annual histories, not a newly updated manuscript. Counts are unique calendar years represented by at least one qualifying pair, not years common to every station and not independent replicates. Station-level modeling is required before these can replace annual averages.

'''+ '\n'.join(rows)+'''

Exact years and identities are in mpa_decisions.csv and selected_pair_zone_years.csv. “Removed” means outside this proposed measured-shallow recovery analysis, not deleted from the repository or declared bad biological data. Other responses/programs at these MPAs remain separate.

## Reference choices

'''+ '\n'.join(refs)+'''

At Abalone Cove, 3 Palms East and Whites Point each have four qualifying years; Point Fermin has two. They fail the proposed recovery-support rule rather than an outcome test. Bunker Point is about 3.6 km from Abalone Cove.

At Gull, both protected Gull stations and both Forney references remain. These references are approximately 13.7–14.8 km away: they are inherited best-available comparisons, not demonstrably close habitat matches. Yellowbanks remains excluded. Gull 2004 has no measured depths; 2015 lacks qualifying same-campaign shallow support. The proposed set recovers 2020 and 2022.

Point Dume retains LECHUZA/NICHOLAS_CANYON_W and the actual 2013, 2019–2023 surveys. Leo Carrillo has only one qualifying year. The mapping lists two coordinates for LECHUZA, so its distance is left unresolved rather than selected by row order. Same-year pairing does not cure the six-year gap, winter-to-fall sampling shift or depth change. Point Dume is a candidate for an adjusted recovery fit, not a validated general recovery estimate simply because six years exist.

Point Vicente retains its existing two protected stations and nine mapped references where they meet the rule. No new outside reference is introduced; years differ among station pairs. UNDERWATER_ARCH's unknown 2023 kelp measurement remains excluded, not recovered. At South Point, the two South Point stations pair with Cluster Point North/South and Johnsons Lee South West. Chickasaw/Trancion and Johnsons Lee North combinations have only 2–3 paired years. Johnsons Lee South East has longer coverage but lacks overlapping paired-record stratum-center depth ranges; it remains a depth-limited sensitivity rather than being called ecologically invalid.

## What the depths tell us

The earlier single-pair differences were partly consequences of narrowing the reference group. Across the proposed eligible pairs, Gull's typical protected-minus-reference depth differences are about +0.8 m in INNER and +1.1 m in MID. South Point's typical INNER difference is about −1.1 m, rather than treating the earlier −2.3 m single-pair value as representative of every available reference.

Depth is not completely stable through time. In the six-year Point Dume record, the protected-minus-reference MID difference changes by roughly 0.29 m/year in a descriptive linear summary. Some Point Vicente pair/stratum differences also drift. These are measurements of sampled habitat, not biological effects. A station intercept alone cannot remove a changing depth difference; the proposed outcome model needs actual depth and seasonal terms, with estimability, influence and extrapolation checks before adoption.

One nominal MID transect at POINT_FERMIN in May 2016 records 27.4 m. It lies outside the shallow target; the original value is preserved and not called a typo. Its partner transect and other years do not justify silently assigning it a nominal depth. The coordinate conflict and missing depths are also retained as explicit metadata questions.

Hard trimming every transect to the exact empirical protected/reference depth-range intersection was tested and rejected as the primary rule. With two transects, losing one just outside a sampled endpoint can remove an entire zone-year. That threshold reflects sample extrema, not an ecological discontinuity. The proposed model instead uses recorded depths and restricts standardized predictions to overlapping support. Common_depth remains a conservative coverage sensitivity; this is not evidence that its excluded observations were wrong.

## Sensitivity of coverage

'''+ '\n'.join(sens)+'''

A zero means no pair passes the five-year/five-span criterion, not no surveys exist. The three-year alternatives and all candidate reasons are retained in pair_decisions.csv. The retained pair and year sets are identical for the 15, 16.25 and 20 m shallow ceilings under the 90-day/five-year rule; that choice is not driving the current keep/drop list.

## What is and is not complete

This establishes a concrete selection proposal and its consequences before refitting outcomes. It does not reuse the previous single-pair pooled estimates as results for this new multi-reference design. All previous biological results remain frozen.

Across the broader 169-history recovery review, a five-year/five-span duration screen alone flags 21 existing histories across Farnsworth, Matlahuayl, Point Dume SMCA/SMR and Swamis. The recovered four Point Dume SMCA swath responses could regain eligibility; that is not a final 152-history primary roster. Other programs, fish and size-based biomass have not received this new depth/date screen, and their sampled populations and size linkage need verification. Before–after eligibility remains a separate decision.

The next implementation is one coordinated station/stratum recovery refit with constant prediction support, shared-reference dependence, actual calendar dates and seasonal/depth diagnostics. Where seasonal or depth changes cannot be separated from recovery without extrapolation, retain that history as supporting evidence. Especially at Point Dume, do not equate a coverage pass with identification of a protection effect. Present the outcome of that refit for review before updating the paper.

## Reproduce and figure

Run code/tools/audit_depth_reference_policy.py, verify_depth_reference_policy.py and report_depth_reference_policy.py. PLAN.md records the candidate rule and coverage-stage refinement. source_hashes.json pins the source inputs; verification.json reports independent raw-depth, date, membership and protection checks. Every eligible pair and nonqualifying alternative is retained.

depth_differences.png/pdf: each teal point is one eligible pair's median protected-minus-reference depth difference; black diamonds are medians across those pair summaries. These are descriptive distributions, not confidence intervals or independent replicates. INNER and MID refer to original survey strata. Positive values mean protected surveys were deeper. No outcome trend is plotted.
'''
(O/'README.md').write_text(text)
print('\n'.join(rows));print('\n'.join(sens))
