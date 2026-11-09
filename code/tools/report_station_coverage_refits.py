from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path('outputs/station_coverage_refits_20261101');e=pd.read_csv(O/'effects.csv');z=pd.read_csv(O/'station_zone_summary.csv');i=pd.read_csv(O/'reference_omission.csv')
wide=e.pivot(index=['series_id','MPA','Resp'],columns='scenario',values='estimate').reset_index();wide['station_adjustment']=wide.Log_station_adjusted-wide.Log_station_unadjusted;wide['reference_subset_change']=wide.Persistent_reference-wide.Approved;wide.to_csv(O/'comparison_summary.csv',index=False)
selected=['Anacapa Island SMR 2003','Point Dume SMCA','Scorpion SMR','South Point SMR'];fig,axs=plt.subplots(1,2,figsize=(10,4.5),sharey=True)
for ax,resp in zip(axs,['Bio','Den']):
 for y,mpa in enumerate(selected):
  for scenario,dy,color,face in [('Approved',-.12,'#777777','white'),('Persistent_reference',.12,'#007b79','#007b79')]:
   q=e[(e.MPA==mpa)&(e.Resp==resp)&(e.scenario==scenario)].iloc[0];ax.errorbar(q.estimate,y+dy,xerr=[[q.estimate-q.low],[q.high-q.estimate]],fmt='o',mfc=face,color=color,ms=5,capsize=2)
 ax.axvline(0,color='#bbbbbb',lw=.8);ax.set_yticks(range(4),['Anacapa','Point Dume SMCA','Scorpion','South Point']);ax.text(.98,.97,{'Bio':'Biomass','Den':'Density'}[resp],transform=ax.transAxes,ha='right',va='top');ax.set_xlabel('Annual change in log protected/reference ratio');ax.spines[['top','right']].set_visible(False);ax.set_ylim(-.5,3.5)
fig.tight_layout();fig.savefig(O/'reference_comparison.png',dpi=180);fig.savefig(O/'reference_comparison.pdf');plt.close(fig)
caption='Gray open points use available reference stations; teal filled points use reference stations present throughout the approved history. Protected observations, post-protection years and ratio offsets are identical within each comparison. Bars are 95% calendar-gap HAC3 t intervals. The four assignments with changing reference membership and a persistent subset are shown. Abalone Cove and Santa Barbara Island have no full-history persistent reference subset. Positive slopes indicate a rising protected/reference ratio, not necessarily an increase inside alone.'
(O/'FIGURE_CAPTION.md').write_text(caption+'\n')
text='''# Station coverage and reference influence

November 1, 2026. PISCO sheephead diagnostic only; paper and canonical outputs unchanged. This follows the approved CRANE/repeat/complete-pass rebuild and the user's coverage decisions.

## What was compared

All 34 sheephead response histories were refitted over their actual post-protection years. First, replace the available reference average with stations present in every approved year; hold protected means, offsets and fitted years constant. This is possible for 30 histories; Abalone Cove and Santa Barbara Island have no persistent reference station over the full approved record. The approved record includes before-years where available; the fitted recovery slopes use only post-years. No station was chosen for its result.

Second, omit each reference station in turn and compare both treatments on exactly the same remaining years. Eighty-two such comparisons show influence, not justification for removing a station.

Third, fit weighted log-station observations with year intercepts and a protected-by-calendar-time interaction, with and without station intercepts on identical rows. Each role/year has equal total weight. This checks whether persistent station differences explain an apparent trend when membership changes. It is a separate estimand: an average of log station values is not the log of an arithmetic mean. Differences between these models and the approved ratio model therefore cannot be attributed solely to station adjustment. Only the adjusted-minus-unadjusted contrast within this model isolates station intercept adjustment. Constant station intercepts do not remove habitat-specific trends or seasonal/depth changes.

## Findings

Point Dume SMCA is most sensitive to station-composition adjustment. Using Nicholas Canyon West throughout the same five observed post-protection years changes sheephead biomass slope from −0.0502 to −0.1733 log units/year and density from +0.0250 to −0.0870. All four 95% intervals span zero. On the separate log-station scale, adding station intercepts changes biomass +0.0784→−0.1442 and density +0.1076→−0.0641. These diagnostics support treating its recovery estimate as sensitive to station coverage; they do not establish a protection effect or justify exclusion by direction.

The persistent reference subset changes Anacapa biomass −0.0023→−0.0350 and density −0.0603→−0.0640; Scorpion biomass −0.0065→−0.0037 and density −0.0066→−0.0080; South Point biomass −0.0188→−0.0156 and density −0.0075→+0.0035. These use identical years within each response. Review the intervals in effects.csv, not just signs.

Reference identity also matters even when membership is constant. Omitting Cuyler West at Harris Point changes biomass slope −0.1349→+0.0898 and density −0.0564→+0.0499 on the same 17 post-years. Omitting Bunker Point reduces Abalone's biomass slope +0.3423→+0.2126 and density +0.1744→+0.0878 on the same 11 post-years. These are influence checks; no researcher-approved reference is removed. Reference-omission intervals and exact years are preserved in reference_omission.csv. They are many exploratory checks, not multiplicity-adjusted claims.

## Shared zones: correction to the broader screen

The preceding role-level screen found shared zones in every fish year. Checking every contributing station is stricter: Abalone has no zone shared across all stations in 2010 or 2014, so a full-history constant common zone is impossible under unchanged membership. Those are coverage findings, not bad-data labels. Other assignments have at least one zone shared at every station throughout, but often a narrower set: Anacapa, Point Vicente and South Point reduce to OUTER; Santa Barbara Island to INNER. Other station-zone sets are in station_zone_summary.csv. This does not endorse an OUTER-only primary analysis: zone definitions and measured depth must be verified, and narrowing zones changes the sampled habitat. No zone-restricted outcome refit is included here.

## Recommendation

Retain the established reference assignments and report these sensitivities. Do not require every historical reference every year. Prioritize a common station/zone prediction population and explicit station adjustment, especially at Point Dume, while checking that seasonal and habitat changes can be separated from recovery. Keep same-year persistent-reference and leave-one-reference checks as diagnostics. Do not select the version with the preferred sign. Before any primary replacement, resolve zone comparability, baseline station continuity and identification, then refit other taxa and pooled/pathway/resilience analyses together.

## Validation and limits

All 34 approved recovery slopes and HAC3 SEs reproduce to <1e-8. R independently reproduces all 68 weighted log-station coefficients to <2e-14. The 162 reported fits include the approved, persistent-reference, matched-year control and two log-station models. All 487 protected manuscript/result files remain unchanged. Year-summed Bartlett HAC3 scores preserve calendar gaps in the exploratory station models; their year-based t reference and finite-sample correction are approximate, not CR2 or validated small-cluster inference. The existing ratio fits retain their original HAC3/t treatment. No regional pooling, before–after refit or claim of causal identification is made. KFM/LTER, satellite, other PISCO taxa and expanded size support remain outside this numerical pass. No emails, commit or push.
'''
(O/'README.md').write_text(text)
cols=['MPA','Resp','Approved','Persistent_reference','Log_station_unadjusted','Log_station_adjusted'];table=wide[wide.MPA.isin(selected+['Abalone Cove SMCA','Santa Barbara Island SMR'])][cols].round(4)
page='<meta charset="utf-8"><title>Station coverage and reference influence</title><style>body{font:16px system-ui;margin:30px;max-width:1250px}img{max-width:100%}td,th{padding:8px;text-align:left;border-bottom:1px solid #ddd}.scroll{overflow:auto}p{max-width:1000px}</style><h1>Station coverage and reference influence</h1><p>PISCO sheephead; same-year diagnostic comparisons. Paper unchanged.</p><p><a href="README.md">Full findings and limitations</a> · <a href="reference_omission.csv">Each reference omitted in turn</a> · <a href="station_zone_summary.csv">Zones shared by every station</a> · <a href="effects.csv">Estimates and intervals</a></p><img src="reference_comparison.png"><p>'+caption+'</p><h2>Comparison estimates</h2><p>Approved versus persistent reference uses log arithmetic role means. The two log-station columns use a different scale; compare them to each other to isolate station adjustment. Missing persistent subsets are unavailable, not zero.</p><div class="scroll">'+table.to_html(index=False,na_rep='Unavailable')+'</div>'
(O/'index.html').write_text(page)
paths=[p for p in O.iterdir() if p.is_file() and p.name!='manifest.json']+[Path('code/tools')/p for p in ['audit_station_coverage_refits.py','verify_station_coverage_refits.R','report_station_coverage_refits.py']];(O/'manifest.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},indent=2))
