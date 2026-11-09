from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,html
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
O=Path('outputs/all_program_station_refits_20261102');e=pd.read_csv(O/'effects.csv');i=pd.read_csv(O/'reference_omission.csv');w=pd.read_csv(O/'comparison_summary.csv');c=pd.read_csv(O/'coverage_summary.csv');a=pd.read_csv('outputs/pisco_fish_rules_20261101/approved_annual_input.csv')
counts=a.groupby('series_id').t.agg(n_post=lambda x:int((x>0).sum()),n_before=lambda x:int((x<0).sum()));c=c.merge(counts,on='series_id',validate='one_to_one');c['recovery_fit_status']=np.where(c.n_post>=3,'Fitted','Fewer than three observed post-protection years');c.to_csv(O/'coverage_summary.csv',index=False)
summary=c.groupby('Source').agg(histories=('series_id','size'),changing_membership_full=('reference_membership_changes','sum'),changing_membership_post=('post_reference_membership_changes','sum'),single_reference=('n_reference_stations',lambda x:int((x==1).sum())),no_persistent_full=('n_persistent_references',lambda x:int((x==0).sum())));summary.to_csv(O/'program_summary.csv')
labels={'M. pyrifera':'kelp','S. pulcher':'sheephead','M. franciscanus':'red urchin','S. purpuratus':'purple urchin','P. interruptus':'lobster'}
# Display the largest omission shifts per program among newly examined responses.
new=i[~(i.Source.eq('PISCO')&i.Taxa.eq('S. pulcher'))].copy();new['abs_change']=new.difference.abs();pick=new.sort_values('abs_change',ascending=False).groupby('Source').head(2).sort_values(['Source','MPA','Taxa','Resp']);pick.to_csv(O/'displayed_reference_omissions.csv',index=False)
fig,ax=plt.subplots(figsize=(10.5,7));ys=[]
for y,q in enumerate(pick.itertuples()):
 ys.append(f'{q.Source}: {q.MPA.replace(" SMCA", "").replace(" SMR", "")}\n{labels[q.Taxa]} {"biomass" if q.Resp=="Bio" else "density"}')
 for arm,dy,col,face in [('before',-.13,'#777777','white'),('after',.13,'#007b79','#007b79')]:
  val=getattr(q,arm);lo=getattr(q,arm+'_low');hi=getattr(q,arm+'_high');ax.errorbar(val,y+dy,xerr=[[val-lo],[hi-val]],fmt='o',mfc=face,color=col,ms=5,capsize=2)
ax.set_yticks(range(len(ys)),ys);ax.axvline(0,color='#aaa',lw=.8);ax.set_xlabel('Annual change in log protected/reference ratio');ax.spines[['top','right']].set_visible(False);fig.tight_layout();fig.savefig(O/'reference_influence.png',dpi=180);fig.savefig(O/'reference_influence.pdf');plt.close(fig)
caption='Gray open points: all available references. Teal filled points: one reference omitted, on identical post-protection years. Bars: 95% calendar-gap HAC3 t intervals. Shown are the two largest absolute omission shifts per program among newly examined responses (previous PISCO sheephead results excluded from display selection). This selection highlights influence, not significant effects or recommended exclusions. Omitted references and exact years are listed in displayed_reference_omissions.csv; all447 checks are in reference_omission.csv. References share observations and these checks are not independent tests.'
(O/'FIGURE_CAPTION.md').write_text(caption+'\n')
text='''# Station coverage and reference checks across programs

November 2, 2026. Review only. Researcher exclusions and the approved 175-history review cohort are preserved; no manuscript, supplement, canonical results or source data were changed.

## Scope and reconstruction

Rebuilt the existing station-level PISCO/KFM/LTER responses in an isolated copy of the processing pipeline, retaining cached biomass calculations and calibrated coefficients. Preserved the approved complete-pass/repeat/CRANE sheephead reconstruction and the separately reconstructed Point Dume membership. Recovered satellite annual values at the delivered polygon level after the documented empty-reference removals. All175 histories reproduce the approved annual ratios to <1e-7 with a single nonnegative additive correction per history. Corrections were recovered algebraically from the existing ratios, held fixed and checked over every year; none was tuned to a desired sensitivity result. This validates the aggregation/reconstruction, not every historical survey convention or size donor.

There are35 KFM,18 LTER,17 satellite and105 PISCO histories. Histories are response/program/MPA combinations, not independent sites. Six PISCO histories lack three observed post-protection years;169 recovery slopes are reproducible. Current counts/size support and inherited seasonal/protocol filters are retained so the comparisons isolate reference handling, rather than mixing new measurement policies.

## What changed

Reference membership changes during the post-protection record in58 PISCO,7 KFM,1 satellite and0 LTER histories. Across the entire before-and-after record the corresponding counts are62,7,3 and2. Requiring a reference to be present before protection as well as after can therefore discard a perfectly continuous recovery-period reference. Both full-record and post-only persistent-reference alternatives are supplied, separately labeled.

PISCO Scorpion kelp is the largest newly identified persistent-reference sensitivity: the annual protected/reference slope changes0.2357→0.0044 on the same20 post-protection years. The intervals are[0.1093,0.3620] and[−0.0790,0.0877]. However, adding station intercepts on the separate log-station scale changes0.1829→0.1824. The narrow reference subset represents different reference locations; its result alone does not establish that changing station membership caused the original trend. The separate KFM Scorpion kelp slope remains0.3884 under these reference checks.

At Campus Point, LTER sheephead reference membership is stable during recovery. Requiring full-record continuity removes Goleta Beach because of baseline coverage and changes biomass slope−0.4624→−0.2441. Requiring only recovery-period continuity leaves−0.4624 unchanged. This is a reference-population choice, not evidence of turnover during recovery.

At Santa Barbara Island, KFM lobster density changes0.1284→0.2049 when restricted to the full-record persistent reference set. Omitting Arch Point in a separate same-year comparison changes0.1392→0.0018 across18 post-years. Those previous slopes differ because the omission comparison is restricted to its own common years. Do not conflate year changes and reference effects.

## Reference identity, even with stable coverage

Removing Mohawk from Campus Point LTER reverses the red-urchin biomass point estimate from−0.1557 to+0.1302 on the same11 years; density changes−0.1235→+0.0117. At Santa Barbara Island KFM, removing Arch Point changes purple-urchin biomass−0.0681→+0.1053 on17 years and density−0.0472→+0.1259 on18 years. These signs describe exploratory point estimates, not new declarations of significance.

Satellite Naples is sensitive to which of its two reference polygons represents outside canopy. Omitting reference polygon1 changes the slope+0.0092→+0.4442 on the same9 post-years. Swamis also differs substantially between references. Polygon totals describe different spatial footprints; these are existing delivered references, not interchangeable diver stations or independent replicates. Empty masks remain unavailable, not zero canopy. No polygon boundary or annual canopy reconstruction is changed.

The source/response table gives every check, including stable cases. No influential reference has been removed and no effect direction determines eligibility.

## Interpretation and recommendation

Keep the researcher-supported assignments. Carry station/polygon identity through analysis and distinguish changing membership from differences in the locations represented. A permanent-reference subset is a useful sensitivity, not automatically a better primary comparison. Existing independent programs can disagree because they sample different footprints; do not merge their station identities or force agreement.

For the main analysis, these results support explicit reference dependence and transparent influence checks. They do not yet select a new primary station model. Log-station models average transformed station values and answer a different question from log ratios of arithmetic means; only their adjusted-minus-unadjusted contrast isolates the station-intercept change. Constant intercepts cannot remove changing habitat, season or station-specific biological trends.

## Biomass, missingness and other limits

Reproducing biomass does not certify sparse or shared size samples. The expanded count-to-size linkage/uncertainty inventory is still required; no new biomass was fabricated for unsupported KFM predators. This check retains existing measurement processing, including unresolved swath visit weighting and source missingness questions, rather than silently fixing them along with reference selection. Different sample sizes for biomass and density remain explicit.

52 histories have only one reference unit, so a reference-omission comparison is unavailable.24 lack any reference present throughout their entire approved record; these remain in the available-reference analysis, not excluded. The post-only persistent subset has broader coverage. No depth/timing rule, before–after estimator or primary cohort was changed. These fits do not establish whether protection caused the observed differences.

## Validation

All169 existing recovery slopes and calendar-gap HAC3 SEs reproduce to <5e-14. R independently reproduces338 log-station coefficients to <1e-12 and447 reference-omission coefficients and interval endpoints to <1e-12. There are950 diagnostic fits plus447 paired omission checks. All487 protected manuscript/result files remain unchanged. The log-station models use equal role/year weights, station/year intercepts and year-summed HAC3 scores; their small-sample intervals remain approximate. No new regional pooled, pathway, resilience or before–after results are claimed.

## Reproduce

Run stage_all_program_station_refits.py; run Rscript run.R inside its ignored stage; run prepare_all_program_stations.py, refit_all_program_stations.py, verify_all_program_stations.R, verify_all_program_omissions.R and report_all_program_stations.py from the analysis root. Source hashes, annual-reconstruction gates, exact model rows, role values, coverage and omissions are saved beside this report. Staging contains source symlinks plus private copies of caches; original source files are not edited. No emails, commit or push.
'''
(O/'README.md').write_text(text)
cols=['Source','MPA','Taxa','Resp','omitted_reference','n_years','before','after'];show=pick[cols].round(4)
page='<meta charset="utf-8"><title>All-program station and reference checks</title><style>body{font:16px system-ui;margin:30px;max-width:1350px}p{max-width:1050px}img{max-width:100%}td,th{padding:8px;text-align:left;border-bottom:1px solid #ddd}.scroll{overflow:auto}</style><h1>Station coverage and reference checks</h1><p>All175 retained review histories reconstructed.950 diagnostic fits and447 same-year reference omissions. Researcher exclusions retained; paper unchanged.</p><p><a href="README.md">Findings and limitations</a> · <a href="coverage_summary.csv">Every history and its reference support</a> · <a href="comparison_summary.csv">All model comparisons</a> · <a href="reference_omission.csv">All reference omissions and intervals</a> · <a href="reconstruction_checks.csv">Annual reconstruction checks</a></p><h2>Reference coverage</h2>'+summary.to_html()+'<p>Full-record and post-protection membership changes are different. Histories are not independent sites.</p><img src="reference_influence.png"><p>'+caption+'</p><div class="scroll">'+show.to_html(index=False)+'</div>'
(O/'index.html').write_text(page)
print(summary.to_string())
