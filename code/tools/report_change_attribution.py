from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path('outputs/change_attribution_20261101');d=pd.read_csv(p/'pooled.csv');old=pd.read_csv(p/'A0_current.csv').rename(columns={'ci_low':'low','ci_high':'high'})
focal=[('S. purpuratus','Den','Purple urchin density'),('S. purpuratus','Bio','Purple urchin biomass'),('M. franciscanus','Den','Red urchin density'),('M. franciscanus','Bio','Red urchin biomass'),('P. interruptus','Den','Lobster density'),('P. interruptus','Bio','Lobster biomass')]
def val(frame,stage,metric,tx,r):
 x=frame[(frame.stage==stage)&(frame.metric==metric)&(frame.Taxa==tx)&(frame.Resp==r)]
 return None if x.empty else x.iloc[0]
def cell(x):
 if x is None:return 'Not estimable'
 v=f'{x.estimate:.3f}'
 return v+(f' [{x.low:.3f}, {x.high:.3f}]; q={x.q:.3f}' if np.isfinite(x.low) else '; no reliable pooled interval')
def table(stages,metric,frame=d):
 lines=['| Response | '+' | '.join(l for s,l in stages)+' |','|---|'+'---|'*len(stages)]
 for tx,r,l in focal:lines.append('| '+l+' | '+' | '.join(cell(val(frame,s,metric,tx,r)) for s,label in stages)+' |')
 return '\n'.join(lines)
sets=[('pooling_uncertainty','Legacy',[('A0_current','Current manuscript'),('A1_pool_structure','Separate-response pooling'),('A2_cluster_uncertainty','Cluster-robust intervals'),('A3_annual_HAC','Annual HAC covariance')]),('definition_window','Post_difference',[('B0_matched_legacy','Legacy contrast'),('B1_common_definition_same_observations','Mean After ratio'),('B2_calendar_designation_rule','Strict calendar periods'),('B3_nine_year_window','After years 1–9'),('B4_target_observations','After observations only')]),('membership','Post_difference',[('M0_old_roster','Old eligible roster'),('M1_removals_only','Removals only'),('M2_additions_only','Additions only'),('M3_final_roster','Final roster')])]
allframe=pd.concat([old,d],ignore_index=True)
for name,metric,stages in sets:
 fig,axs=plt.subplots(3,2,figsize=(15,10));fig.subplots_adjust(left=.19,right=.98,bottom=.08,top=.97,hspace=.5,wspace=.95)
 for ax,(tx,r,label) in zip(axs.flat,focal):
  for j,(s,_) in enumerate(stages):
   x=val(allframe,s,metric,tx,r)
   if x is None:continue
   if np.isfinite(x.low):ax.plot([x.low,x.high],[j,j],color='#286785',lw=1.5)
   ax.plot(x.estimate,j,'o',color='#286785',mfc='#286785' if np.isfinite(x.low) else 'white',ms=6)
  ax.axvline(0,color='.6',lw=.8,ls='--');ax.set_yticks(range(len(stages)),[label for s,label in stages],fontsize=9);ax.invert_yaxis();ax.set_ylabel(label,fontsize=10);ax.spines[['top','right']].set_visible(False);ax.grid(axis='x',alpha=.15)
  ax.set_xlabel('Log contrast' if metric=='Legacy' else 'Log protected/reference ratio',fontsize=10)
 fig.savefig(p/(name+'.png'),dpi=170);fig.savefig(p/(name+'.pdf'));plt.close(fig)
text='''# What changes the urchin and lobster results?

November 1, 2026. Review only; the manuscript, supplement and canonical outputs remain unchanged. No eligibility tolerance selected.

## Main finding

The changes have different causes for different responses. Purple-urchin density and lobster density lose support when the legacy mixture of means, steps and fitted trends is replaced with an average post-protection difference and a common nine-year window. They do not lose support merely because cluster-robust or annual-correlation uncertainty is added. Red-urchin biomass is especially sensitive to the pooling structure and the added histories. Purple-urchin biomass also loses support under a consistent post-protection target before new sites are added. Lobster biomass and red-urchin density remain uncertain throughout the comparisons below.

These are controlled, sequential and factorial comparisons, not a universal percentage decomposition. Changes in pooling weights and estimated heterogeneity interact with membership and covariance. BH-adjusted p values depend on the full nine-response family, so q can change even when a response's estimate and raw p are unchanged. Use the point estimates and raw intervals alongside q.

All tables below report log-scale estimates, 95% intervals and BH q; zero denotes no difference/change. Legacy contrasts and post-protection differences answer different questions despite sharing a scale.

## 1. Pooling and uncertainty, exactly the same 143 legacy histories

The original joint model shares heterogeneity across taxa within biomass/density models and has MPA and program random effects. Separate-response pooling estimates heterogeneity for each response, with MPA and repeated series components. The next step changes only inference to CR2 with shared-reference clusters. The last changes only individual covariance to calendar HAC3; point estimates remain fixed individually, but their weights can change the pooled estimate.

'''+table(sets[0][2],'Legacy',allframe)+'''

Red-urchin biomass changes substantially at the pooling-structure step, before any histories or annual observations are changed. Purple-urchin density remains supported after both uncertainty changes. Lobster density and purple-urchin biomass actually receive stronger adjusted support under these alternative pooling specifications. Robust uncertainty is not uniformly wider or more conservative.

`A0b_current_pool_HAC_floor.csv` separately reproduces the manuscript's annual-correlation convention: original joint pooling and max(original SE, HAC SE). This differs from the proposed raw-HAC convention; do not conflate them.

## 2. Effect definition and follow-up window, fixed intersecting roster

First retain only legacy histories estimable under the proposed target, then compare models without any further membership changes. For post differences this is 137 histories. The first two columns use identical annual observations for every history. Common models use a flat Before level and a post-period slope (or a post-period line when no Before observations exist). The target is the average fitted After ratio. Designation-year removal/calendar classification, nine-year restriction, and removal of unused Before observations then occur in distinct steps. All steps use the same separate-response pooling, HAC3 covariance and CR2 inference.

'''+table(sets[1][2],'Post_difference')+'''

For purple-urchin density, changing the target on the same observations shifts the log contrast from −1.379 to −0.902; restricting follow-up shifts it further to −0.433. Adding annual-correlation uncertainty alone did not erase the old finding. For lobster density the corresponding sequence is 1.161 → 1.010 → 0.737. Both intervals include zero after standardizing the follow-up window.

## 3. Uncertainty alone on the same new-model observations

All cells below use the same 137-history roster, first-nine-year observations and effect definitions. IID versus HAC changes only individual covariance; model versus CR2 changes only pooled inference. CR2 leaves the fitted pooled point estimate exactly unchanged for any fixed covariance input.

'''+table([('U1_same_data_IID_model','IID / model interval'),('U2_same_data_IID_CR2','IID / CR2'),('U3_same_data_HAC_model','HAC / model interval'),('U4_same_data_HAC_CR2','HAC / CR2')],'Post_difference')+'''

Purple-urchin density and lobster density are already uncertain using conventional IID/model intervals under the standardized target and window. Thus loss of support cannot be attributed simply to the robust uncertainty method. Red-urchin biomass is borderline across uncertainty choices on this intersecting roster.

## 4. Membership alone: removals and additions separated

Every retained history has exactly the same observations, individual effect, covariance, pooling specification and inference across these four columns. The old roster contains 141 post histories after the fixed nine-year coverage requirement removes the two Farnsworth sheephead histories. Researcher removals then remove four Gull KFM histories. Additions supply 25 histories. The final roster has 162 post histories. No baseline-trend tolerance screen is applied.

'''+table(sets[2][2],'Post_difference')+'''

Added histories substantially reduce the red-urchin biomass contrast (0.708 → 0.245 on the old roster plus additions). Removing Gull KFM alone changes it less (0.708 → 0.607). For purple-urchin density, additions move the estimate toward zero (−0.430 → −0.262), but it was already uncertain without them. No new lobster-density histories are added: its numerical change (0.684 → 0.737) comes from removing Gull KFM, not additions. Changes in q across other membership arms can reflect the nine-response BH adjustment.

## 5. Before–after changes: membership alone

These are distinct from the post-protection comparisons above. The old eligible roster supports 41 before–after histories; researcher removals leave37, additions alone yield58, and both together yield54. Several responses have too few independent clusters/effective degrees of freedom for reliable intervals. No before–after lobster biomass history is eligible.

'''+table(sets[2][2],'Change')+'''

Purple-urchin density had only one eligible legacy before–after history. Adding histories provides broader coverage but not a precise pooled decline. Purple-urchin biomass had no eligible legacy before–after histories; its four-history estimate is new, with inadequate independent information for reliable pooled inference. Removing Gull KFM reduces already sparse red-urchin biomass and lobster-density before–after replication. Missing robust intervals are not evidence of no biological effect.

## Source reconstruction check

C1–C3 in pooled.csv separate common-year input values from reconstructed survey support. For the six focal urchin/lobster responses, reconstructed source values/support do not change individual contrasts on the retained intersecting roster. KFM kelp reconstruction affects the kelp result and can indirectly change family-adjusted q. Source changes were not silently attributed to site additions.

## Checks, limits and reproduction

Current Table2 and the proposed endpoint were independently matched to saved outputs. The R reconstruction checks all143 legacy effects and SEs and requires convergence of accepted nonlinear fits. Three warnings arose from unsuccessful nonlinear starting fits; accepted final fits converged and reproduced the saved estimates. New common estimands are full-rank linear fits. This is conditional on existing legacy model selection, not a new model-selection exercise or proof of causal validity.

Independent Python checks reconstruct common individual contrasts/HAC covariances, pooled GLS point estimates, interval arithmetic and nine-response BH correction. They assert unchanged cohort membership through the B model steps, unchanged per-history values across membership arms, and original observation counts for the same-observation estimator comparison. Existing protected-file hashes are checked. See verification.json and endpoint_verification.csv. Sampling/size uncertainty and causal identification limitations remain; this comparison does not resolve them.

```sh
Rscript code/tools/attribute_analysis_changes.R
python3 code/tools/check_change_attribution.py
python3 code/tools/report_change_attribution.py
```

## Figure captions

- pooling_uncertainty.png/pdf: same143 legacy histories; points and95%intervals after each pooling/uncertainty step. The current joint model is followed by separate-response pooling, CR2 inference and calendarHAC3 individual covariance. The axis shows the legacy mixed contrast; zero is no contrast.
- definition_window.png/pdf: fixed intersecting legacy roster. Effect definition, designation-year treatment, follow-up window and removal of Before observations change sequentially. Points/intervals use the same pooling and robust covariance convention; their biological interpretation changes from the legacy mixture to mean post-protection differences.
- membership.png/pdf: fixed new estimators and individual inputs, with old eligible roster, removals alone, additions alone and both. Points are log protected/reference ratios; bars are95%CR2 intervals. Open points would indicate unsupported pooled intervals. No embedded legends, titles or caption prose.

Machine-readable effects.csv, observations.csv and pooled.csv retain all nine responses and all comparison steps. The paper has not been updated.
'''
(p/'README.md').write_text(text)
print('Reports and three figures written')
