#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
out=Path('outputs/pooled_common_estimand_20261101');r=pd.read_csv(out/'pooled.csv');e=pd.read_csv(out/'effects.csv');c=pd.read_csv(out/'coverage.csv');old=pd.read_csv('tables/table_02_meta_analysis.csv')
names={'M. pyrifera':'Kelp','S. pulcher':'Sheephead','P. interruptus':'Lobster','M. franciscanus':'Red urchin','S. purpuratus':'Purple urchin'}
primary=r.query("scenario=='Primary' and view=='All'")
primary.to_csv(out/'primary_pooled.csv',index=False)
old['Resp']=old.Response.map({'Biomass':'Bio','Density':'Den'})
comparison=old[['Taxa','Resp','Estimate','CI_lower','CI_upper','pval_fdr']].merge(primary,on=['Taxa','Resp'],how='left');comparison.to_csv(out/'old_vs_new.csv',index=False)
text='''# Satellite and diver comparisons on consistent scales

November 1, 2026. Review analysis only. The manuscript, supplement, source data and canonical primary outputs remain unchanged.

## Recommended structure

Separate two pooled questions. They cannot be collapsed into one protection-effect coefficient:

1. **Before–After change:** how much did the protected/reference ratio change after designation? Require usable baseline data and estimate the same mean change for every program. This is the primary analysis for change associated with designation, with reference-site and counterfactual qualifications.
2. **After-protection difference:** how different were protected and reference areas during the post-designation period? This retains the broader monitoring coverage, but it cannot attribute existing differences to protection.

Do not convert a baseline-poor record into an estimated protection effect by fitting its later slope. A post-designation trend answers a third question and can remain a separate supporting analysis.

Primary timing is strictly prior Before observations and observed years 1–9 after designation. The nine-year horizon is available in the 2012 satellite cohort through 2021; it avoids comparing nine-year and eighteen-year responses as the same duration. Omit designation year from both periods for all programs, while retaining original source labels unchanged. Earlier annual-label conventions, all observed After years, the last ten baseline years, two alternate HAC lags and changing-baseline models are explicit sensitivities. Missing years are not interpolated. Observed-year mixtures still differ within the common horizon; exact same-year source comparisons are provided below.

Scope is the existing 129 retained diver histories plus the independently audited 18 satellite histories. Historical omitted diver histories were not silently reinstated. The six fixed-station PISCO sheephead histories have verified branch-level baseline eligibility; other retained PISCO histories contribute to post differences, not newly asserted baseline designs. The existing KFM/LTER assignments are retained, with their documented station-continuity limits. This resolves effect-definition mixing within the analyzed cohort; it is not a completed outcome-independent audit of every past diver exclusion.

## Coverage

Of 147 candidate histories, 45 support the common Before–After change: 27 diver and 18 satellite histories. Another 102 do not meet that design/window requirement: 94 have fewer than three strictly prior observations, six lack a verified PISCO fixed-baseline design, and two have fewer than three post years in the common horizon. Histories without baseline are not discarded from the study.

The post-difference analysis has 145 histories. Farnsworth PISCO sheephead biomass and density each have only two observations in years 1–9; they return in the full-observed-post sensitivity. Lobster biomass and purple-urchin biomass have no eligible Before–After histories in this retained cohort. Purple-urchin density has just one; there is no pooled change inference for it.

## What changes

Kelp's old primary estimate was 1.138 (95% CI 0.613–1.662), but it mixed mean differences with fitted changes. The new common Before–After kelp estimate is **0.309 (−0.159–0.778)** from 23 program histories at 19 MPAs. The new post-protection difference is **0.326 (−0.141–0.792)** from 35 histories at 21 MPAs. These are new, clearly defined questions; the old and new numbers are not alternative estimates of an identical quantity.

A staged check separates the main causes:

| Kelp comparison | Histories | Estimate | 95% interval |
|---|---:|---:|---:|
| Original effects, new per-response pooling and shared-reference robust uncertainty | 31 | 1.087 | 0.673 to 1.501 |
| Original effects restricted to histories eligible for the new change analysis | 19 | 1.146 | 0.613 to 1.679 |
| Same 19 histories, common mean-change definition and time window | 19 | 0.470 | −0.057 to 0.997 |
| Expanded cohort, common mean change | 23 | 0.309 | −0.159 to 0.778 |

Thus the result is not explained simply by robust uncertainty or adding two sites. Resolving the mixed questions and time support makes the largest difference in this sequence; adding the four satellite histories absent from the current 14-satellite production cohort reduces the estimate further. Sequential attribution depends on the sequence and is not an independent decomposition of causal contributions.

The common post-protection analysis supports higher sheephead biomass, sheephead density and red-urchin biomass after adjustment across nine responses (q≈0.025 each). Kelp and purple-urchin density no longer have clear pooled differences in this common window. Positive post differences do not establish increases caused by protection. None of the Before–After coefficients with usable robust inference survives that adjustment, while several other responses have too few independent baseline comparisons for reliable inference. This is not evidence that their true effects are zero.

## Primary pooled estimates

Values are natural-log ratio units. Change is a difference in log ratios; post difference is a log ratio. Intervals use reference-group/MPA CR2 inference, not the old model-based errors. Missing intervals mean insufficient independent information, never zero uncertainty.

| Response | Quantity | Histories / MPAs | Estimate | 95% interval | Adjusted p |
|---|---|---:|---:|---|---:|
'''
for _,g in primary.iterrows():
 label=names[g.Taxa]+' '+('biomass' if g.Resp=='Bio' else 'density');metric='Before–After change' if g.metric=='Change' else 'Post difference';ci=f'{g.low:.2f} to {g.high:.2f}' if pd.notna(g.low) else 'Insufficient independent information';q=f'{g.q:.3f}' if pd.notna(g.q) else '—'
 text+=f'| {label} | {metric} | {g.k} / {g.n_mpa} | {g.estimate:.2f} | {ci} | {q} |\n'
text+='''
## How sensitive is the kelp conclusion?

| Specification | Estimate | 95% interval | Adjusted p |
|---|---:|---|---:|
'''
for scenario in ['Primary','Full_after','Baseline10','Original_timing','HAC1','HAC5','Trend']:
 g=r.query("scenario==@scenario and metric=='Change' and view=='All' and Taxa=='M. pyrifera'").iloc[0];text+=f'| {scenario} | {g.estimate:.3f} | {g.low:.3f} to {g.high:.3f} | {g.q:.3f} |\n'
text+='''
Using the full observed post history gives a larger kelp change, with a pointwise interval above zero, but q=0.128 across the nine-response family. Using the last ten Before years also gives a pointwise positive interval, but q=0.410. The common-window result is similar under alternate annual timing and HAC lags. Trend adjustment remains uncertain. Do not choose the horizon or baseline because it produces the desired interval.

The common-window satellite-only change is 0.271 (−0.208–0.749). Diver-only change is approximately 0.984, but the five available diver histories provide fewer than two effective CR2 degrees of freedom; no robust inferential interval is reported. The canopy and diver measurements remain different biological measurements even after standardizing their contrast. Keep these source-separated estimates visible.

## Exact same-year kelp comparison

These four MPAs have both eligible satellite and diver kelp records. The following means use exactly the same observed Before and After calendar years within each MPA. This removes time coverage as the explanation for differences within each row, but does not equalize footprints or canopy versus diver measurements. Values are descriptive mean changes, not independent pooled estimates.

| MPA | Diver program | Diver change | Satellite change |
|---|---|---:|---:|
'''
m=pd.read_csv(out/'matched_kelp_years.csv')
for mpa,g in m.groupby('MPA'):
 sat=g[g.Source=='Landsat'].iloc[0];div=g[g.Source!='Landsat'].iloc[0]
 for field in ['before_years','after_years']:assert set(sat[field].split(';'))==set(div[field].split(';'))
 text+=f'| {mpa} | {div.Source} | {div.change:.3f} | {sat.change:.3f} |\n'
text+='''
Naples remains the clearest discrepancy: diver kelp increases relative to its reference while satellite canopy declines relative to its own reference, even on identical years. Pooling alone cannot resolve a footprint, reference or measurement difference. Scorpion contributes to the diver-only baseline set but lacks an included satellite comparison.

## Estimation and dependence

Within-history effects come from flat-Before/linear-After models, with calendar-gap HAC covariance. Their average contrasts equal direct observed period means; the trend sensitivity allows a separate earlier slope. Fit each taxon and biomass/density outcome separately using REML. Random MPA intercepts account for repeated programs at an MPA; an additional series heterogeneity term is used only where those repeated measurements exist. This avoids borrowing the same variance across different taxa merely because they occur in one model.

Shared MPA and known shared-reference connections form CR2 clusters, including Naples/Campus satellite pixel overlap and retained diver reference identities. Satterthwaite degrees of freedom govern intervals. Fewer than three groups or df<2 yields no inferential interval/p; df<4 is flagged as fragile. Model-based intervals and variance components are exported as diagnostics, not substituted when robust inference is weak. Some fitted variance components approach zero; see fit_log.csv and pooled.csv rather than assuming all random components are identified. Known reference clustering does not eliminate unmeasured regional dependence, and annual HAC errors are not a rerun of the full length-to-biomass uncertainty propagation.

## Figure and outputs

[Comparison figure](common_estimands.png): colored points are individual program-history estimates (satellite orange, KFM blue, LTER green, PISCO purple). Filled black diamonds and bars show pooled estimates and pointwise 95% CR2 intervals. Hollow diamonds denote descriptive pooled point estimates with insufficient information for an inferential interval. A lone history has no pooled diamond. Blank rows in the change panel mean no eligible Before–After data, not no effect. The left and right panels answer different questions. Colors and explanations are deliberately outside the image.

- primary_pooled.csv: primary results.
- pooled.csv: all 238 primary/sensitivity/source-specific result rows.
- coverage.csv: explicit eligibility and reason for every history, scenario and question.
- effects.csv and model_rows.csv: 1,333 history-level fits and their exact annual observations.
- legacy_bridge.csv and kelp_same_cohort_bridge.csv: staged comparison.
- matched_kelp_years.csv: exact source-matched years and descriptive means.
- reproduced_current_pooled.csv: independent reproduction of the existing nine-response Table 2.
- PLAN.md and freeze_manifest.json: decisions and inputs fixed before fits.

## Reproduce and checks

```sh
Rscript code/tools/refit_common_estimand.R
Rscript code/tools/bridge_common_estimand.R
python3 code/tools/check_common_estimand.py
python3 code/tools/summarize_common_estimand.py
```

All 1,333 within-history estimates/HAC standard errors independently reproduce in NumPy (maximum error <7e-13). Pooled GLS coefficients independently reconstruct from exported variance components (<6e-15 maximum error); CR2 interval arithmetic and BH adjustments also pass. The existing nine Table 2 estimates and adjusted p values reproduce. Tests require unique scenario/quantity/history keys, and all 487 protected table/plot/input/document fingerprints remain unchanged. The final comparison figure was visually checked. Initial implementation issues were corrected before accepting these results; only the passing final run is reported here.

Recommendation for the paper: separate post-protection differences from Before–After evidence, use one quantity per pool, retain source-separated kelp results, and describe the baseline limitations plainly. The present analysis does not support carrying forward the old pooled coefficient as a general protection effect. Before promotion, resolve whether to retain this nine-year target or the full observed-history target on biological grounds, audit the remaining historical diver omissions under the same outcome-independent logic, and update affected size-uncertainty/pathway/resilience dependencies together. No main-text/supplement edits or full downstream rebuild has been performed in this review-only turn.
'''
(out/'README.md').write_text(text)
print('Saved report and old/new comparison; exact-source calendar matches checked')
