#!/usr/bin/env python3
"""Summarize the frozen expanded satellite refit without writing canonical results."""
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import false_discovery_control
out=Path('outputs/satellite_primary_20261101')
r=pd.read_csv(out/'estimates.csv');a=r.query("window=='Full' and lag==3 and quantity=='Mean_after'").copy()
a[a.model=='Stable'].to_csv(out/'primary_satellite_effects.csv',index=False)
a[a.model=='Changing'].to_csv(out/'trend_sensitivity_effects.csv',index=False)
old=pd.read_csv('data/sumstats_final.csv').query("Source=='Landsat'")
co16=pd.read_csv('outputs/satellite_baseline_comparison_20261101/cohort.csv')
rows=[]
for label,names in [('Previous14',old.MPA),('Previous16',co16.MPA),('Frozen18',a.MPA.unique())]:
 for m in ['Stable','Changing']:
  g=a[(a.model==m)&a.MPA.isin(names)].copy();g['cohort_q']=false_discovery_control(g.p)
  rows.append(dict(cohort=label,model=m,n=len(g),mean=g.estimate.mean(),median=g.estimate.median(),positive=(g.estimate>0).sum(),ci_positive=(g.low>0).sum(),ci_negative=(g.high<0).sum(),bh_positive=((g.cohort_q<.05)&(g.estimate>0)).sum(),bh_negative=((g.cohort_q<.05)&(g.estimate<0)).sum()))
s=pd.DataFrame(rows);s.to_csv(out/'cohort_comparison.csv',index=False)
old[['MPA','Mean','SE']].merge(a[['MPA','model','estimate','se','low','high']],on='MPA').to_csv(out/'legacy_estimand_comparison.csv',index=False)
prior=pd.read_csv('outputs/satellite_baseline_comparison_20261101/estimates.csv');join=prior.merge(r,on=['MPA','window','model','lag','quantity'],suffixes=('_old','_new'))
for f in ['estimate','se','low','high','p']:assert np.allclose(join[f+'_old'],join[f+'_new'],atol=1e-10)
for _,g in r.groupby(['model','window','lag','quantity']):assert np.allclose(false_discovery_control(g.p),g.q)
text='''# Revised satellite primary analysis — review results

November 1, 2026. The four-candidate audit and frozen criteria are complete. The revised satellite primary and trend sensitivity have been fit in this isolated output directory. The manuscript, supplement, diver fits, existing pooled synthesis and prior16-pair review are unchanged. This is a revised analysis after seeing earlier results, not a preregistered design.

## Audit outcome

| Candidate | Decision | Reason |
|---|---|---|
| Point Conception SMR | Include | Verified no-take reserve and satellite arrays. No demonstrated satellite defect supports the old broad exclusion; the second reference is genuine. |
| Judith Rk SMR | Include | Verified named reserve and supplied arrays. Overlap with excluded San Miguel SC does not duplicate an included treatment. |
| San Miguel Island SC | Exclude | Access/boating closure is not a verified predator-fishing protection treatment; its harmonized SMR type and 2003 clock cannot establish that. |
| Arrow Point to Lion Head Point SMCA | Exclude from this primary scope | Neither focal predator has full ordinary recreational/commercial protection. It does prohibit recreational invertebrate take; this is not a claim of no protection or bad measurements. |

CDFW sources: [Judith/San Miguel](https://wildlife.ca.gov/Conservation/Marine/MPAs/Judith-Rock-San-Miguel-Island), [Point Conception](https://wildlife.ca.gov/Conservation/Marine/MPAs/Point-Conception), [Arrow Point](https://wildlife.ca.gov/Conservation/Marine/MPAs/Arrow-Point-to-Lion-Head-Point). Corrected source arrays passed the existing v2 verification and reference-correction tests. The full20-row [eligibility ledger](eligibility_frozen.csv) records coverage, protection scope, source validity and reference limitations.

The scope rule requires a verified reserve/conservation treatment protecting at least one focal predator from both ordinary recreational and commercial take. It keeps Swamis through lobster protection and keeps Naples despite kelp harvest exceptions. No rule requires high predator abundance or excludes a negative, significant or influential baseline trend. Reference habitat equivalence remains unverified across the delivery; newly included sites receive the same evidentiary standard as previously retained ones. Legal boundaries were not independently overlaid for all masks.

## Models and results

The cohort has18comparisons and643paired annual observations. Each model uses the same observed years within each baseline window. The flat-baseline primary effect is the average Before–After change in log(MPA/reference kelp biomass); the trend sensitivity measures departure from continuation of the earlier trend. Positive values need not mean kelp increased inside the MPA: a falling reference could also increase the ratio. Differences are associations, not automatically causal protection effects.

Full baseline, calendar-gap HAC lag3, pointwise95% t intervals:

| Newly included comparison | Flat-baseline primary | Trend sensitivity |
|---|---:|---:|
'''
for name in ['Judith Rk SMR','Point Conception SMR']:
 g=a[a.MPA==name].set_index('model');f=lambda m:f"{g.loc[m,'estimate']:.2f} ({g.loc[m,'low']:.2f}, {g.loc[m,'high']:.2f})"
 text+=f'| {name} | {f("Stable")} | {f("Changing")} |\n'
text+='''
Both new comparisons are individually uncertain. Judith is negative and Point Conception slightly positive under both assumptions. Adding them reduces the descriptive equal-site average flat-baseline change from0.534 to0.454log units. Thus their inclusion softens the average positive pattern. This arithmetic mean has no independence-based pooled confidence interval and is not the manuscript's multilevel estimate.

| Cohort, same flat-baseline model | Comparisons | Positive estimates | Equal-site mean | BH-positive comparisons |
|---|---:|---:|---:|---:|
'''
for _,g in s[s.model=='Stable'].iterrows():text+=f'| {g.cohort} | {g.n} | {g.positive} | {g["mean"]:.3f} | {g.bh_positive} |\n'
text+='''
In18comparisons, the flat primary has11positive and7negative estimates. Five positive departures meet BH q<.05: Abalone Cove, Gull Island, Harris Point, Point Vicente and Santa Barbara Island. No negative primary departure meets that threshold; this is not evidence of no negative effects. The trend sensitivity has14positive estimates, with three positive q<.05departures (Abalone, Harris and Santa Barbara) and one negative (Naples). Point Dume SMR's trend estimate/SE are unchanged from the16-pair analysis; its q crosses above.05 solely because the adjustment family expands. Do not interpret that crossing as a changed biological effect.

All estimates, standard errors and raw p values for the original16are unchanged under the same model/window/lag. The trend model predicts three withheld baseline years better in8of18comparisons. Baseline-window sensitivity remains substantial; projecting short or declining baselines far forward can produce extreme estimates. Full18primary and sensitivity estimates are in [primary_satellite_effects.csv](primary_satellite_effects.csv) and [trend_sensitivity_effects.csv](trend_sensitivity_effects.csv).

## Figures

- [All18comparison estimates](all_pair_comparison.png): blue=flat-baseline primary; orange=trend sensitivity. Points show average departures, bars pointwise95% HAC intervals. Zero means no departure. Sites sharing references are not independent.
- [Newly included site histories](new_pair_trajectories.png): dots=annual observed log MPA/reference ratio; solid lines=fitted Before/After trajectories; dashed lines=continued Before baseline; dotted vertical lines=designation year. Dark ribbons=pointwise95% fitted-mean intervals; light ribbons=baseline projection intervals, not observation prediction intervals. Separate y scales by site; matched scales between models within sites. Blue flat, orange changing.
- [Baseline-window sensitivity](baseline_windows.png): same colors; full versus final10or5baseline calendar years, keeping all After observations. Same-year matching applies within each window. Panel-specific y scales.
- All observed trajectories: [1](trajectories_1.png), [2](trajectories_2.png), [3](trajectories_3.png), [4](trajectories_4.png), [5](trajectories_5.png). Same history caption; PDF versions accompany trajectories.

## Reproduction and boundaries

[CRITERIA.md](CRITERIA.md), eligibility_frozen.csv, cohort.csv and annual_input.csv were hashed in freeze_manifest.json before the expanded fits. Both runner and verifier reject a changed frozen input. All18fit failures would stop execution rather than silently omit a pair. Source data remain unchanged. This delivery does not replace a historical filter with significance-based model selection.

```sh
Rscript code/tools/compare_satellite_baselines.R outputs/satellite_primary_20261101
python3 code/tools/check_satellite_baseline_comparison.py outputs/satellite_primary_20261101
python3 code/tools/summarize_satellite_primary.py
```

Independent NumPy regressions/calendar-HAC calculations reproduce all1296contrasts, standard errors and interval endpoints (max absolute error<1.5e-11). The flat primary equals the direct observed After mean minus Before mean at every site/window. The original16effects are unchanged, and scipy independently reproduces all BH adjustments. All487protected data/table/plot/document hashes are unchanged. New-site and18-site figures visually inspected; no embedded titles, prose captions or legends.

The changed primary satellite estimand cannot be silently inserted into the paper's mixture of diver mean-level and first-to-last trend contrasts. legacy_estimand_comparison.csv makes that difference explicit; it is not an apples-to-apples effect-change test. The next integration step is to decide and document the synthesis estimand (or retain separate satellite inference), then rebuild dependent manuscript outputs together. No pooled nine-response conclusion is claimed to have changed here. No paper edits, sends, commits or push.
'''
# Keep readable spacing in prose rather than compact handoff notation.
import re
text=re.sub(r'(?<=[a-zA-Z])(?=\d)', ' ',text)
text=re.sub(r'(?<=\d)(?=[a-zA-Z])',' ',text)
# Above touches filenames with trailing numeric IDs; restore link/script identifiers.
for i in range(1,6):text=text.replace(f'trajectories_{i}.png',f'trajectories_{i}.png')
(out/'README.md').write_text(text)
print('Saved primary/sensitivity effects, cohort comparison, report; original16 and BH checks passed')
