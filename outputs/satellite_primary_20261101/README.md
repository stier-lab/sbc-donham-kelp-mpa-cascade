# Revised satellite primary analysis — review results

November 1, 2026. The four-candidate audit and frozen criteria are complete. The revised satellite primary and trend sensitivity have been fit in this isolated output directory. The manuscript, supplement, diver fits, existing pooled synthesis and prior 16-pair review are unchanged. This is a revised analysis after seeing earlier results, not a preregistered design.

## Audit outcome

| Candidate | Decision | Reason |
|---|---|---|
| Point Conception SMR | Include | Verified no-take reserve and satellite arrays. No demonstrated satellite defect supports the old broad exclusion; the second reference is genuine. |
| Judith Rk SMR | Include | Verified named reserve and supplied arrays. Overlap with excluded San Miguel SC does not duplicate an included treatment. |
| San Miguel Island SC | Exclude | Access/boating closure is not a verified predator-fishing protection treatment; its harmonized SMR type and 2003 clock cannot establish that. |
| Arrow Point to Lion Head Point SMCA | Exclude from this primary scope | Neither focal predator has full ordinary recreational/commercial protection. It does prohibit recreational invertebrate take; this is not a claim of no protection or bad measurements. |

CDFW sources: [Judith/San Miguel](https://wildlife.ca.gov/Conservation/Marine/MPAs/Judith-Rock-San-Miguel-Island), [Point Conception](https://wildlife.ca.gov/Conservation/Marine/MPAs/Point-Conception), [Arrow Point](https://wildlife.ca.gov/Conservation/Marine/MPAs/Arrow-Point-to-Lion-Head-Point). Corrected source arrays passed the existing v 2 verification and reference-correction tests. The full 20-row [eligibility ledger](eligibility_frozen.csv) records coverage, protection scope, source validity and reference limitations.

The scope rule requires a verified reserve/conservation treatment protecting at least one focal predator from both ordinary recreational and commercial take. It keeps Swamis through lobster protection and keeps Naples despite kelp harvest exceptions. No rule requires high predator abundance or excludes a negative, significant or influential baseline trend. Reference habitat equivalence remains unverified across the delivery; newly included sites receive the same evidentiary standard as previously retained ones. Legal boundaries were not independently overlaid for all masks.

## Models and results

The cohort has 18 comparisons and 643 paired annual observations. Each model uses the same observed years within each baseline window. The flat-baseline primary effect is the average Before–After change in log(MPA/reference kelp biomass); the trend sensitivity measures departure from continuation of the earlier trend. Positive values need not mean kelp increased inside the MPA: a falling reference could also increase the ratio. Differences are associations, not automatically causal protection effects.

Full baseline, calendar-gap HAC lag 3, pointwise 95% t intervals:

| Newly included comparison | Flat-baseline primary | Trend sensitivity |
|---|---:|---:|
| Judith Rk SMR | -0.56 (-1.17, 0.05) | -1.13 (-2.90, 0.63) |
| Point Conception SMR | 0.18 (-0.39, 0.76) | 0.47 (-0.56, 1.50) |

Both new comparisons are individually uncertain. Judith is negative and Point Conception slightly positive under both assumptions. Adding them reduces the descriptive equal-site average flat-baseline change from 0.534 to 0.454 log units. Thus their inclusion softens the average positive pattern. This arithmetic mean has no independence-based pooled confidence interval and is not the manuscript's multilevel estimate.

| Cohort, same flat-baseline model | Comparisons | Positive estimates | Equal-site mean | BH-positive comparisons |
|---|---:|---:|---:|---:|
| Previous 14 | 14 | 9 | 0.508 | 4 |
| Previous 16 | 16 | 10 | 0.534 | 5 |
| Frozen 18 | 18 | 11 | 0.454 | 5 |

In 18 comparisons, the flat primary has 11 positive and 7 negative estimates. Five positive departures meet BH q<.05: Abalone Cove, Gull Island, Harris Point, Point Vicente and Santa Barbara Island. No negative primary departure meets that threshold; this is not evidence of no negative effects. The trend sensitivity has 14 positive estimates, with three positive q<.05 departures (Abalone, Harris and Santa Barbara) and one negative (Naples). Point Dume SMR's trend estimate/SE are unchanged from the 16-pair analysis; its q crosses above.05 solely because the adjustment family expands. Do not interpret that crossing as a changed biological effect.

All estimates, standard errors and raw p values for the original 16 are unchanged under the same model/window/lag. The trend model predicts three withheld baseline years better in 8 of 18 comparisons. Baseline-window sensitivity remains substantial; projecting short or declining baselines far forward can produce extreme estimates. Full 18 primary and sensitivity estimates are in [primary_satellite_effects.csv](primary_satellite_effects.csv) and [trend_sensitivity_effects.csv](trend_sensitivity_effects.csv).

## Figures

- [All 18 comparison estimates](all_pair_comparison.png): blue=flat-baseline primary; orange=trend sensitivity. Points show average departures, bars pointwise 95% HAC intervals. Zero means no departure. Sites sharing references are not independent.
- [Newly included site histories](new_pair_trajectories.png): dots=annual observed log MPA/reference ratio; solid lines=fitted Before/After trajectories; dashed lines=continued Before baseline; dotted vertical lines=designation year. Dark ribbons=pointwise 95% fitted-mean intervals; light ribbons=baseline projection intervals, not observation prediction intervals. Separate y scales by site; matched scales between models within sites. Blue flat, orange changing.
- [Baseline-window sensitivity](baseline_windows.png): same colors; full versus final 10 or 5 baseline calendar years, keeping all After observations. Same-year matching applies within each window. Panel-specific y scales.
- All observed trajectories: [1](trajectories_1.png), [2](trajectories_2.png), [3](trajectories_3.png), [4](trajectories_4.png), [5](trajectories_5.png). Same history caption; PDF versions accompany trajectories.

## Reproduction and boundaries

[CRITERIA.md](CRITERIA.md), eligibility_frozen.csv, cohort.csv and annual_input.csv were hashed in freeze_manifest.json before the expanded fits. Both runner and verifier reject a changed frozen input. All 18 fit failures would stop execution rather than silently omit a pair. Source data remain unchanged. This delivery does not replace a historical filter with significance-based model selection.

```sh
Rscript code/tools/compare_satellite_baselines.R outputs/satellite_primary_20261101
python 3 code/tools/check_satellite_baseline_comparison.py outputs/satellite_primary_20261101
python 3 code/tools/summarize_satellite_primary.py
```

Independent NumPy regressions/calendar-HAC calculations reproduce all 1296 contrasts, standard errors and interval endpoints (max absolute error<1.5 e-11). The flat primary equals the direct observed After mean minus Before mean at every site/window. The original 16 effects are unchanged, and scipy independently reproduces all BH adjustments. All 487 protected data/table/plot/document hashes are unchanged. New-site and 18-site figures visually inspected; no embedded titles, prose captions or legends.

The changed primary satellite estimand cannot be silently inserted into the paper's mixture of diver mean-level and first-to-last trend contrasts. legacy_estimand_comparison.csv makes that difference explicit; it is not an apples-to-apples effect-change test. The next integration step is to decide and document the synthesis estimand (or retain separate satellite inference), then rebuild dependent manuscript outputs together. No pooled nine-response conclusion is claimed to have changed here. No paper edits, sends, commits or push.
