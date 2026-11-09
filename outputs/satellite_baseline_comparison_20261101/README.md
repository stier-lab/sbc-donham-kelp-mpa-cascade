# Satellite baseline comparison: results for review

November 1, 2026. The manuscript, supplement and primary results have not been changed.

All 16 requested comparisons have enough paired observations for these fits: 569 annual observations, 15–28 Before years and 9–18 After years per pair. Missing years were left missing. Eligibility is provisional for this comparison: the table documents protection level, coverage, supplied reference identity, current inclusion and unresolved ecological matching. It does not certify every reference as a valid counterfactual. Naples and Campus Point share reference pixels. The four other harmonized satellite candidates were outside this requested analysis.

## What the comparison asks

Did the kelp biomass ratio between the MPA and its reference depart from its earlier pattern after designation? Blue assumes a flat earlier baseline; orange allows a linear earlier trend. Both models allow a level change and a post-designation slope, use exactly the same years within each baseline window, and fit the same post-designation line. Their difference is the baseline against which that line is assessed.

The reported quantity is the average post-designation departure from continuation of the earlier baseline, in natural-log ratio units. Positive means a higher MPA/reference ratio than projected; it does not necessarily mean that kelp increased inside the MPA. Negative means lower than projected, not necessarily less kelp inside than outside. These quantities are not the existing primary pBACIPS contrasts and must not be substituted into its pooled analysis without an explicit estimand decision. The current primary model also lacks the level-change term included in both diagnostic models here.

## What changes

With the full baseline, 10 of 16 flat-baseline estimates and 13 of 16 changing-baseline estimates are positive. At lag-3 HAC uncertainty and after adjustment within each 16-site family, five flat-baseline positive departures and four changing-baseline positive departures have q<.05; Naples has a negative changing-baseline departure. These counts describe this dataset, not independent replications or a network-level causal test.

- Santa Barbara Island remains positive under both assumptions. Its estimate increases from 1.97 to 6.19 log units when the previous decline is projected forward. The very large adjusted value depends on extending that decline for another 18 years; it should not be presented as a reliable magnitude of protection benefit.
- Swamis changes from −0.53 to +0.81 with full-baseline adjustment, but both intervals span zero. The adjusted estimate becomes negative with the last 10 or 5 baseline years. Its direction is not robust.
- Naples changes from −0.07 to −1.02. Its earlier MPA/reference ratio was rising; the later ratio falls below that projection. The adjusted estimate remains negative across windows, but its interval spans zero with only five baseline years.
- Gull Island's positive departure largely disappears after allowing its earlier increase. Point Vicente's estimate barely changes but uncertainty widens. Harris Point becomes more positive because its earlier relative decline is projected forward.
- Cabrillo and Farnsworth also change sign under full-baseline adjustment, but their adjusted intervals span zero.

## Checks that qualify interpretation

The changing-baseline model predicts the last three withheld baseline observations better at only 6 of 16 sites. That limited holdout checks near-term prediction, not an 18-year extrapolation. Using only the final 10 baseline years changes the number of positive trend-adjusted estimates from 13 to 6; using five gives 8. Five-year trend extrapolations can produce extreme estimates (for example, Carrington and Cabrillo). Shortening the baseline is not a generally reliable remedy.

Changing the calendar-HAC lag from three to one reduces the number of positive trend-adjusted q<.05 departures from four to two; Naples remains negative. Exported intervals describe uncertainty conditional on each model. They do not include uncertainty about whether a historical trend would continue, reference suitability or spatial dependence. HAC precision with a short baseline is itself uncertain.

Recommendation: keep data/design eligibility separate from baseline-model suitability. These checks do not justify excluding a series because its slope is negative or its baseline test is significant. They also do not justify a blanket replacement of flat baselines with extrapolated linear trends. Show both assumptions and identify conclusions that depend on them before deciding which satellite estimand belongs in the paper. No model was selected by which result favored protection.

## Full-baseline departures

Values are estimates and pointwise 95% calendar-HAC intervals (lag 3); multiplicity-adjusted q values are in estimates.csv.

| Comparison | Flat baseline | Changing baseline |
|---|---:|---:|
| Abalone Cove SMCA | 1.54 (0.79, 2.29) | 1.46 (0.37, 2.55) |
| Cabrillo SMR | -0.49 (-2.62, 1.65) | 0.82 (-2.59, 4.23) |
| Campus Point SMCA | 0.42 (-0.08, 0.92) | 0.51 (-0.16, 1.18) |
| Carrington Pt SMR | 1.51 (0.19, 2.83) | 1.33 (-2.20, 4.87) |
| Farnsworth Onshore SMCA | -0.90 (-1.66, -0.13) | 0.07 (-1.43, 1.58) |
| Gull Island SMR | 1.08 (0.49, 1.66) | 0.17 (-0.90, 1.24) |
| Harris Point SMR | 1.78 (0.76, 2.79) | 3.98 (1.92, 6.03) |
| Naples SMCA | -0.07 (-0.55, 0.41) | -1.02 (-1.56, -0.47) |
| Point Dume SMCA | 0.09 (-0.15, 0.33) | 0.33 (-0.05, 0.70) |
| Point Dume SMR | 0.67 (-0.31, 1.64) | 1.49 (0.31, 2.67) |
| Point Vicente SMCA | 1.16 (0.69, 1.63) | 1.16 (-0.26, 2.59) |
| Santa Barbara Island SMR | 1.97 (0.49, 3.44) | 6.19 (4.83, 7.56) |
| Skunk Pt SMR | -0.11 (-1.33, 1.10) | -1.00 (-3.78, 1.78) |
| South La Jolla SMR | 0.97 (-0.07, 2.00) | 1.16 (-1.20, 3.51) |
| South Point SMR | -0.53 (-1.03, -0.03) | -1.32 (-2.77, 0.12) |
| Swamis SMCA | -0.53 (-1.34, 0.28) | 0.81 (-0.91, 2.53) |

## Figures and captions

- [All 16 estimates](all_pair_comparison.png): blue is the flat baseline, orange the changing baseline. Points are mean departures during observed After years; bars are pointwise 95% HAC intervals. Zero denotes no departure from the projected baseline. Rows are not independent because references can overlap.
- [Three focal trajectories](key_trajectories.png): points are annual log MPA/reference kelp biomass ratios. Solid lines are fitted trajectories; dashed lines continue the Before pattern into the After years. Dotted vertical lines mark designation years. Darker ribbons are pointwise 95% intervals for fitted means; lighter ribbons surround projected baseline means, not prediction intervals. The existing annual labeling treats designation years as Before. Y scales vary by site, but are matched across the two models within each site.
- All observed trajectories: [1](trajectories_1.png), [2](trajectories_2.png), [3](trajectories_3.png), [4](trajectories_4.png). Same caption; PDF versions are available beside each PNG.
- [Baseline-window sensitivity](baseline_windows.png): same colors and average-departure quantity, full versus final 10 and 5 calendar baseline years. All After observations are retained. Both models use identical observations within each window. Pointwise 95% HAC intervals; panel-specific y scales.

## Eligibility and reproduction

[Eligibility CSV](eligibility.csv) contains the full metadata and qualifications. The compact table below shows paired annual coverage; Before includes the designation year. “Missing” counts unobserved years between each pair's first and last usable year, not unsampled years outside that span. Protection categories describe legal take allowances, not enforcement or historic fishing intensity. Exact legal dates were not newly reconstructed; the table preserves the audited designation years used by the pipeline.

| Comparison | Protection | Observed span | Paired years | Before | After | Missing |
|---|---|---|---:|---:|---:|---:|
| Abalone Cove SMCA | Limited-take | 1984–2021 | 36 | 27 | 9 | 2 |
| Cabrillo SMR | No-take | 1985–2021 | 34 | 25 | 9 | 3 |
| Campus Point SMCA | No-take | 1985–2021 | 37 | 28 | 9 | 0 |
| Carrington Pt SMR | No-take | 1985–2021 | 37 | 19 | 18 | 0 |
| Farnsworth Onshore SMCA | Limited-take | 1990–2021 | 24 | 15 | 9 | 8 |
| Gull Island SMR | No-take | 1985–2021 | 37 | 19 | 18 | 0 |
| Harris Point SMR | No-take | 1985–2021 | 36 | 18 | 18 | 1 |
| Naples SMCA | Limited-take | 1985–2021 | 37 | 28 | 9 | 0 |
| Point Dume SMCA | Limited-take | 1985–2021 | 37 | 28 | 9 | 0 |
| Point Dume SMR | No-take | 1985–2021 | 37 | 28 | 9 | 0 |
| Point Vicente SMCA | No-take | 1984–2021 | 35 | 26 | 9 | 3 |
| Santa Barbara Island SMR | No-take | 1984–2021 | 37 | 19 | 18 | 1 |
| Skunk Pt SMR | No-take | 1985–2021 | 36 | 18 | 18 | 1 |
| South La Jolla SMR | No-take | 1985–2021 | 35 | 26 | 9 | 2 |
| South Point SMR | No-take | 1985–2021 | 37 | 19 | 18 | 0 |
| Swamis SMCA | Limited-take | 1984–2021 | 37 | 28 | 9 | 1 |

Run from the analysis repository root:

```sh
python3 code/tools/prepare_satellite_baseline_comparison.py
Rscript code/tools/compare_satellite_baselines.R
python3 code/tools/check_satellite_baseline_comparison.py
```

The frozen cohort and protection metadata make this comparison reproducible without reading manuscript files. source_manifest.json records original source hashes. PLAN.md was written before fitting. annual_input.csv and model_input_years.csv record exact observations. estimates.csv includes all four quantities, three windows, two models and three uncertainty lags; coefficients.csv, predictions.csv and baseline_holdout.csv support reconstruction.

Independent NumPy regressions and calendar-gap HAC calculations reproduce all 1,152 estimates, standard errors and interval endpoints (maximum absolute discrepancy <1.5e-11). Separate Before/After regressions verify the model interpretation. All 487 protected input/output/document fingerprints are unchanged. The seven PNG figures were visually inspected, including all 16 histories; no embedded titles, captions or legends. No manuscript change, primary-analysis replacement, send, commit or push.
