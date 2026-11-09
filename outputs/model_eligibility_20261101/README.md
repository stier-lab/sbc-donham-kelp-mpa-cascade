# Baseline eligibility audit for review

November 1, 2026. No manuscript, supplement or canonical analysis outputs were updated. No biological tolerance was selected.

Site inclusion and model eligibility are separate decisions. The researcher/ecological exclusions remain fixed, leaving 175 candidate histories and 162 estimable post-protection comparisons. Of these, 54 previously supported a before–after contrast. The baseline audit flags seven linear trends, nine curvature/level shifts, eight sparse baselines and one constant baseline. The other 29 have no detected departure, which does not establish stability.

Thiault et al. require assessing baseline stationarity and other BACIPS assumptions; they do not prescribe the numerical thresholds used here. Sources: [2017 methods](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12655) and [2019 application](https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecs2.2576).

## Several tolerances, none adopted

The audit asks whether the uncertainty in the baseline slope is contained within a specified change in the protected/reference ratio over nine years. Bounds are reciprocal: factor 1.5 means 0.667–1.5, not ±50%. The 90% slope envelope uses the largest standard error among IID, HC3 and calendar-gap HAC lags 1, 3 and 5. This is an exploratory conservative envelope, not a newly validated inference procedure. Additional safeguards require six Before observations spanning five years, no detected linear/curvature/shift flag, and all leave-one-year-out slope estimates inside the bounds. These are analyst choices, not Thiault's rules or changes to site inclusion.

| Factor | Ratio bounds over nine years | Histories compatible with mean-trend screen |
|---|---|---|
| 1.10 | 0.909–1.10 | 0 |
| 1.25 | 0.80–1.25 | 0 |
| 1.50 | 0.667–1.50 | 1: Campus Point satellite kelp biomass |
| 2.00 | 0.50–2.00 | 5: Campus Point, Abalone Cove, Point Dume SMCA and Point Vicente satellite kelp biomass; Scorpion KFM purple urchin density |

Relaxing the coverage safeguard to four Before observations over three years leaves the factor-1.5 set unchanged. The three-year PISCO histories are unresolved for this stability screen; they are not excluded sites. Constant Scorpion KFM kelp observations cannot estimate natural baseline variation and remain unresolved rather than being declared stable automatically.

Even the five factor-2 histories are not certified pBACIPS fits: Point Vicente has a variance-change flag; Scorpion purple urchins have residual serial correlation and borderline curvature. Sampling comparability and additivity also require the existing design evidence. Multiple exploratory diagnostic tests can produce chance flags. All subset estimates are conditional on a screen that uses the same baseline contributing to the effects; their intervals do not include selection uncertainty.

## What changes in the pooled results

The following estimates reuse the same period-average log-ratio Change effects and years and only change membership. They are not refits of the original nonlinear pBACIPS/AICc curves. A final pBACIPS claim requires that separate step after the eligibility policy is settled.

| Baseline screen | Kelp histories | Pooled change, 95% interval |
|---|---:|---|
| All before–after candidates | 22 | 0.341 [−0.147, 0.829] |
| No detected linear trend (permissive diagnostic only) | 20 | 0.294 [−0.237, 0.825] |
| Longer baseline, no linear/shape flag | 19 | 0.276 [−0.263, 0.815] |
| Factor 1.5 | 1 | 0.421; no pooled inference |
| Factor 2 | 4 | 0.756 [−0.283, 1.796] |

The kelp intervals all include no change. Stricter screens greatly reduce information rather than settling whether protection increased kelp. The factor-2 purple urchin result is one history (0.171), with no pooled inference. Most other responses become too sparse to support pooled inference under the longer-baseline screen.

All post-protection contrasts and pooled inference are unchanged: sheephead biomass and density remain higher inside protected areas, but these differences alone do not attribute that advantage to protection. This audit does not establish a general protection-driven trophic cascade, nor does low power establish that cascades are absent.

## Files and figure captions

- `COMPARISON_DECISIONS.md`: readable decisions for all 54 before–after histories.
- `model_eligibility.csv`: all 175 histories, years, diagnostic statistics and tolerance membership.
- `pooled.csv`: all responses, programs and scenarios, including small-sample inference limits.
- `kelp_model_eligibility.png` / `.pdf`: pooled kelp change under each baseline screen. Points show estimates; bars show 95% cluster-robust intervals where supported. The open point is the single Campus Point history, without pooled inference. Zero is no change in the protected/reference log ratio. No sample-size or color key is embedded in the artwork.
- `baseline_atlas.pdf`: all 54 Before histories. Red indicates a detected linear trend; orange a curvature/shift flag; gray sparse or constant baselines; blue no detected departure, not proven stability. Dashed lines show fitted linear baselines. Species, program and response identifiers distinguish each comparison.

## Reproduction and checks

Run from the analysis repository:

```sh
python3 code/tools/audit_baseline_model_eligibility.py
Rscript code/tools/pool_model_eligibility.R
Rscript code/tools/check_baseline_eligibility.R
python3 code/tools/check_common_estimand.py outputs/model_eligibility_20261101
python3 code/tools/plot_model_eligibility.py
python3 code/tools/report_model_eligibility.py
```

Independent R reconstruction verifies baseline linear slopes and the standard-error envelope (maximum absolute difference 4.23e−10). The contrast/pooling checker verifies 1,432 individual effects and 209 pooled rows; 487 protected files remain unchanged. Membership checks confirm nested tolerance sets, unchanged inherited effects and unchanged post-protection pooling within 1e−8 numerical tolerance in all eight scenarios. All six atlas pages and the pooled figure were visually inspected. Quadratic and shift diagnostics have not received a separate independent numerical reconstruction.
