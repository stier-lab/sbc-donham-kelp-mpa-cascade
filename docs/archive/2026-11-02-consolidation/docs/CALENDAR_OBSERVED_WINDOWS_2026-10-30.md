# Population changes within surveyed calendar years

The primary trend contrasts use elapsed calendar years and two surveyed post-designation endpoints. Each contrast starts at the first post-designation survey and ends at the latest survey no more than eleven years later. All retained observations enter the fit. Baseline observations have calendar time zero in the progressive model; baseline-trend screening uses actual survey years. Mean-level and after-minus-before step contrasts retain their definitions.

This implements the October 30 decision to avoid extrapolating primary population changes beyond the years surveyed. It retains all 145 histories, all 2827 annual observations, the existing program/design assignments and selected model forms. It does not change reference selection, PISCO eligibility policy, the designation-year boundary or the heatwave windows. The 73 trend contrasts span 4–11 years (40 span 11 years); the other 72 are means or steps. Both endpoints are surveyed for every trend. The two four-year contrasts are Matlahuayl PISCO red- and purple-urchin density, 2019–2023.

## Findings

Forty-two population estimates change. All nine pooled directions persist. Higher kelp and sheephead biomass and lower purple-urchin density and biomass remain supported after correction across nine tests. Red-urchin biomass remains positive in direction but its adjusted P changes from 0.040 to 0.080. Its primary log contrast is 0.486 (nominal 95% CI 0.012–0.960). Automatic calendar-HAC uncertainty gives adjusted P 0.047; red-urchin evidence therefore depends on uncertainty assumptions.

The four-MPA lobster–purple-urchin biomass slope is −0.555 (nominal 95% CI −1.032 to−0.078; P 0.038), with exact permutation P 0.0833. This remains limited evidence for a trophic relationship, not demonstration of a complete cascade. The targeted size checks preserve all nine pooled directions and four adjusted classifications across 2500 perturbations and 38 specified stress cases. The pathway's nominal interval includes zero in 38.9% of individual-weight draws and 27.3% of recorded-group draws; these are conditional diagnostics, not total confidence intervals.

## Reproducibility

- `08_effect_sizes.R` fits calendar-time trends and selects actual endpoints. Nonlinear contrast covariance uses Richardson differentiation to avoid finite-difference cancellation. An unavailable candidate p-value cannot shorten the candidate-diagnostic table.
- `primary_uncertainty_utils.R` reconstructs current estimates/SEs and their gradients for downstream HAC, station, young-sheephead, missing-lobster and size-sample checks.
- Module 52 keeps time-window alternatives as explicit sensitivities and uses current observed endpoints in alternative PISCO routes.
- Module 58 independently validates every fitted history's calendar support and exports `table_s_primary_calendar_windows.csv`. It runs in both primary and review runners.
- Modules 31/32 refresh pathway and fitted-history displays. The size-filtered pathway helper also uses observed calendar endpoints; the annual temporal and heatwave analyses retain their existing calendar definitions.
- `code/tools/audit_calendar_transition.py` compares current estimates with the frozen preceding fitted results in `data/external/calendar_policy_20261030/`, checking unchanged identities, sample sizes, design labels and model forms. Those files are historical analysis comparators, not input observations.

Current comparison exports: `tables/table_s_calendar_transition_histories.csv` and `tables/table_s_calendar_transition_pooled.csv`. The window ledger records each trend's actual start, end and duration. Reference registry metadata was regenerated against the same source assignments after the primary-result fingerprint changed; no sites were added or removed.

## Interpretation limits

Observed durations differ and may start well after designation. A fitted change between two post-designation survey years is not a before/after designation effect. Combining means, steps and progressive changes still does not yield a common fold increase or a common-duration causal effect. Fixed historical model assignments do not include model-selection uncertainty; this timing update does not settle the remaining eligibility/reference questions. Seasonal refits remain deferred.

The permutation check now uses the same equal-source MPA means as Table 3 and permutes predator labels while keeping each response estimate with its own sampling uncertainty. All six observed slopes must reproduce Table 3 before shuffling. This repairs a separate diagnostic mismatch uncovered during concordance checks; the focal four-MPA relationship remains unresolved.
