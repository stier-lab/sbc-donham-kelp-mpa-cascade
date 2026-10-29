# Malakhoff and Miller: best-available data update

16 October 2026. Exploratory reconstruction, not an exact reproduction or a
2025 update. No manuscript result or figure has been replaced.

For the subsequent same-data investigation of differences from the main paper,
see `docs/MALAKHOFF_2021_METHOD_BRIDGE.md`. It distinguishes the actual primary
estimands, restores density-only evidence from sites without earlier sizes,
and separates endpoint, covariance, reference and calendar choices.

## Findings

Red urchin biomass increased more at reserve sites than controls between the
before and after periods. Extending the reconstructed series from 2019 through
2022 changed this contrast little. The positive red-urchin contrast persisted
when each supported reserve was omitted in turn and under either period
boundary. Purple biomass showed no clear reserve-control change; its
uncertainty was broad and its estimated direction depended on the boundary.

This does not establish a trophic cascade. The summed-species contrast was
sensitive to the boundary and reserve membership. Neither species showed a
clear change in the biomass of individuals <35 mm when serial correlation
was included. Ignoring serial correlation produced narrower intervals and
apparent declines in small-urchin biomass; those simpler fits should not
supply the conclusion.

## Design and measurements

A fixed subset of **three reserves and eight controls** has at least two
before and two after years with measured biomass for each species under
both boundaries in 1991-2019. This original-window rule fixes membership
before extending the series; it does not require complete annual coverage.
Missing sizes are never interpolated or borrowed.

Excluded sites are Rodes Reef (5), Fry's Harbor (7), and Southeast Sea Lion
Rookery (14). The first two lack size-supported biomass histories; the third
lacks a before-period biomass baseline. Retained reserves are Hare Rock (2),
Gull Island (6), and Scorpion Anchorage (9). This follows the original pooled
control design, not the main manuscript's paired-reference analysis.

Wet biomass is count density times the count-weighted mean of masses calculated
from individual sizes. Small-urchin biomass uses only the mass of individuals
<35 mm, rather than assigning all sizes the same mean mass. That differs from
the original code's count-fraction allocation. The combined response is the
sum of both species, not the paper's stacked-species random-slope model.

The modified BACI model fits log(biomass + 1) to reserve status, period and
their interaction, with a site random intercept and calendar-year AR1 errors
within sites. Missing 2020 observations retain their actual time gap. The
table uses 2005 as the first after year, matching the archived code; a separate
sensitivity uses 2004, following the paper's wording. This covariance structure
is not the original year/site grouping and does not account for common year
shocks across sites.

## Estimates

Contrasts are on the log(wet grams per square metre + 1) scale, with nominal
95% model-based intervals. Positive values mean a larger before-after change
at reserves, not necessarily an absolute increase at every reserve. The +1
offset prevents interpreting the exponentiated contrast as a simple ratio of
unshifted biomass. Intervals are not adjusted for multiple comparisons.

| Response | Through 2019 | Through 2022 |
|---|---|---|
| Red biomass | 1.43 (0.66, 2.20) | 1.32 (0.58, 2.07) |
| Purple biomass | 0.45 (-0.84, 1.73) | 0.41 (-0.87, 1.70) |
| Summed biomass | 0.86 (0.11, 1.61) | 0.76 (0.02, 1.50) |

Through 2022, small-urchin contrasts were -0.09 (-1.03, 0.85) for red and
-0.23 (-1.54, 1.08) for purple urchins. These are size-specific mass estimates,
not direct reproductions of the paper's small-size endpoint.

Using 2004 as the first after year gave red biomass 0.93 (0.18, 1.68), purple
biomass -0.84 (-2.16, 0.48), and summed biomass 0.24 (-0.50, 0.98).
Omitting each reserve in turn from the through-2022, first-after-2005 fits gave
red contrasts of 1.04-1.66, with all three intervals above zero. The summed
contrast's interval included zero after omitting Hare Rock or Gull Island.
These checks are not additional independent replication.

Only two retained sites contribute observations in 2023. The through-2023
fit is therefore a partial-coverage sensitivity, not the principal update.
Red and purple contrasts were 1.36 (0.61, 2.10) and 0.43 (-0.85, 1.70).

## Checks and limits

- All 75 common-support fits converged. The AR1 fits' pooled normalized
  residual lag-one correlations had absolute values <=0.101. This is not a
  formal test or assurance about every site's residuals.
- The red small-urchin leave-Gull-out fit had an almost-zero site-intercept SD
  (0.00029). Its variance decomposition needs caution.
- Before-period slope-contrast intervals included zero for all five responses,
  but were broad. This does not demonstrate parallel trends.
- Residual-versus-fitted and QQ plots are exported for every fit. Convergence
  and numerical checks do not establish publication readiness.
- The rendered through-2022, first-after-2005 red, purple and summed-biomass
  diagnostic pages were visually inspected. QQ plots show heavy tails and
  extreme negative residuals (approximately -4.6 to -6.1). Model-based
  intervals therefore need caution; these checks have not established robust
  uncertainty. Not all 75 diagnostic pages were visually reviewed. Both
  comparison and observed-history figures were inspected for clipping.
- Three reserves, shared environmental shocks, spatial differences, missing
  sizes and fishing effort limit causal interpretation. Evidence remains
  exploratory and observational.
- Kelp stipe-density and kelp-urchin analyses remain unverified. Multiple kelp
  taxon IDs and survey methods lack a verified measurement crosswalk. They
  were not summed or relabelled as stipe density. Per-plant stipe observations
  were not treated as an area-based census.
- The identified 1982-2025 NPS package has restricted files. This update uses
  the local MBON export, not that package.

## Reproduce and inspect

From the analysis repository root:

```sh
Rscript code/replications/malakhoff2021/run.R
Rscript code/replications/malakhoff2021/test.R
```

The runner includes the fixed-cohort extension. To rerun just that extension
after rebuilding the annual inputs, run
`Rscript code/replications/malakhoff2021/common_support.R`.

All outputs below are in `outputs/malakhoff2021/`:

- `common_support_membership.csv` and `common_support_counts.csv`: inclusion
  decisions and species/boundary-specific before-after support.
- `common_support_inputs.csv`: five endpoint series, preserving missingness.
- `common_support_results.csv`: all 75 estimates, intervals, diagnostics and
  omitted-reserve identities, with manuscript eligibility false.
- `common_support_pretrends.csv`: before-period slope contrasts.
- `common_support_trajectories.png`: observed site histories and descriptive
  group means through 2022; breaks mark unsampled years. These are not model
  predictions and group means can have different annual coverage.
- `common_support_comparison.png`: internal audit plot of window-specific
  contrasts and intervals, not a proposed manuscript figure.
- `common_support_diagnostics.pdf`, `common_support_models.rds` and
  `common_support_sessionInfo.txt`: diagnostics, fitted models and runtime.

Source audit and exact-replication requirements:
`docs/MALAKHOFF_2021_REPLICATION.md`.
