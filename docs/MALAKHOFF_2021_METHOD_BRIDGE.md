# Why the cascade evidence differs

16 October 2026. Diagnostic analyses only. The manuscript and its primary
outputs have not been changed. This is not an exact replication of Malakhoff
and Miller (2021), nor a test using the restricted 1982-2025 release.

## Main finding

The analyses do not show opposing demonstrations of a trophic cascade.
They measure different responses and, importantly, different kinds of change.
The network synthesis contains negative purple-urchin density contrasts, but
most are not designation-associated before-after effects. The reconstructed
KFM comparison does not resolve a consistent purple-urchin BACI response.
Neither result by itself establishes predator-mediated kelp recovery.

The new checks demonstrate several reasons the results need not agree.
They do not assign a percentage of the discrepancy to each cause: changes in
sites, responses, transformations and weighting interact.

## 1. The current synthesis mixes three questions

The 13 purple-urchin density estimates in `data/sumstats_final.csv` comprise:

| Current primary estimate | Count | What it measures |
|---|---:|---|
| CI, Linear | 9 | Fitted relative change over 11 model-time intervals |
| CI, Mean | 3 | Mean protected/reference log contrast over retained input years |
| pBACIPS, Step | 1 | Before-after change in the protected/reference contrast |

These are MPA-programme estimates, not 13 independent MPAs. The nine linear
estimates span seven MPAs. All current CI estimates for both urchin species
were reconstructed numerically from their actual annual histories, matching
the stored primary estimates to <1e-8. In particular, a CI linear estimate is
11 times the fitted slope, not the post-protection protected/reference level.
Most diver clocks count survey intervals rather than elapsed calendar years.
Some CI histories retain before-designation observations as well; calling
every CI estimate "after-only" is inaccurate.

For example, Abalone Cove's PISCO purple-density primary estimate is -2.81,
but its mean after-period log contrast is +0.43. South Point's PISCO values
are -1.41 and +0.38, respectively. A declining relative trajectory can coexist
with higher abundance inside a reserve. These are different quantities, not
contradictory measurements of the same quantity.

The selection rule adds another limit: CI linear effects are selected when
their ordinary-regression slope is significant; otherwise the mean effect is
selected. The reported intervals do not account for that selection. None of
the new checks treats the selected estimates as an unbiased set of common
designation effects.

## 2. Density and biomass can tell different biological stories

The balanced comparison holds three reserve and eight control sites fixed
and retains the same 29 calendar years for both species and responses.
The original-window extension compares 1991-2019 with 1991-2022. Years 1993
and 2003 are excluded from this balanced layer because of missing biomass;
2020 has no counts in the local export. True zeros remain in its model fits.

A separate positive-observation check uses the identity
`log(biomass) = log(density) + log(mean individual mass)`.
Within each species, it retains only years with positive, size-supported
observations at all eleven sites. The same observations and weights enter all
three terms. Equal-site before-after reserve-control contrasts through 2022,
with 2005 first after, are:

| Species | Density contribution | Mean-mass contribution | Biomass contrast |
|---|---:|---:|---:|
| Red | +0.283 | +0.953 | +1.236 |
| Purple | -0.264 | +0.667 | +0.404 |

The terms add exactly before rounding. This supports a biological explanation
for the apparent mismatch: a relative increase in individual mass can offset
a decline in numbers. For red urchins, much of the positive biomass point
contrast comes from mass per individual, not more individuals.

This is an accounting identity, not evidence that particular individuals grew
or that protection caused a size change. Positive-year filtering changes the
sample: red retains 10 before and 16 after years; purple retains 9 and 14.
Welch intervals among site-level changes are broad and include zero for all
six contributions; even the red biomass interval is -0.043 to 2.514. Three
reserves do not justify a precise size-mechanism claim. The decomposition is
not a replacement for the zero-inclusive BACI fits.

## 3. Density-only evidence needs its own site cohort

Santa Barbara Island supplies the manuscript's sole purple-density pBACIPS
estimate. Its reserve, site 14, lacks the earlier sizes required for the
biomass reconstruction, but has earlier density observations. Excluding it
from a biomass-supported cohort must not erase its density evidence.

The density-only sensitivity therefore restores all four reserves and ten
controls in the original study, holding the 29-year calendar above fixed.
With calendar AR1 errors, the through-2022 purple-density BACI point estimate
ranges from +0.066 to -0.536 as the zero correction varies from 0.01 to 1
individual per square metre. All three nominal 95% intervals include zero.
Changing the correction changes the transformed response, not just its units.

Reserve responses differ substantially. With correction 0.1, descriptive
before-after contrasts against all ten controls are +1.01 at Hare Rock,
-0.99 at Gull Island, -1.06 at Scorpion and -1.66 at Southeast Sea Lion Rookery.
The last site's contrast against the two mapped Santa Barbara Island controls
is -1.85. These individual descriptive changes are not additional independent
tests, but explain why one strong paired result need not establish a common
four-reserve response.

## 4. Pairing is not itself the explanation

On the same positive observations, assigning the same two reserve and two
control sites into pairs leaves the equally weighted mean log contrast
unchanged, to numerical precision. Averaging log pair ratios and subtracting
the two group mean logs give the same answer. Pairing can affect uncertainty;
this identity concerns the point estimate only.

Changing which controls enter, or taking logs after arithmetic averaging,
does change the answer. In the two-reserve experiment, the red-biomass
before-after contrast is 0.964 against all eight supported controls, 0.712
against the two mapped controls using mean logs, and 0.414 using the log of
arithmetic group means. The same two reserves and species-specific positive
calendar are held fixed across these descriptive calculations.

The two-pair experiment uses Harris Point and Scorpion legacy mappings; it
does not validate their ecological matching. Gull is excluded from this
experiment because its KFM reference assignment is disputed. None of these
calculations resolves that contributor question or changes the manuscript's
Gull membership.

## 5. Autocorrelation matters, but does not explain away the pooled result

The 23 currently selected CI linear fits for both urchin species were refitted
with calendar-year AR1 errors, retaining their annual observations and the
same 11-model-interval contrast. Each GLS fit tries multiple starting
correlations and optimizers using the existing repository helper.

For purple density, nine of nine original IID intervals exclude zero in the
negative direction; five of nine AR1 intervals do so. These are conditional
intervals for already-selected models, not a corrected selection procedure.
Two of the nine AR1 correlations approach the boundary, including a Scorpion
fit with an invalid approximate parameter covariance. They cannot be treated
as validated replacement fits.

Replacing the selected urchin linear estimates and SEs in an otherwise
unchanged joint synthesis retains a negative purple-density coefficient:
-1.847 (-3.020, -0.675), compared with -1.990 (-3.128, -0.852) in the current
primary model. A separate check retaining all original point estimates and
only increasing the corresponding SEs gives -1.855 (-2.968, -0.742).
The original joint coefficients were reproduced before making substitutions.
The broader density pattern therefore does not disappear under these checks.

These are partial covariance sensitivities, not publication-ready corrections:
other model classes and taxa are held fixed, boundary fits remain flagged,
and the mixture of estimands and model-selection uncertainty is unresolved.

## 6. Timing and geography remain limited explanations

Era-specific KFM models distinguish 2005-2013, 2014-2016, 2017-2019 and
2021-2022 from the before period. The later purple-density point contrasts
are more negative than the first after-period contrast, but every era's
interval includes zero. This does not establish that averaging across eras
hid a well-supported purple-urchin response.

A PISCO-only comparison requiring every year in 2013-2022 leaves four
purple-density MPAs: one island and three mainland, with shared reference
sites reducing replication further. Biomass has no island MPA meeting that
strict calendar rule. Thus the available records cannot cleanly separate
geographic coverage from programme, calendar coverage and reef selection.
The earlier, less restrictive 2013-2023 geography check is retained as context,
not substituted when the strict comparison becomes sparse.

## Consequences for the paper

The defensible contrast is between a broader network pattern of relative
urchin trajectories/levels and a narrower designation-focused BACI test.
We should not say our method reveals a cascade that Malakhoff and Miller
missed, or explain the difference solely by statistical power or geography.

The main text should distinguish the pooled purple-density pattern from a
network-wide causal decline following protection. The Methods and supplement
should continue identifying means, fitted trajectory changes and before-after
steps separately. A full reassessment would fit a common estimand across the
same source histories, with model selection and shared-control dependence
handled explicitly, rather than relabelling the existing pooled coefficient.

The trophic mechanism remains a separate question: the kelp stipe-density
endpoint is still unverified in this reconstruction, and none of these
urchin-only checks demonstrates predator-driven kelp recovery. Both studies'
mechanistic limits need to remain visible.

## Reproduction and audit files

`run.R` now calls both `common_support.R` and `method_bridge.R` in the standalone
replication workflow. Run from the analysis repository root:

```sh
Rscript code/replications/malakhoff2021/run.R
Rscript code/replications/malakhoff2021/test.R
```

To rerun just the bridge after preparing annual inputs, use
`Rscript code/replications/malakhoff2021/method_bridge.R`.
Dependencies include data.table, nlme, ggplot2 and metafor. All new results
are under `outputs/malakhoff2021/bridge/`:

- `primary_estimand_inventory.csv`, `primary_estimate_reconstruction.csv` and
  `primary_annual_input.csv`: the actual quantities entering the current paper.
- `calendar_support.csv`, `balanced_input.csv`, `balanced_model_results.csv`:
  80 controlled response/allometry, covariance, boundary and window fits.
- `positive_calendar_support.csv`, `size_density_site_changes.csv`,
  `size_density_decomposition.csv`, `decomposition_identity.csv`: exact
  accounting and its positive-observation support.
- `density_cohort_results.csv`, `density_reserve_contrasts.csv`: 48 additional
  all-four-reserve density fits and descriptive reserve-specific contrasts.
- `reference_experiment_pairs.csv`, `reference_annual_contrasts.csv`,
  `reference_comparison.csv`: membership and reference/aggregation algebra.
- `era_model_results.csv`: four models and their 16 temporal interactions.
- `primary_linear_covariance_check.csv`, `primary_linear_fit_attempts.csv`:
  46 fits, all optimizer attempts, boundary and covariance flags.
- `joint_sensitivity_inputs.csv`, `joint_sensitivity_results.csv`: three
  otherwise matched synthesis scenarios. Model RDS files are retained.
- `geography_complete_calendar_coverage.csv`, `geography_complete_calendar_pairs.csv`,
  `geography_complete_calendar_summary.csv`: strict-calendar exclusions,
  MPA-level values and shared-reference-aware descriptive summaries.
- `red_size_density_histories.png` and `purple_size_density_histories.png`:
  observed site histories, with thick arithmetic group means. Mean-mass
  coverage can differ by one site-year; lines are not fitted predictions.
  Missing years break the trajectories. Figures have no titles or subtitles.
- `diagnostics.pdf`, `diagnostic_index.csv`, model RDS files, `sessionInfo.txt`
  and `input_manifest.csv`: diagnostics, runtime and unchanged-input checks.

The observed-history figures and three through-2022, first-after-2005
diagnostic pages (red biomass, purple density with offset 0.1, purple biomass)
were rendered and visually inspected. Residual QQ plots have heavy tails and
extreme observations; the nominal intervals should not be treated as robust.
The remaining diagnostic pages were exported, not individually reviewed.
`model_audit.csv` separately records correlation boundaries, approximate
covariance validity, near-zero site-variance fractions and residual extremes
for every controlled-fit and era model.
Of these 132 models, five have non-positive-definite approximate parameter
covariances and two have near-zero site-variance fractions, with one overlap.
These include the four-reserve purple-density fit with offset 0.1 through
2022. They remain flagged diagnostic results, not validated replacements.

Tests check matching site/year support, preserved zeros, both log identities,
primary-estimate reconstruction, unchanged synthesis membership, variance-only
substitution rules and the strict geographic calendar. Passing these checks
does not clear the analyses for manuscript use or validate model assumptions.
