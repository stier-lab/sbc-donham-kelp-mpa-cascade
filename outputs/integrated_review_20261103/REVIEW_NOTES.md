# Reading the integrated review

This is a review of the corrected data and frozen rules. The manuscript and supplement have not changed.

## Coverage and the seven-day decision

The candidate inventory retains all 31 researcher/ecological exclusions. Of 175 retained histories, 62 support before–after models, another 100 support post-protection trends only, and 13 lack enough supported post-protection years. Before–after histories also have separately estimated post-protection trends. Missing sizes do not erase valid density observations.

Allowing same-site/year KFM urchin and kelp sizes within seven days recovers 16 comparison-years relative to exact-date matching. These are additional observations within histories, not 16 new sites. See [the recovered years](seven_day_recovered_years.csv) and [all fitted sensitivities](pooled_coverage_sensitivities.csv).

## What the patterns say

- After protection, the kelp protected/reference ratio increases by +16.9% (95% interval +4.9% to +30.3%) per year.
- Purple-urchin density becomes lower relative to reference sites: -9.7% (95% interval -16.9% to -1.7%) per year.
- Red-urchin density also has a negative estimated trend, but its interval endpoint is close enough to zero to change sides between random seeds. Do not turn that small numerical distinction into a different biological conclusion.
- The pooled sheephead and lobster post-protection intervals include both increases and decreases. Individual programs and histories can differ; use the source-specific results rather than interpreting the pooled interval as no response anywhere.

These trends describe change relative to reference sites. They do not establish that protection caused those changes or that a predator-driven cascade occurred.

## Before–after interpretation

The flat-baseline model estimates a kelp ratio change of +279.7% (95% interval +55.5% to +826.7%) over observed post-protection years. Allowing the before-protection trend to continue gives +9.5% (95% interval -93.7% to +1792.8%). That difference matters: the inference depends on the assumed trajectory without protection.

The corresponding flat-baseline lobster density estimate is +229.3% (95% interval -39.6% to +1695.4%). Its positive estimate is imprecise. Before–after estimates and post-protection trends answer different questions and remain separate.

## Numerical and release checks

Across 134 pooled summaries, changing the measurement-draw seed moves point estimates by at most 0.11 standard errors. The earlier large KFM lobster instability was traced to finite-draw correlations between independent survey years. Known-zero measurement covariances are now enforced; biological temporal dependence is still estimated. Variance fitting uses multiple starts and a checked temporal-correlation profile.

Intervals remain conditional on the represented sampling units, conversion equations and fitted covariance. Group-mean sizes, allometry, detection, sparse/zero samples and missing habitat information limit what the files can establish. Boundary flags remain in the diagnostics. The [annual input decomposition](input_change_decomposition.csv) separates coverage and measurement changes using the same simple diagnostic trend; those slopes are not substitutes for the joint model results.

Before updating the paper, review these results and the remaining offset, satellite-timing and habitat limitations, then regenerate downstream pathway/resilience analyses, figures and concordance checks. KFM redistribution permission remains a separate requirement. The current tables do not certify a final paper release.
