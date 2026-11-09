# Primary and Temporal Uncertainty Checks

Implemented 12 October 2026. These checks supersede the pooled Durbin-Watson
variance-inflation reconstruction as evidence about primary uncertainty. The
older reconstruction remains exported for historical comparison, not validation.

## Exact Primary Cohort

Module 36 reconstructs all 145 retained histories, including three nonlinear
fits, using the selected model, native model clock, zero correction and primary
observations. Every contrast and independent-error SE matches the primary
export; reusing original variances reproduces Table 2's joint crossed-MPA/source
REML model, coefficient t tests and nine-test Benjamini-Hochberg adjustment.
No input, selection, primary point estimate or exclusion changes.

For each selected model, `sandwich` 3.1-3 estimates calendar-lag HAC covariance
with Bartlett weights and lags 1, 3 and max(1, floor(4(N/100)^(2/9))). Missing
calendar years have zero estimating scores, not imputed responses. We use the
N/(N-p) finite-sample factor and propagate the full covariance through the
primary contrast gradient. Each synthesis variance is at least the original
variance. Unfloored SEs are also exported. The same joint synthesis is refitted
for each bandwidth, including its variance components and weights.

All nine pooled directions and BH detection classifications persist under all
three bandwidths. With automatic lag, the median individual synthesis-SE ratio
is 1.0366 and the maximum is 1.5309; 64 of 145 variances are floored. Pooled
intervals do not all widen because REML components and weights are reestimated:
the automatic pooled interval-width ratios range from 0.9143 to 1.0311
(median 0.9541). This is a conditional sensitivity check, not guaranteed coverage
after model selection and not a solution for shared references, zero offsets,
short histories, endpoint extrapolation or observational attribution.

## Figure 4 Cohort and Structure

The former biomass fit used 978 rows from 74 histories, including 147 rows in
11 histories outside the retained primary cohort. Both fitting and plotting now
use the same retained diver histories: 831 biomass rows, 63 histories and
15 MPAs. The separate density fit uses 919 rows. Satellite data do not enter
Figure 4. Models use all post-designation years (1-20); panels display 723
biomass observations through year 15 from the same 63 histories. Thin lines
identify MPA-source-taxon histories rather than joining
different source programmes within an MPA.

The biomass and density fits both retain correlated MPA intercepts/slopes and
crossed source intercepts, converge and are nonsingular. A recorded
rich-to-simple hierarchy chooses the first converged nonsingular structure,
not the first significant result or minimum AIC. All five candidate structures
and their failures/boundaries are exported. Slopes and every pointwise mean
prediction have their own linear-contrast Satterthwaite df. The plotted bands
remain model-based intervals under independent annual errors. Residual plots
now diagnose the actual Figure 4 biomass model, not the pooled response fit.

## Same-Row Calendar AR1 Comparison

Module 37 uses Gaussian REML in `glmmTMB` 1.1.14 on the same 831 rows. Explicit
year levels retain calendar gaps; AR1 processes are separate for each retained
MPA-source-taxon history, with an additional independent annual-error term.
The correlated-MPA-slope AR1 candidate fails
convergence and Hessian checks and is rejected. Independent MPA slopes with
crossed source intercepts converge with a positive-definite Hessian. The
estimated AR1 process correlation is 0.6771, not the total residual correlation.
Candidate exports record process/error SDs and the resulting within-history
lag-one correlation conditional on MPA/source effects. Source SD is near zero (0.00053), with
only three contributing programmes; this limitation is retained in the SI.

All five slope signs persist, but purple-urchin biomass is no longer clearly
negative: slope -0.09090, conditional Wald 95% interval -0.18358 to 0.00178.
Kelp remains positive: slope 0.16525, interval 0.09221 to 0.23829. These AR1
intervals are model-based Wald intervals, not small-MPA calibrated intervals.
Figure 4 describes annual trajectories and does not independently establish
a causal predator-grazer sequence. The primary pooled contrast synthesis and
the temporal slope comparison answer different questions.

## Reproduce and Inspect

```sh
Rscript --vanilla code/R/run_all.R
Rscript --vanilla code/tests/test_primary_uncertainty.R
Rscript --vanilla code/tests/test_temporal_uncertainty.R
python3 -B code/audit_primary_uncertainty.py --before /path/to/preserved/before
```

The full pipeline ran successfully in 8.85 minutes in
`logs/pipeline_20261012_144853.log`. Targeted replay regenerated Figure 4 and its
biomass residual diagnostic after reporting/legend cleanup. The independent
audit records unchanged input/reference bytes, all primary effects, Tables 2/3,
modern/historical kelp results, cohort membership, reconstruction, BH arithmetic,
model-selection validity and temporal interval arithmetic. R unit checks also
verify Newey-West equivalence on regular years, calendar gaps, positive
semidefiniteness, duplicate rejection and analytic contrast means/SEs.

Primary outputs are `table_s_primary_hac_histories.csv`,
`table_s_primary_hac_pooled.csv`, `table_s_primary_hac_inputs.csv`,
`table_s_temporal_candidates_bio.csv`, `table_s_temporal_candidates_den.csv`,
`table_s_temporal_predictions.csv`, and `table_s_temporal_ar1_*.csv` in `tables/`.
The same-row observations and former-cohort comparison are exported separately.
Original dirty files were preserved under
`local_archive/uncertainty_checks_2026-10-12/before/`; they were not reverted.
Current manuscript/SI construction remains in the private manuscript repository.
# Current-input follow-up

The 18 October singleton-biomass correction supersedes the numerical conclusions
below for current inputs. All nine pooled directions persist, but purple-urchin
biomass crosses the adjusted threshold under automatic and three-year HAC
bandwidths, not under the primary or one-year specification. See
`BIOMASS_SINGLETON_REFIT_2026-10-18.md` and its machine-readable audit. The
12 October record below remains a historical check of its then-current inputs.

