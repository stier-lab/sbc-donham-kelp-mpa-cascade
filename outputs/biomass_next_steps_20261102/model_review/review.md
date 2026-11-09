# Revised biomass model re-review — 2026-11-02

## Verdict

**The revised split resolves the prior latent-posterior-covariance recycling error. No newly demonstrated numerical blocker to completing the current, explicitly exploratory full run was found.** This is a bounded code/synthetic review, not verification of the full run or confidence-interval coverage. Keep the finite-year latent functional and the direct-only working process model separate in every table/caption/interpretation.

Reviewed the current 226-line `code/tools/biomass_observation_model.py`, frozen here as `reviewed_model.py`, SHA256 `d57d0fdd3640c2486105bf1b85d2a130c3a89165c45d87ba34050035a0d20510`. All line references are to this revised snapshot, not the first review. Repository and continuity records untouched; parent session `9901638869734200a6633d0f88ec31ef` retained. No sends/push/subagents. The parent's full 800/40/500 run and reported lobster eligibility counts were not independently rerun here.

## Revised split: statistically defensible within the stated scope

1. **Latent finite-year slope (174–179):** each column of the joint annual latent draws is regressed on centered actual calendar year. This preserves calendar gaps and within/between-year draw dependence, uses one fixed support set per series, and adds no second process covariance. Independent centered-cross-product slopes matched saved mean and lower quantile within 6.25e-17. Call its interval a *conditional model-draw interval for the slope of the observed-year latent trajectory*, not a process confidence interval or exact posterior credible interval. First-stage temporal shrinkage and approximate bootstrap variance-parameter mixing remain part of its assumptions.

2. **Separate direct-only process diagnostic (180–220):** only post years whose entire visit roster has direct same-visit sizes or observed zero counts enter. Size sampling replicates are centered on original arithmetic sample means via `mean * exp(-v/2 + sqrt(v)*Z)`, independent of latent size draws. Their annual log-ratio covariance, combined with empirical count weights, supplies a plug-in measurement covariance for the unsmoothed ratio point estimates. No latent posterior covariance enters `fit_trend`. The Gaussian likelihood and fixed covariance remain working approximations; the explicit interval label is appropriate. “Independent” applies to size/count simulation mechanisms, not all annual estimates: preserve shared count-unit/sample dependence.

3. **Eligibility is genuinely separate:** synthetic fixtures with one unsupported year and one borrowed year yielded five latent supported years but only four direct process years. Direct point estimates reproduced the independently known zero log-ratio exactly. No borrowed year entered the process fit. Current actual lobster 8/5-year eligibility was reported by the parent, not established by these fixtures.

4. **Log-mean adjustment:** `log(mean)+v/2` is now consistently used in training/main size models. Held-out observed-mean predictions transform back with `-vte/2` and add sampling variance. Direct sampling replicates use `-v/2` to preserve the original arithmetic mean. A 100,000-draw scalar check with mean 100 and v=.7 gave mean 99.8023, consistent with the intended expectation. This fixes the lognormal moment-matching offset; it does not validate the lognormal approximation for tiny samples.

## Failure/status checks: repaired, with reporting cleanup needed

- Base size and held-out fit failures now stop; no silent base-fit acceptance (105,116).
- Every parameter chunk now performs an actual bootstrap refit (119–126), rather than including the original fit as the first chunk.
- A forced first-attempt failure was retried and saved as `attempt_converged=False`, `retried=True`, `converged=True`, `used_fallback=False`. No fallback substitution occurred.
- Forcing both attempts to fail raised the intended `RuntimeError`; no final policy file was produced in a fresh output directory.
- Forcing a process fit failure stopped via assertion; no final policy file was produced. Process/bootstrap/LOO success checks are present (208,211,216).
- Rank-deficient single-zone and constant-year designs now raise explicit errors; 1–3-year process fits are rejected. Empty prediction targets return empty arrays.
- Missing canonical visit rows and an unsupported positive visit with finite mass both failed the guards.

**Required reporting cleanup:** line 133 still calls the first-attempt-failure count `bootstrap_fallbacks`, although these are retries and no fallback is used. The test produced `bootstrap_fallbacks=1` alongside `used_fallback=False` everywhere. Rename to `bootstrap_initial_failures` or `bootstrap_retries`; do not describe the current field as fallback use. Save final retry optimizer messages as well as initial messages if practical.

**Run-completion integrity:** failures leave partial outputs. In the forced process failure, `annual_observations.csv` existed but final policy did not. Reusing an output directory can leave old policy/trend files that appear complete after a failed new run. Require successful process exit plus an externally recorded run status/source hash and exact expected row counts; prefer a fresh directory or atomic completion manifest. Do not interpret the existence of a CSV/policy alone as success. Numerical exceptions still abort rather than retry, which is fail-closed; it is a robustness/diagnostic limitation rather than silent invalid inference. Assertions must remain enabled (do not execute with `python -O`).

## Remaining limitations and small repairs

### Held-out coverage is conditional on held-out dispersion

Training pooled CV and model fits exclude the held-out site (101–105), as intended. However, `vte` uses `te.relative_ss`, `te.df`, and `te.n` as well as the training pool (103). Thus interval widths and centers use dispersion estimated from the held-out size observations. This is not leakage into the latent fit or training pool, but it is not a completely outcome-blind predictive interval. Describe it as an **observed-sample prediction check conditional on estimated held-out sampling variance**. If claiming prediction before observing held-out sizes, use training-pool dispersion alone for the evaluation sampling model (with observed/design n) and report this as a separate check. Coverage is plug-in, conditional, and clustered within held-out sites; no protected-site validation follows from the reference results. All-singleton/zero-dispersion training pools remain undefined and stop; report the affected fold explicitly rather than interpreting it as biological failure.

### Direct working likelihood still needs explicit approximation limits

Unsmoothed empirical log ratios are not exactly unbiased Gaussian measurements: the nonlinear log/offset transformation can introduce sampling bias, particularly with sparse counts/high CV. The code uses the original ratio as its point, not the replicate mean, and uses replicate covariance only. This is acceptable as the advertised working diagnostic but should be checked by exporting `mean(replicates)-observed_point`, skewness, and covariance conditioning if feasible. The process bootstrap conditions on estimated S; it does not propagate uncertainty in S, pooled CV, allometry, offsets, or missing-year selection. Normalized-exponential count weights are empirical Bayesian-bootstrap weights, not a fully parametric count model. With one count unit or all observed zeros, their count uncertainty degenerates. Describe the combined direct sampling simulation accurately.

`biomass_point` supplies empirical direct points (200), while direct sampling masses are independently reconstructed from `samples` (185). Add a direct-row identity check that `biomass_point == density * raw matched mean mass` and zero rows have zero biomass. The code asserts sample raw mean identity but does not itself assert this cross-input point identity. No mismatch was seen in the synthetic fixtures; current real-data equality is not certified by this review.

### Remaining edge cases outside the reported full-run cohort

- The finite-year functional uses a pseudoinverse without a minimum-year guard (177). One supported year produces an artificial slope exactly zero with a degenerate interval; reproduced independently. Require at least two distinct years or emit an unavailable slope. This does not affect the reported current multi-year histories.
- Zero/one total annual observations still break `np.cov(Y)`/matrix indexing (169–172). Current 44+ annual observations avoid this, but subset/failed-support runs should have explicit handling.
- Unsupported mass guards do not explicitly reject positive infinity, negative/nonfinite density, or inconsistent zero-count flags. Add finite/nonnegative checks as input hygiene; there is no evidence these occur in current inputs.
- Process variance optimization has an artificial upper bound and only reports the lower boundary. Add an upper-bound warning/expansion and covariance condition checks before trusting extreme cases; no upper-bound failure was demonstrated here.
- Record observed first/last year and preferably the exact year list in the finite-year slope output. Current `n_years` alone does not identify the finite-year target.
- Names `annual_observation_covariance` and `observation_sd` still describe latent draw uncertainty in the first branch. Clarify labels/captions so downstream code cannot mistake them for the direct measurement covariance. The revised code itself does not make that mistake.

## Checks actually executed

`test_revised.py`: importlib import of this revised guarded module; synthetic complete `run()` with one protected and three reference sites, six years, two zones, measured/borrowed/unsupported cases, 40 draws, two parameter refits per model and ten process bootstrap draws per rho. It includes three held-out-reference fits, two full-data size models, 48 count reconstruction checks, five latent annual estimates and three process scenarios on four direct years. Checks slope arithmetic, support separation, first-failure retry provenance, rank/df guards, empty predictions, canonical missing-row/unknown-mass rejection, and lognormal mean centering. Results: `results.json`; synthetic artifacts: `pipeline/`.

`test_failure_paths.py`: independently forces exhausted size-bootstrap retry and process convergence failure; both stop without a final policy in fresh directories. Results: `failure_results.json`. Partial failure artifacts are intentional test evidence, not successful outputs.

No production ecological output or biological claim was validated by these bounded tests. The old combined-interval objection is resolved in this code; the remaining direct process result is explicitly a small-sample conditional working-model diagnostic, separate from the supported latent trajectory analysis.

## Concurrent update rechecked before completion

The parent added direct sampling covariance export and model-source hashing while this review was running. I inspected the exact diff, refreshed `reviewed_model.py` to the current 230-line version (SHA256 `8dda9d210c3d0ad25dfd9272255c5a263619e1f7cbce23d69c64098abc9ad338`), and reran both synthetic suites successfully against it. Exported direct covariance diagonal also matches the direct annual sampling SD squared. The first-read snapshot is retained as `reviewed_model_first_read.py`; earlier line references refer to that file (the final additions shift later lines by 2–4). No inference/failure-policy changes occurred in this final diff. Note that hashing source at run completion can misidentify executed code if the source changes during a run: record the launch-time hash or execute an immutable snapshot for the full run.
