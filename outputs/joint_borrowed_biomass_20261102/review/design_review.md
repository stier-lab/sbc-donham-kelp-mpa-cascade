# Independent joint biomass design review — 2026-11-02

**Verdict: the proposed stack is a coherent conditional, linearized Gaussian mixed model, provided the signed cross-covariance is retained. It resolves the latent-posterior-covariance reuse problem. It is not yet validated inference for borrowed biomass.** The main gates are covariance-component identifiability, nonlinear approximation error, and calibration of the complete fitting procedure. No new joint implementation was present in the inspected model file; this reviews the proposal against its predecessor and saved inputs.

## Correct likelihood and interpretation

Write z = X_t alpha + u_t, D = diag(v), S = fixed count-only sampling covariance, Q = temporal-process covariance, and G_ab = Cov(u_a,u_b), using one shared random-effect construction for observed and target sizes. Under independence of sampling errors and process/size random effects:

```
r = [y_s; c]
A = [ X_s       0
     -W X_t     X_recovery ]
gamma = [alpha; beta]

V = [ G_ss + D          -G_st W'
     -W G_ts            W G_tt W' + Q + S ]
```

This follows directly from c + Wz = X_recovery beta + process + count_error. The negative cross-block is essential; making the stack block diagonal discards shared size information. Construct V through signed random-effect incidence matrices plus measurement covariance, rather than independently assembled blocks. Do not add posterior target-size covariance, alpha's estimated covariance, or a second donor residual to V: joint estimation already accounts for these pathways.

Estimate covariance parameters by REML, profiling **both alpha and beta together**; then estimate fixed effects by GLS. “Jointly REML beta” is shorthand, not a literal optimization of beta in the restricted likelihood. The criterion includes log|V|, log|A'V^-1 A| and the GLS residual quadratic. The existing `fit_reml` accepts diagonal known variance only; the new implementation needs the full known block containing S. Fixed-effect comparisons require ML rather than comparing REML values across different A. Standard mixed-model estimation context: [Bates et al., 2015](https://www.jstatsoft.org/article/view/v067i01).

The +v/2 size correction is correct for the working model log(sample arithmetic mean mass) ~ N(log(true arithmetic mean mass) − v/2, v). It is moment matching, not an exact small-sample sampling law. Here v is sampling variance of the log sample mean, not biological log-mass variance. Transform individual sizes through the allometry before averaging, as the current code does.

For role r, define B_r(z) = sum_i a_i d_i exp(z_i), with a_i reproducing equal visits within zone, zones within site, and sites within role. Then

```
h(z) = log(B_mpa(z)+offset) − log(B_reference(z)+offset)
W_i = +a_i d_i exp(z0_i)/(B_mpa(z0)+offset)       [mpa]
W_i = −a_i d_i exp(z0_i)/(B_reference(z0)+offset) [reference]
c = h(z0) − W z0
```

The offset is fixed, not differentiated as mass. Zero observed density contributes zero derivative; exclude zero-mass placeholders from logmass coordinates. Unsupported positive counts remain unsupported. W's signed row sum generally is **not zero** with the offset; do not normalize it to zero or force role sums to ±1.

## Findings requiring action

1. **Remove exact variance aliases; quantify weak identification.** All three Abalone protected size groups have one site. Its random site effect has the same joint incidence as that group's fixed intercept, including the −W target block, so joint REML cannot recover its variance. Retain the current zero/omitted component treatment, or specify an explicit cross-group variance-sharing assumption; borrowed target uses do not create independent sites. Point Vicente has only two protected sites. Numerous existing variance estimates are on boundaries. This is weak evidence for separately estimating every site/site-year/event variance, even when algebraic rank is adequate.

   Check rank and condition of the full A and covariance kernels after projecting out A. For fitted V, inspect the variance-component information matrix I_jk = tr(P K_j P K_k)/2, where P = V^-1 − V^-1 A(A'V^-1 A)^-1 A'V^-1. A full-rank fixed-effect design alone is insufficient. Profile process variance against size site-year/event variance. At rho=0 these can be nearly indistinguishable in the annual block; shared size observations can separate them only to the extent those observations supply replication. Do not fit an unstructured six-history process covariance to these short histories by default.

2. **Fixing expansion is a legitimate approximation, not a cure for generated-data dependence.** Empirical density and donor-based mass affect c, W and z0. Conditional fitting treats these realized quantities as fixed. It therefore ignores their estimation dependence and higher-order count–mass interactions. Compute count-only S from resampled h(d*, z0), keeping the ORIGINAL W and z0 fixed; subtracting the fixed Wz0 leaves the same covariance. Recomputing W in each count replicate estimates a different statistic's covariance. Never use the existing direct sampling covariance: it already includes size sampling. Never use annual latent covariance.

   Quantify log-ratio sampling bias and skewness at fixed mass. Covariance alone does not correct a biased log observation. For mass curvature, compare exact h(z) with c+Wz across plausible joint target deviations and alternative expansion points. The leading error is 0.5 delta' H delta; its expectation need not be zero. No second +v/2 correction fixes this aggregation curvature. Offsets and sparse/zero roles can make it substantial.

3. **Borrowing is now a model assumption, not merely a donor lookup.** The current size design is intercept + MID + linear calendar trend, separately by species/MPA/role. Shared site-year effects correlate all same-year dates; the event kernel shares only identical site/year/zone/date cells. This does not enforce the empirical donor gap limits, seasonal effects, or specific donor weights. A donor mass used as z0 anchors an expansion, not a second observation or a prior that forces the fitted target mass to equal that donor estimate. Preserve the support mask, expose which size observations inform each target, and distinguish eligibility policy from model borrowing. Sparse protected validation remains unresolved. A recovery model can feed back into size estimates in a joint likelihood; compare those shifts with size-only fits rather than assuming feedback is harmless.

4. **Sampling independence needs an explicit scope.** The displayed V assumes Cov(size sampling error, count sampling error)=0 and independence of temporal process from size random effects. Shared physical transect resampling preserves count dependence across taxa and reused units, but does not establish either independence assumption. If count/size observations share units or size affects detectability, evaluate a joint-unit resampling sensitivity where identifiers permit; otherwise state the unresolved dependence. If R_sc is nonzero, the cross-block becomes −G_st W' + R_sc. Preserve cross-history entries of S when fitting the six histories jointly; separate series fits are not equivalent to that joint likelihood.

5. **Sampling support and process claims remain conditional.** The present normalized-exponential count weights are empirical Bayesian-bootstrap weights, not ordinary multinomial resampling or a fitted count distribution. Name the scheme and check its variance scale against the intended sampling design. A one-transect or all-zero stratum yields zero empirical count variance; observed zero is not proof of a zero population mean. Pooling size dispersion does not repair absent count replication. Keep count method averaging, physical identity, and reuse across taxa; do not infer repeated-transect identity across years without evidence. Allometry, offsets, pooled sampling variance, missingness/eligibility and seasonal/depth representativeness remain outside the stated conditional uncertainty.

## Bootstrap and required acceptance tests

For each fixed rho scenario, simulate complete r* ~ N(A gamma_hat, V_hat) and refit the **entire** joint model, including all estimable variance components and both fixed-effect vectors. Equivalently simulate shared size random effects once, then independent process and sampling errors with their specified dependence. Do not simulate y_s* and c* independently, condition on original size observations, draw beta additionally, or re-estimate S/v/z0/W in this conditional bootstrap. Hold the exact year/site/support roster fixed across rho=0, .3, .6; use rho^abs(calendar-year difference), including gaps. A stationary variance q gives Q_ij=q rho^gap; distinguish q from AR innovation variance.

Specify the confidence interval construction (e.g. basic or studentized, with its actual refit statistic), rather than labeling any simulated slope range a CI. Report failed refits, retry rules, boundary frequencies, upper-bound hits and Monte Carlo stability. Valid zero variance estimates are not failed fits. Do not replace failures with the original fit or silently discard them. A fitted-model bootstrap checks conditional uncertainty; it cannot certify its own model assumptions. Boundary and weak-variance cases particularly require coverage experiments; see [lme4 singular-fit guidance](https://lme4.github.io/lme4/reference/isSingular.html).

Before interpreting recovery intervals, require:

- **Algebra/identity:** analytic W versus finite differences; h(z0)=c+Wz0; unequal visit/zone/site counts; offset-dominated and zero-role cases; correct signs and physical deduplication; direct-only and zero-size-uncertainty limiting cases.
- **Gaussian reference:** independent dense implementation of A,V, REML objective, GLS estimates and covariance; PSD and cross-block sign checks; simulated empirical joint covariance matches V. Compare joint fitting with the mathematically equivalent conditional Gaussian factorization.
- **Identification/optimization:** full projected-kernel diagnostics above, multiple starts, boundary fits, profile intervals, scaling and jitter sensitivity. Cholesky success or optimizer success alone is insufficient.
- **Calibration:** repeated datasets under the linearized model, preserving actual sparse design, to measure beta bias, 95% coverage, interval width and failure rate for each history/rho. Include small and zero process/site/event variances and high size uncertainty. Report Monte Carlo uncertainty on coverage.
- **Approximation robustness:** a separate nonlinear observation simulation/reconstruction that recomputes the empirical expansion and sampling covariance, plus alternative expansion points, to measure errors intentionally excluded by the conditional bootstrap. Compare exact and linearized annual responses and slope contrasts relative to their reported SEs.
- **Influence/model fit:** leave-one-year and leave-one-donor-site analyses; whitened residuals and time patterns in both blocks; size-only versus joint size fits; pooling-strength sensitivity. Assess nonlinear size trends/seasonality without silently changing the admitted years. Report each history's slope/interval across all three rho scenarios, not the most favorable one.

The estimand is the modeled trend in the offset biomass log ratio on the admitted survey design. Changing sites/visits or missing years may change its ecological meaning. Joint uncertainty propagation alone does not establish an MPA causal effect, protected-site transportability, or a population-wide trend.

## Evidence actually checked

Read-only inspection of `code/tools/biomass_observation_model.py` (249 lines; SHA256 `ba3c595e0ce94f5277c3f918436dbb223b2a6ba2bb8ff506924e2eecafd658c5`), the next-steps README, model-review/disposition records, policy, size diagnostics, size registry, visit support, annual slope roster and fixed offsets. Existing test success statements were read, not rerun or independently endorsed here.

Independent in-memory checks, with no repository writes:

- Size registry: 682 rows and 682 unique physical sample IDs.
- Projected size-kernel ranks: 2 for each single-site Abalone protected group (site kernel norm effectively zero); 3 for each other group. These are size-only structural checks, not proof of identification of the proposed full joint model.
- Reconstructed aggregation derivatives for all 45 post-protection years under the existing seven-day diagnostic support: maximum absolute finite-difference error 4.875e-11. Signed row sums ranged −0.97695 to 0.98830, confirming material offset dependence.
- Existing post-year counts are Abalone red/lobster/purple 8/10/7 and Vicente red/lobster/purple 5/9/6. These short histories warrant the calibration tests above. Seven-day eligibility is a diagnostic policy, not a newly approved primary cutoff.

No joint model was fitted and no bootstrap coverage experiment was performed in this review. Only this `/tmp` report was written; no activity session, repository edits, primary promotion or external messages.
