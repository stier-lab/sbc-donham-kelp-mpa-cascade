# Concrete joint implementation review — 2026-11-02

**No algebraic blocker found in the current within-history joint likelihood.** The independent tests verify signed covariance, aggregation derivatives, shared latent/physical-unit reuse, GLS estimates and inverse, REML objective/gradient, conditional latent predictions, input guards and optimizer failure handling. This is implementation verification, not certification of interval coverage or scientific transportability.

Current reviewed source: `code/tools/joint_borrowed_biomass.py`, SHA256 `970bfd18683faae175cf4edf48711878ff9331aabb24191c3abee1c631b0e864`.

## Reusable test artifact

```
python3 -B /tmp/donham_joint_review_20261102/test_joint.py \
  --repo /Users/adrianstier/sbc-donham-kelp-mpa-cascade
```

Latest run: **19 tests, 16 passed, 3 skipped, no expected failures**. Log: `test_results.txt` in this directory. Skips are the completed-main-output checks because the parent rerun had removed `completed.json`. These tests automatically run when a completed main output is available; they use its configured bootstrap count and seed. They check for output mutation during verification and source mutation during the suite. No production runner is invoked and bytecode writes are disabled.

## Independent checks executed

- Each size covariance component equals an independently constructed signed incidence product B B'. The observed/annual cross-block is negative, and target propagation uses the same physical random effects. Components are PSD to numerical tolerance.
- Analytic Jacobians match central finite differences, including unequal visit/zone/site replication, zero density and offset/no-offset cases. Exact aggregation is independently reconstructed by successive grouped means.
- Whitened QR/least-squares calculations reproduce beta, beta covariance and REML objective. The saved inverse reproduces the identity; numerical objective derivatives reproduce the analytic gradient.
- An augmented universal-kriging linear system reproduces `latent()` means and covariance. Re-requesting the same latent coordinate gives identical predictions with zero variance for their difference, without adding an observation.
- Independently reconstructed physical-unit exponential weights reproduce count replicates. Reused annual aggregates share draws exactly. Subtracting a fixed expansion constant does not change their covariance.
- Duplicate size observations, duplicate visits, duplicate years, invalid rho, unsupported biomass, nonfinite positive-count mass, missing count units, rank-deficient X and insufficient residual df are rejected. Single-site protected site variance is disabled.
- Forced failure of both L-BFGS-B attempts invokes a real SLSQP optimization; the returned fit matches independent GLS and records three optimizer attempts. Forced failure of all three attempts raises RuntimeError. No original-fit substitution occurs.
- The preserved smoke NPZ agrees with an independent Gaussian solve.

Earlier, before the parent restarted the output directory, all **18 original saved fits** reproduced under independent Gaussian solves, and the **3,582 original bootstrap rows** reproduced their reported studentized interval endpoints. This applies to the initial 199-bootstrap run, not the active 999-bootstrap rerun. A subsequent mixed-generation interval mismatch was traced to concurrent output replacement while the completion marker was absent; the test now guards against that situation. Real-data rebuilds and all 45 annual Jacobians also passed during that intermediate pass, but that was not a completed-run certification.

## Findings and disposition

1. **Resolved:** initial code regenerated count covariance for each rho. Current code gives `build()` the same per-series count RNG seed across rho, while keeping bootstrap RNG streams separate. The saved-output equality test is now a required pass, not an expected failure; it awaits completed current outputs.
2. **Resolved:** duplicate size/visit/year, rho and unsupported-biomass guards now pass required tests. The earlier test that deliberately duplicated a visit was replaced by a latent-coordinate reuse test, since duplicated visit observations should be rejected.
3. **Resolved within tested scope:** optimizer retries are genuine refits, including SLSQP, and exhaustion stops. This does not prove a global optimum; parent identification/profile diagnostics remain necessary. Recording optimizer method/messages and projected gradient alongside attempt count would strengthen provenance but is not a demonstrated blocker.
4. **Scope restriction:** `run()` still fits the six histories separately. It jointly models size and annual count contributions within each history, not a single six-history likelihood preserving cross-species count covariance. Individual-history marginal intervals can be used under that scope. Do not treat histories as independent in pooled or cross-history contrasts without carrying the omitted covariance.
5. **Remaining inference limits:** fixed sampling covariance/size pooling/allometry/offset/expansion, conditional count–size independence, first-order aggregation and donor transportability remain assumptions. Parent is handling information/rank profiles, donor influence and sparse-design simulation. Twenty simulated datasets per setting are exploratory calibration evidence, not a precise nominal-coverage certificate.

Repository and activity records were not edited. Only scratch tests, their log and this review were written. No external messages sent.
