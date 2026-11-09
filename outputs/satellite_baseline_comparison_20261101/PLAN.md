# Satellite baseline comparison — review analysis, November 1, 2026

Scope fixed before fitting: the 16 Landsat histories in the archived 145-history primary cohort. No candidate is selected by response sign or baseline significance. Four other harmonized satellite candidates remain outside this requested comparison; this analysis does not validate their historical exclusions.

For each pair fit ordinary least squares to annual log MPA/reference kelp ratios:
- Stable baseline: intercept + post-designation indicator + elapsed post-designation years.
- Changing baseline: the same terms plus calendar years centered at designation.
Both allow a level change and a subsequent slope change, and use identical finite observed years. These are new diagnostic models, not exact replicas of the existing primary linear pBACIPS model (which has no level-change term). Preserve existing designation-year Before labels for satellite data, including partial implementation years; do not quietly change annual timing.

Primary diagnostic quantity: average fitted departure during observed After years from continuation of the fitted Before trajectory. Also export departure in the last observed After year, level change and slope change. Neither quantity is automatically a causal effect or directly interchangeable with the primary pooled effect. Extrapolating a linear baseline for decades is a substantive assumption.

Uncertainty: calendar-gap Bartlett HAC covariance, lag 3 years, finite-sample n/(n-p) correction and t intervals with residual degrees of freedom. Export lag 1 and 5 sensitivity. No interpolated response observations. Individual pointwise intervals are exploratory; report BH-adjusted tests within each 16-site/model/window/lag family as well. No pooled inferential interval, because spatial/shared-reference dependence is not resolved.

Fit full available baseline, then the final 10 and 5 calendar years through designation, retaining all observed After years. Both models always use identical years within each window. Minimum three observed baseline and three After years, unique pair/year, finite ratios, full-rank design. Baseline slopes do not enter eligibility. Evaluate baseline prediction by withholding the last three observed baseline years from the full baseline and predicting them from the remaining baseline (flat versus linear). This tests near-term prediction only, not the long-term counterfactual.

Export eligibility, exact input rows, estimates/diagnostics, predictive checks and trajectories. Independently verify regression coefficients, contrasts and HAC matrices; fingerprint canonical primary outputs and manuscript/supplement before and after. User reviews results before any manuscript or primary-cohort change.
