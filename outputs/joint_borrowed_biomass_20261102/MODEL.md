# Joint working model

This model asks how the protected/reference biomass ratio changes through observed calendar years, allowing the animal sizes used to calculate biomass to be uncertain. It is not a before–after estimate of protection effects.

For a measured size sample, let y be log(mean individual mass) plus half its estimated log sampling variance. Under the working lognormal mean approximation, y measures the latent log arithmetic mean mass z. Size sampling variance is derived from size-frequency records and partially pooled within taxon, comparison and protected/reference role; singleton variance therefore comes from the pool, not an invented within-sample variance.

The size model is z = X_size alpha + u. It includes a zone contrast and calendar trend, with shared site, site-year and visit effects. A single measured protected station cannot identify site variance separately from its intercept, so that component is omitted. Measured samples enter once. Targets at the same site, year or visit share the corresponding effects. Eligibility still comes from the audited donor ledger, but **prediction now uses the hierarchical size model rather than reproducing a particular donor's arithmetic mean**. All eligible measured samples within the comparison and role can inform that model. The temporal trend is linear; it does not encode a validated seasonal decay in similarity. This is a substantive borrowing assumption to review.

For each year, arithmetic biomass is aggregated by equal visits within zones, equal zones within stations and equal stations within role. Define h(z) as the log ratio of these role biomasses with the existing fixed offset. Around the empirical mass vector z0:

    h(z) approximately c + W z
    c = h(z0) - W z0

W is the exact first derivative at z0. A protected visit contributes positively and a reference visit negatively; zero observed counts have zero mass derivative. Offsets mean the derivative weights need not sum to plus/minus one. Counts have their own sampling covariance S_count, estimated by resampling physical survey units with masses fixed. It contains no size uncertainty.

The recovery process is h(z) = X_recovery beta + process variation. Combining it with the measured sizes gives a single observation model:

    [measured log size] = [ X_size     0 ] [alpha] + size effects + size sampling
    [c               ]   [-W X_target X_recovery] [beta ] - W size effects + process + count sampling

The covariance between the two observation blocks is negative G_size,target W'. This shares the size information rather than treating borrowed estimates as new independent measurements. No conditional posterior size covariance is inserted as a new measurement likelihood.

Variance components are estimated by REML; all size and recovery fixed effects are jointly profiled by GLS. The process covariance uses variance times rho raised to the calendar-year gap. Rho=0, 0.3 and 0.6 are reported scenarios. The same count covariance and year roster are used across them. Species/comparisons are fitted separately; this report does not estimate cross-taxon pooled inference.

Parametric bootstrap datasets simulate the entire joint vector, including its cross-block covariance, then refit all variance components and fixed effects. Studentized intervals use the bootstrap slope error divided by its refitted standard error. Basic bootstrap alternatives remain in the results table. Optimizer retries run actual alternative optimizations; exhausted failures stop the run. No failed replicate is replaced by an original estimate.

## Conditional limits

Sampling variances and their partial-pooling strength, count covariance, allometry, offsets and first-order expansion are conditional inputs. They are not newly estimated in each bootstrap. Counts and size sampling errors are assumed independent conditional on the model. Observed zero count visits are held at zero; detection error and unsampled rare animals are not modeled. Within-year visit effects are exchangeable beyond the shared site-year effect. Permanent transect identities across years have not been established.

Alternative pooling, random seed, strict timing eligibility, calendar correlation, expansion center, omitted count years and omitted size-donor sites are checked separately. Simulations under this working model do not validate its nonlinear approximation or transportability to protected sites. Independent Gaussian algebra tests validate implementation; coverage simulations evaluate a different question and may reveal inadequate intervals. Neither result licenses a total-uncertainty claim.
