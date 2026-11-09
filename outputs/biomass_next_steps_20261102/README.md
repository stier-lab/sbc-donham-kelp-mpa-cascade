> Earlier review stage. The [joint borrowed-biomass model review](../joint_borrowed_biomass_20261102/README.md) now extends the uncertainty work below; retain this report as dated evidence.

# Implemented biomass and coverage follow-up

November 2, 2026. These are review analyses, not changes to the manuscript or primary release. Researcher exclusions remain intact.

## What changed

**Density now has its own baseline roster.** Six prior comparison-years and 33 station-year uses are restored. All 63 post-protection density ratios remain unchanged. The current density input has 84 comparison-year rows; three earlier shared-zone ratios remain unknown, with their full-zone observations retained alongside. Point Vicente has four prior years with its full observed zones, but only three with shared INNER/MID coverage: ROCKY_POINT_S was DEEP-only in 2008. Retaining the count evidence does not certify a comparable before–after baseline.

**Size matching now keeps both attempted paths.** The ledger retains local-other-date and other-site candidates, failed guards, dates, sample counts, calibration sites/years, and stable physical keys. Missing sizes remain distinct from absent validation. Positive counts without supported sizes remain unknown biomass; no reference stations are dropped silently.

**The date boundary is tested consistently.** Zero- and three-day margins retain 44 of 63 post biomass comparison-years; seven and 14 days retain 45. Seven days restores Abalone red 2019 through Point Fermin. It also changes two Abalone purple 2019 visit records to local rather than cross-site donors. Red slope changes from −0.180509 to −0.188102 through coverage; purple changes from −0.098868 to −0.101968 through donor choice. The other four recovery slopes are unchanged. These are conditional temporal fits; no final primary tolerance is selected.

**Uncertainty is developed without recycling a smoothed estimate as new data.** Twelve species/MPA/role size models use measured means, explicit partially pooled sampling dispersion, site/site-year/event variation and the lognormal sampling-mean correction. Shared physical count units receive the same resampling weights across taxa. Forty covariance-parameter bootstrap refits per model feed 800 conditional draws. Support is not expanded: the seven-day scenario still has 45 post comparison-years.

The resulting slope-draw ranges describe uncertainty in change over the observed years. They are not confidence intervals for a general recovery process or protection effect. The rejected shortcut of using latent posterior covariance as independent measurement covariance is absent from generated results.

A separate observation-plus-process working model uses **unsmoothed empirical ratios and only directly measured sizes**. Only the two lobster histories have at least four complete post-protection years (8 Abalone, 5 Vicente). Sampling covariance is estimated from count-unit and measured-size replicates; process variance is then estimated conditional on that covariance. Calendar correlations 0, 0.3 and 0.6 are explicit scenarios. All six conditional bootstrap intervals span zero. The linked decomposition distinguishes coverage losses from GLS weighting changes; this subset must not replace the full density analysis.

## What still cannot be called resolved

- Original recent KFM sizes remain unavailable for 54 comparison-years across six biomass histories. Forty-five historical years have diagnostic replacements; one verified zero-reference year needs no size donor. Density remains retained. The exact sender/date inbox check found no new Scott reply; the request was already sent.
- Protected-site borrowing is not independently validated. Reference-site model checks do not establish transportability to protected stations. Conditional held-out size intervals cover roughly 88–98% of reference sample means, depending on comparison; observations within sites are dependent. Evaluation intervals condition on the held-out sample’s estimated dispersion; these are not outcome-blind predictive coverage estimates.
- Size-model parameter mixtures are an approximation, not an exact posterior. Pooling strengths 2 and 10 and a second random seed are compared. The pooled sampling variance itself, allometric coefficients, detection, count–size dependence and persistent survey-unit identity are not fully estimated.
- Fully joint process inference for histories containing borrowed biomass remains unfinished. Neither latent-slope ranges nor direct-only working intervals are total uncertainty. No additional animal histories were excluded to hide these limits.
- KFM kelp single-date visit variation and all-program numeric depth equivalence remain data limitations. Current available-reference assignments and program-specific methods remain; fixed-reference and depth alternatives are sensitivities.

## Current tables

- [Current density annual input](current_density_annual.csv), [density baseline](density/baseline_coverage_diagnostic.csv) and [restored station-years](density/recovered_station_year_uses.csv).
- [Detailed visit support](borrowing/visit_support.csv), [both donor attempts](borrowing/donor_attempts.csv), [timing tolerance](borrowing/summary.csv), [effects and decomposition](borrowing/effects.csv).
- [KFM year holds](kfm_year_holds.csv), [program rules](program_coverage_rules.csv), [implementation status](implementation_status.csv).
- [Size/count registry and diagnostics](uncertainty/size_model_diagnostics.csv), [validation summary](size_validation_summary.csv), [observed-year slope draws](uncertainty/observed_year_slope_draws.csv).
- [Direct-process admission](uncertainty/direct_process_admission.csv), [process results](uncertainty/direct_trend_scenarios.csv), [coverage versus estimator changes](direct_process_change_decomposition.csv), [pooling/seed checks](uncertainty_sensitivity.csv).
- [Current consolidated guide](../dataset_inventory_20261101/README.md).

## Reproduction and validation

Run `python 3 code/tools/run_biomass_next_steps.py` from the analysis repository. It rebuilds density, donor matching/tolerance, conditional uncertainty and sensitivity variants, verifies them, then updates the guide. Raw sources, prior dated scientific outputs and the manuscript are not overwritten. The live email check is a dated record and is not repeated by this runner.

Density has 16 independent source/aggregation/edge checks; borrowing has 5 regression tests plus complete-reference, zero/unknown, source-identity, calibration and hash assertions. Size covariance is checked against independent Gaussian precision algebra; visit identity, singleton pooling, fit provenance, PSD covariance and saved slope draws are tested. R independently reproduces all 6 direct-only GLS/REML fits. All 487 protected manuscript/canonical files match their original hashes. Validation and limitations are recorded separately.
