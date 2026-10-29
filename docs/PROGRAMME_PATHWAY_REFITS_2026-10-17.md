# Programme-separated pathways on matched coverage

Module `code/R/42_programme_pathway_refits.R` answers the coauthor question about
averaging programmes within MPAs. It is a secondary comparison, not a replacement
for Table 3, a test of survey-method causation, or a change to primary membership.

## Comparison

For each of the six primary endpoint pairs, consider KFM/PISCO, LTER/PISCO and
KFM/LTER separately. Retain only post-designation calendar years with finite
predictor and outcome log MPA/reference ratios in both programmes at the same
MPA. Require at least four common years per MPA. All programmes and their average
then use identical MPAs and years within a comparison; different links can have
different coverage. Years need not be consecutive. Existing harmonised ratios
and their upstream zero-offset conventions are unchanged; no new offset or
imputation is introduced.

Two scopes are kept separate:

- `retained_primary_histories`: both endpoints in each programme must have a
  retained primary history. Matching is performed after applying that rule.
- `available_diver_histories`: available harmonised annual records at MPAs in
  the primary synthesis, including histories without a retained primary effect.
  This broader diagnostic scope does not promote those records into the primary
  synthesis. It also does not recover programmes or histories absent from the
  harmonised input. Landsat and expanded reconstructed lobster biomass are absent.

Within each programme, take each MPA's mean annual post-designation log ratio.
For the average, first average the two programmes' annual log ratios, then take
the MPA mean. Estimate annual-mean variances with the existing calendar-gap HAC
helper at three-year lags, and repeat at one-year lags. Averaging before the
variance calculation retains observed contemporaneous cofluctuation between the
two programme records, rather than assuming their variances can be divided as
though independent. It is not full propagation of survey or spatial uncertainty.

Fit REML meta-regressions with Knapp-Hartung intervals and outcome variances,
requiring four MPAs and positive outcome variances. Predictor uncertainty is
retained in output but not fitted. No route-selection test chooses between levels
and trends: every estimate here is a post-designation mean level. These quantities
therefore differ from the mixture of primary selected-model contrasts.

Use a six-link BH family separately for each programme/average, source pair,
scope and HAC lag, including unestimable links in the planned family size. Repeat
the fits after each MPA and each connected set of verified declared reference
identities is omitted. Omission BH families are separate six-link families within
each omission. Single-MPA and reference-group omissions can be identical and are
not independent evidence. Reference identities do not establish physically
identical transects or complete overlap mapping for non-primary histories.

## Findings

On the broader matched KFM/PISCO scope, Anacapa, Gull Island, Harris Point,
Scorpion and South Point provide five MPAs and 76 paired MPA-years for each of
four density-based links. At the three-year HAC lag, slopes are:

| Predictor to outcome | KFM slope [95% CI] | PISCO slope [95% CI] | Equal-programme mean slope [95% CI] |
| --- | --- | --- | --- |
| Purple-urchin density to kelp biomass | -0.268 [-2.071, 1.536] | -0.784 [-2.844, 1.277] | -0.363 [-2.128, 1.402] |
| Red-urchin density to kelp biomass | -0.364 [-2.083, 1.355] | -0.757 [-2.690, 1.177] | -0.506 [-1.970, 0.958] |
| Lobster density to purple-urchin density | -0.382 [-2.216, 1.451] | -0.432 [-1.722, 0.858] | -0.461 [-2.006, 1.083] |
| Sheephead density to purple-urchin density | -0.725 [-10.349, 8.899] | 0.483 [-1.530, 2.496] | 1.118 [-4.434, 6.671] |

Every full-cohort interval includes zero, at both HAC lags. Opposite sheephead
slope signs are descriptive, not a statistically established programme difference.
The small overlap does not establish equivalence or demonstrate that pooling
could not obscure a cascade. It also cannot isolate survey protocol from the
different reefs, references, detection or within-year sampling used by programmes.

Strict primary-history matching supports only the red-urchin-density/kelp-biomass
comparison at four KFM/PISCO MPAs (61 common MPA-years). Its KFM, PISCO and mean
intervals include zero. The other density links have two MPAs. KFM has no predator
biomass records; LTER/PISCO overlap supplies at most two MPAs. These sparse cases
are explicitly unestimable, not zero-effect findings.

Omissions are site-sensitive. Removing Anacapa yields a negative KFM
red-urchin/kelp slope with a nominal interval excluding zero at lag three
(-2.583 [-4.116, -1.050], P = 0.0185), but it does not pass the six-link correction.
Removing Harris Point instead makes that slope positive. These are sensitivity
results, not an additional confirmatory cascade test.

## Rebuild and audit

```
Rscript --vanilla code/R/42_programme_pathway_refits.R
Rscript --vanilla code/tests/test_programme_pathway_refits.R
```

The module is wired after module 41 in `run_all.R`. It records input hashes,
package versions, explicit fit/failure statuses, model cache, annual membership,
per-MPA estimates/variances, fitted slopes and diagnostics, omitted sets, coverage
and descriptive slope differences in `tables/table_s_programme_pathway_*`.
Primary effects, Tables 2/3 and main figures are protected by before/after hashes.
Tests cover endpoint/year matching, broader versus retained membership, identical
programme covariance, sparse/invalid fits, BH family sizes and numerical model
reproduction. No raw downloads are required.
