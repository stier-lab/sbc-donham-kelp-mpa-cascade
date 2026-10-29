# PISCO reference-control sensitivity

## Question

Bob Miller raised the possibility that designated PISCO reference reefs are poor
counterfactuals for their paired MPAs. If so, a response ratio can be biased even
when the underlying PISCO survey data are accurate.

## What the current harmonized data can test

`code/R/26_pisco_control_sensitivity.R` runs two dependence checks using the
effect-size snapshot used by the primary meta-analysis:

1. **Remove PISCO entirely.** The joint meta-analysis is refit after dropping
   all PISCO MPA--reference effects. This asks whether independent monitoring
   programs support the same result.
2. **Leave one PISCO pairing out.** Each PISCO MPA's PISCO effects are removed
   in turn. This asks whether one designated PISCO reference pairing drives a
   result. Any observations from other programs at the same MPA remain.

The outputs are `tables/table_s_pisco_control_source_removal.csv` and
`tables/table_s_pisco_control_leave_pair_out.csv`. The analysis is included in
`run_all.R` and is reproducible directly with:

```r
source("code/R/26_pisco_control_sensitivity.R")
```

## Results from the current effect-size snapshot

The current cohort has 145 effect sizes, 89 from PISCO. Use the generated tables for current estimates; prior prose summaries referred to an older cohort.

With PISCO removed, the kelp contrast remains positive (log contrast 1.617, nominal 95% CI 0.195 to 3.039) and purple-urchin density remains negative (-2.483, -4.010 to -0.957). Lobster density remains positive but its interval includes zero (1.003, -0.467 to 2.473). These are the table's nominal intervals, not a new multiple-test-adjusted decision set.

Sheephead biomass changes sign on the two remaining non-PISCO histories (-3.638, -6.352 to -0.925). That sparse, differently distributed comparison does not justify a claim of cross-program consistency. Source omission changes sites, coverage and precision as well as monitoring program. Inspect the per-pair deletion table separately rather than treating source removal as a test of reference suitability.

## What this cannot establish

These analyses test **dependence on PISCO references**, not whether a given
reference reef is truly well matched. The tracked harmonized files retain only
the paired MPA/reference ratio (`mpa`, `reference`) and not the individual
reference-site observations, depth, habitat, fishing pressure, or
transect-level replication. The current registry supplies declared identities
and PISCO survey footprints, but not an independent ecological matching audit. A zero/weak result after removing PISCO is not
evidence that its controls are poor; it may simply reflect loss of precision.

## Definitive follow-up with raw PISCO data

The appropriate next analysis requires the raw site- and transect-level PISCO
tables plus the MPA-to-reference pairing key. Before looking at post-protection
outcomes, create a pair-balance table for depth, substrate/kelp habitat,
latitude/region, exposure, baseline abundance, and pre-period trend. Then:

1. pre-specify match thresholds (or estimate propensity/matching weights);
2. refit the CI/pBACIPS effects using only well-balanced pairs and using
   inverse-distance/propensity weights;
3. report the primary, matched-only, and re-paired/alternative-reference
   estimates together; and
4. preserve all candidate controls and pair-selection decisions in an audit
   table.

That audit would directly address whether the controls are poor rather than
only whether the aggregate conclusion depends on them.
