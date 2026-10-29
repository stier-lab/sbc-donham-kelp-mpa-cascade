# Biomass reconciliation and singleton correction

## Biological conclusions

The corrected primary synthesis retains all 145 histories, all nine pooled
directions and their primary adjusted classifications. Five individual
biomass contrasts change: PISCO lobster at Abalone Cove, Matlahuayl and Point
Vicente; PISCO red urchins at Campus Point; and LTER red urchins at Naples.
The other 140 individual contrasts retain their estimates and SEs.

Two qualifications change. The negative lobster-purple-urchin biomass
relationship now has a nominal interval excluding zero (slope -0.5616,
95% CI -1.0577 to -0.0655; P = 0.03965). Its four MPAs provide only two
residual degrees of freedom; an independently enumerated 24-label test gives
P = 0.125. Campus Point and Naples share references. This is an association,
not evidence of a complete predator-mediated cascade.

Purple-urchin biomass remains inconclusive in the primary pooled fit
(adjusted P = 0.0940), but is negative under automatic and three-year
calendar-HAC specifications (adjusted P = 0.0460 and 0.0371). The one-year
specification remains inconclusive. All pooled directions persist; the earlier
claim that every adjusted conclusion survives annual-dependence choices must
not be used with these corrected inputs. These pooled contrasts are distinct
from the separate post-designation temporal slopes in Figure 4.

Modern diver-kelp resistance/recovery and historical satellite comparisons
are numerically unchanged. The correction does not resolve reference habitat,
unverified shared observations, historical protection or satellite masks.

## Error and its scope

The inherited bootstrap expanded size-class counts into a vector `a`, then
called `sample(a, n, replace = TRUE)`. In R a single positive numeric value
is interpreted as an integer upper bound, not a one-element population.
A sole 60 mm lobster could therefore be resampled as fictitious lengths
1-60 mm. The correction uses `a[sample.int(length(a), n, replace = TRUE)]`.
All draws from a one-animal pool now retain its measured length.

The source cache repair replaces only deterministic singleton totals; all
nonsingleton cached draws, record keys and missingness are preserved. It
changes 14 lobster, 50 PISCO-urchin and 26 LTER-urchin cached biomass cells.
Many repaired cells lie outside exported or primary histories. No KFM
singleton exposure was found, and no retained primary purple-urchin history
has a positive-count singleton exposure. One source row with count greater
than one is not a singleton; neither is a constant-size pool containing
several measured animals.

The source-data promotions change six exported lobster and 17 exported
urchin annual/status biomass rows, followed by 31 lobster and 48 urchin
annual-ratio rows. They preserve densities, kelp, sheephead, Landsat,
membership and all other non-target record bytes. Earlier monthly-kelp and
sheephead markers are retained as parents of the two scoped correction
markers. Routine input guards verify that chain rather than rewriting history.

## Remaining primary-secondary differences

The secondary size analyses formerly rounded the lobster allometry. They now
use the same SBC-LTER PAIN coefficients as primary processing:
`mass(g) = 0.001352821 * CL(mm)^2.913963113`. The source is the version-3
coefficient table, not Barsky; source DOI
https://doi.org/10.6073/pasta/0fe9233dabe35df5d61fb3b07f8fb51e.

The original analysis-directory bootstrap caches differed from the caches
that actually produced the harmonized data. Corrected authoritative source
caches were therefore copied into the analysis directory after preserving the
older copies. Do not mix old cache draws with the checked current release.

Module 43 reconstructs the primary VRG hierarchy: cached zone totals divided
by their transect denominator, equal zone means within sites, equal site means
within MPA/status, then division by 60 m2. It reproduces overlapping primary
biomass to a maximum absolute difference of 4.3e-14 g/m2. Using the same
full-year size pools, deterministic expected masses agree with that hierarchy
to floating-point tolerance. Their ratios to the stochastic primary estimates
range from 0.9685 to 1.0144 across 47 positive overlapping records, consistent
with finite Monte Carlo draws. Degenerate pools now have zero arithmetic error;
that does not establish biological certainty about population size structure.

Secondary size comparisons deliberately use May-October size pools. With
unrounded coefficients their ratios to corrected primary biomass range from
0.8498 to 1.5680 in those overlapping records. Aligning the size-pool season
removes these larger discrepancies. Expanded UCSB records also reconstruct
mass from pooled sizes rather than primary direct transect lengths. Expanded
MPAs/years remain secondary; reconciliation does not authorize new primary
membership or establish equivalence of sampling methods.

## Reproduction and verification

- Source workspace: `code/R/correct_lobster_singleton_cache.R`,
  `code/R/correct_urchin_singleton_cache.R`, the fixed `bootstrap_biomass()`
  helper and singleton regression tests. Repair scripts require separate
  backups and refuse repeated destructive application.
- Run source harmonization into `DONHAM_HARMONIZED_OUTPUT`, then use
  `code/promote_lobster_singletons.py` (lobster and `--urchins` phases) against
  preserved baseline directories. The existing scoped-promotion machinery
  rejects changed keys, metadata and unrelated values.
- Refresh the reference registry after corrected primary contrasts are
  generated. Analysis `run_all.R` includes the updated size checks and module
  43. The singleton helper is the correction for fresh source rebuilds;
  targeted cache repair preserves the previous nonsingleton realization.
- `code/audit_biomass_singleton_refit.py` compares current inputs/results
  against the preserved initial checkpoint. Its JSON report lists every
  changed input and individual contrast, pooled classifications and the two
  HAC detection changes.
- `code/tests/test_lobster_biomass_reconciliation.R` checks full-precision
  masses, degenerate-pool arithmetic and primary/secondary reconstruction.
- Independent source investigations are preserved in
  `LOBSTER_SINGLETON_INDEPENDENT_AUDIT_2026-10-18.md` and
  `SINGLETON_NONLOBSTER_AUDIT_2026-10-18.md`. These are dated source audits,
  not substitutes for the final refit verification.

Frozen pre-correction material and cell-level cache ledgers remain in ignored
`local_archive/lobster_reconciliation_2026-10-18/`. New manuscript/SI prose
and coauthor correspondence belong in the separate manuscript workspace.
