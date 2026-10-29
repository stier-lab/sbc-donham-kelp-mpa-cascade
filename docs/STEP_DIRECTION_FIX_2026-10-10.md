# Step Contrast Direction: Independent Confirmation And Focused Fix

## Confirmed Bug

Read-only verification of all six retained primary Step histories confirmed the actual `add_step_effect_size()` output, not merely a manual reconstruction. Their BA columns are character; `lm(lnDiff ~ BA)` orders levels `After, Before`. Consequently `pairs(em, reverse = TRUE)` explicitly labels its contrast **Before - After**. The helper reproduces all six stored primary means, each exactly the negative of the direct mean(After) - mean(Before) contrast. This is a pre-existing direction bug, separate from the COUNT correction.

## Preserved Baselines

The original phase remains in `local_archive/lter_count_refit_2026-10-10/before`. The parent preserved the full COUNT-only phase in `local_archive/lter_count_refit_2026-10-10/count_only`; immediately before this helper edit, all 298 saved result/cache/figure files matched the independently audited live outputs byte-for-byte. See `docs/LTER_COUNT_REFIT_2026-10-10.md` for COUNT-only numerical results. None of those results is superseded silently here.

## Expected Six Reversals

The following are verified helper outputs on the actual retained histories. CI is the stored **half-width**, not a lower or upper endpoint. SE, SD, CI half-width, N and retained model are unchanged by the contrast-direction fix.

| MPA | Taxon | Source | Response | Original and COUNT-only mean | Correct After - Before | SE | CI half-width | N |
|---|---|---|---|---:|---:|---:|---:|---:|
| Gull Island SMR | M. franciscanus | KFM | Den | 0.430712882940 | -0.430712882940 | 0.365251878038 | 0.742281123661 | 36 |
| Naples SMCA | M. franciscanus | LTER | Den | 1.189470103710 | -1.189470103710 | 0.395236662858 | 0.821939636024 | 23 |
| Santa Barbara Island SMR | S. purpuratus | KFM | Den | 1.992445112450 | -1.992445112450 | 0.303334766379 | 0.613552477687 | 41 |
| Scorpion SMR | S. pulcher | KFM | Den | -1.543103140442 | 1.543103140442 | 0.329899948959 | 0.669068107407 | 38 |
| Scorpion SMR | M. franciscanus | KFM | Den | -0.622103592406 | 0.622103592406 | 0.205535197728 | 0.416084254758 | 40 |
| Scorpion SMR | M. franciscanus | KFM | Bio | -1.253036353820 | 1.253036353820 | 0.270680502420 | 0.549510633985 | 37 |

## Authorized Fix And Tests

Only the contrast construction in the Step helper in `code/R/08_effect_sizes.R` changed: it now matches `c(After = 1, Before = -1)` to the actual BA labels in the emmeans grid, then requests a named `After - Before` contrast. The model formula, annual histories, clocks, inclusion rules, diagnostics, model selection and uncertainty calculations are untouched. The helper does not depend on alphabetical order or a caller's factor-level order.

`Rscript --vanilla code/tests/test_step_effect_direction.R` passes. Tests extract only the original/current helper definitions, without sourcing module 08, and cover character BA, both factor-level orders, reversed rows, direct arithmetic means, independent pooled-model SE/CI calculations, and all six actual retained histories. Fingerprints confirm no production annual or primary CSV writes. The parent owns all full pipeline reruns.

## Combined Refit Gate

After the parent's full rerun, expect **six Step means plus the Campus sheephead/LTER/Bio numerical change relative to the original primary set**. Relative to COUNT-only, expect exactly the six Step mean sign reversals, with all other primary fields unchanged. The Step fix does not change the raw/RR inputs or any model membership, and should not be used to admit additional histories.

The pooled biomass and density fits, their nine-test FDR adjustment, pathways and primary-effect-dependent sensitivities must be recomputed by the parent; their combined-refit numerical consequences have not yet been audited. In particular, the COUNT-only statements that density coefficients and q values are unchanged do **not** describe the combined refit. New module-32 caches and display fingerprints must reflect the corrected primary table. No combined-refit numbers are asserted in this pre-rerun record.
