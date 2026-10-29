# Naples Sheephead Biomass Inclusion

Adrian approved inclusion on 12 October 2026 after the earlier omission could not
be explained by its recorded before-trend rule. This is a membership decision,
not a change to fish conversion, offsets, raw observations or reference assignments.
The primary cohort now contains 145 histories: 108 CI and 37 pBACIPS. Models are
68 Mean, 68 Linear, six Step and three Sigmoid.

## Eligibility and Limits

Naples LTER sheephead biomass has 23 paired annual observations (2001-2023),
twelve before/designation-year and eleven subsequent observations. The corrected
before-trend P is 0.098899, above the stated 0.05 exclusion threshold. It does
not establish parallel trajectories. AICc gives the linear candidate the largest
weight (51.031%). The history retains the existing observation-index clock and
model-time-11-minus-0 contrast: -4.081200, SE 1.580579. It is a relative fitted
contrast, not an absolute biomass decline or demonstrated protection effect.
The last observed after-period index is 10, so the conventional horizon at 11
extends one index beyond observed support. The fitted-history audit records
this extrapolation explicitly; the atlas draws curves only over observed years.

Reference biomass is zero in seven before-period years and one subsequent year.
The existing half-minimum offset is retained. Linear residual dependence remains
visible (DW 1.4125, P 0.044). Neither the offset nor the failed alternatives were
changed to improve the result. Naples density and all other omissions remain unchanged.
The original exclusion's rationale is still unknown; the new inclusion is an
explicit approved decision, not a retrospective claim about that rationale.

## Refit and Independent Audit

The full pipeline completed with exit 0 in 6.63 minutes; log
`logs/pipeline_20261012_140404.log`. The exact production LTER biomass block was
then replayed to export the added eligibility/diagnostic metadata; it did not
rewrite primary result files. Optional modules 18 and 20 were refreshed too.
The registry was explicitly re-exported and now has 861 reference records.
The original 857 records retain their identities, footprints, exclusions and years;
four records belong to the newly admitted history.

`audit_naples_inclusion.py` independently reconstructs the linear contrast and SE,
verifies that all 144 preceding estimates, sample counts, models and designs are
unchanged, checks byte-preserved annual/satellite/site inputs, and reconstructs
the nine-test BH adjustment. Its report is `NAPLES_BIOMASS_INCLUSION_2026-10-12.json`.

The exponentiated pooled sheephead biomass contrast changes from 3.4366 to 3.3271
(95% CI 2.1436-5.1642; adjusted P 0.000005739). All nine pooled directions and
detection classifications remain. Joint biomass coefficients and their intervals
shift; density coefficients differ only at numerical optimizer tolerance, while
their BH values also reflect the changed biomass family.

The sheephead-biomass/purple-urchin-biomass slope changes from 0.2835 to 0.1395
(P 0.453735; nine MPAs) and remains inconclusive. The five other pathway
regressions and all six MPA memberships are unchanged. Modern kelp
resistance/recovery and historical Landsat comparisons are unchanged.
The matched-design subset now contains 60 source-specific histories.

Module 34 refits the same joint REML models with Naples omitted. Its exclusion
reproduces the preserved 144-history pooled coefficients and BH values within
1e-7. This is a consequence check, not a criterion chosen to retain a desirable
effect. It is exported in `table_s_naples_biomass_inclusion.csv`; history and
eligibility exports identify the sample and diagnostics. Existing Gull omission
scenarios also retain all nine primary directions/detection classifications.

The residual-test ledger has 49 detectable tests among 138 tested fits, median
DW 1.636. The separate nested-source variance-inflation reconstruction has median
rho 0.163 and median interval-width ratio 1.024, with no sign/detection flips.
The manuscript's old width ratio 1.021 was stale even relative to the preceding
allometry outputs (1.023); it is now refreshed to the current value. This remains
an exploratory reconstruction, not a temporally corrected primary fit.

## Reproduction

```sh
Rscript --vanilla code/R/run_all.R
Rscript --vanilla code/tests/test_naples_inclusion.R
python3 -B code/audit_naples_inclusion.py --before /path/to/pre-inclusion/checkpoint
```

The new R test evaluates the actual production guard and rejects altered coverage,
duplicate years, missing responses, invalid pre-trend screens and a changed selected
model. The ordinary pipeline requires only tracked harmonized inputs and registry
files. Re-export the registry explicitly after any authorised membership change.
Raw data and release markers were not changed. No push or circulation was authorised.
