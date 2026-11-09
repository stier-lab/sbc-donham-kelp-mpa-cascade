# Empty Landsat reference correction

Six documented empty second-reference slots were removed before averaging.
The source CSV is preserved. Observed zeros in nonempty slots remain.

Paired satellite MPA-years: 719 -> 718.
South La Jolla 1991 loses its artificial zero reference and is no longer a paired observation.
The 145 primary histories and MPA membership are unchanged. All diver inputs
and diver primary effect estimates are unchanged. Six Landsat estimates change.
Shared variance components cause small changes in other pooled biomass estimates.
All nine pooled signs and BH classifications persist.

Kelp's exponentiated pooled contrast is 2.24 (95% CI 1.43-3.49; adjusted P = 0.0020).
This is an exponentiated synthesis of design-specific contrasts, not a common physical abundance ratio.
The modern diver resistance/recovery results are unchanged.

The secondary matched before/after kelp check now meets its adjusted threshold: log change 0.693, 95% CI 0.168-1.218, adjusted P = 0.046, 17 MPAs (previous adjusted P = 0.054).
Its post-only contrast remains uncertain (adjusted P = 0.316).
This secondary change does not establish a predator-mediated cascade or erase placement confounding.

No-take rules and inherited programme-specific inclusion are documented separately.
Gull Island KFM remains primary; its omission sensitivity remains.
The missing habitat-specific classification list and exact extraction masks/coverage rules
are still unresolved. This correction does not reclassify other histories.

Rebuild satellite inputs explicitly through `code/R/export_landsat_annual_index.R`
with the source workspace path. After the primary refit, explicitly rebuild the
reference registry with `code/R/export_reference_registry.R` and the same source
workspace path, then rerun `Rscript code/R/run_all.R` to validate all dependencies.
The correction marker is validated before analyses and cached figure regeneration.
The before/after comparison and focused tests are reproducible from this checkpoint
with `python3 code/audit_landsat_reference_refit.py --before <checkpoint>`.
