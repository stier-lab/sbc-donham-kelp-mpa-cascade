# Corrected Landsat masks and input verification

The 20 October 2026 corrected delivery resolves the reference-mask discrepancies recorded in `TOM_LANDSAT_SOURCE_AUDIT_2026-10-20.md`. The annual data and fitted results are unchanged.

## Checks completed

- Downloaded all nine files from Tom Bell's corrected delivery. Original files remain local under `data/external/tom_bell_2026-10-21/`; checksums are recorded in `data/harmonized/landsat_v2_verification.json`.
- Reviewed and ran the supplied Python reconstruction against `MPA_Final_2022.mat`. All nonempty series matched across 1984-2021; empty groups were reported separately.
- Independently reconstructed annual biomass from the cached `CAkelpCanopyEnv_2021.nc` using R and exact coordinate keys rounded to six decimal places. Duplicate or unmatched pixel keys fail the check. All 24 MPA and 34 nonempty reference series match the delivered workbook with zero absolute difference and identical missingness. The 14 empty groups contain only missing/zero placeholders in the workbook and are not biological zero observations.
- Matched all 60 existing input rows by site centroid, status and reference replicate. All match the delivered workbook. The six previously removed second-reference slots are still empty. The 718 paired MPA-years and fitted results therefore remain unchanged.

The earlier masks had reference-label swaps at South Point, Farnsworth, Naples and Campus Point; four reassigned pixels plus swapped labels at Harris Point; three omitted pixels at Skunk Point; twelve reassigned pixels at San Miguel; and absent second-reference masks at Point Conception, Abalone Cove and Point Vicente. The three second references are valid (220, 141 and 26 pixels respectively). Their annual columns were retained in the primary input.

## Reproduction and pipeline check

From the analysis root:

```sh
Rscript code/tools/verify_landsat_v2.R
Rscript code/tests/test_landsat_v2.R
Rscript code/tests/test_landsat_reference_correction.R
```

The full reconstruction requires the local delivery, cached NetCDF, previous delivery's site-centroid table and sibling data repository's empty-slot registry. `LANDSAT_V2_DIR` and `LANDSAT_NETCDF` override the first two paths. It writes the two compact comparison tables in `outputs/landsat_v2_verification/` and a fingerprint manifest only after the comparisons pass. It does not write primary data or fit models.

Module 40, called during harmonized-data loading, now requires that manifest and checks the harmonized inputs and comparison tables before fitting. Public checkouts can use this recorded verification without the locally held raw files. When local sources are present, their fingerprints must also match. A cached pass is distinct from a fresh raw-data reconstruction. Regression tests reject altered inputs, altered reports and a missing manifest.

## Scope

The supplier recovered some clustering splits by matching the delivered arrays. Exact reconstruction verifies the supplied masks' agreement with the measurements; it does not establish that the recovered split is a unique record of historical clustering. No new stochastic clustering was performed.

This verification leaves the existing equal-available-reference procedure intact. A single preferred-reference sensitivity remains unperformed. The fixed pre-2012 correlation criterion is not uniformly pre-designation for island MPAs. Historical footprint and inclusion decisions, including Judith Rock, remain separate contributor questions. The MATLAB script's approximate 1-3 km boundary selection is not an exact geodesic buffer. No claim is made to have rerun the original complete MATLAB extraction from a clean workspace.
