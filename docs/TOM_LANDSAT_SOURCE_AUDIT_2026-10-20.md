# Landsat source conventions and reconciliation

Tom Bell's 19 October reply and shared folder were inspected on 20 October 2026. The local originals are in `data/external/tom_bell_2026-10-20/` and excluded from redistribution. The audit does not replace the primary delivered input, choose a new reference, change site membership or run primary models.

Update 21 October: the corrected delivery resolves the mask discrepancies below and supplies the MATLAB workspace. All 24 MPA and 34 nonempty reference arrays now reconstruct exactly. See `LANDSAT_V2_VERIFICATION_2026-10-21.md`; the following records the earlier delivery, not an ongoing discrepancy.

## Verified conventions

The script reads the 1984–2021, 152-quarter release. Its MPA quarterly biomass is the sum of pixel wet biomass, with NaNs omitted. Quarters with more than 25% missing pixels are invalid; annual maxima are invalid when more than one quarter is missing. Observed zeros are retained. Single-reference sites have no meaningful second-reference array. All 60 rows in the currently used `MPA_Runs_new.csv`, including their missingness patterns, match the newly supplied MPAs/MPAC1/MPAC2 CSVs exactly at 1e-6 tolerance. This verifies delivery identity, not every footprint assignment.

Candidate references are selected outside MPA polygons using a nominal 1–3 km boundary band, then sometimes split into two groups. The implementation searches within 3 km of polygon vertices in latitude/longitude coordinates, uses `p_poly_dist` for the 1 km exclusion, and excludes pixels within the 24 selected polygons. It is not an exact geodesic polygon buffer. Cluster code 1 runs unseeded k-means on mean interpixel distances; codes 2 and 3 use latitude 34.0335 and longitude -118.5. Saved masks, not a new stochastic split, define the intended reconstruction.

Tom identifies the stronger pre-2012 correlate as his preferred reference. The script calculates quarterly and annual correlation diagnostics but does not replace both exported arrays with a selected single series. Our primary harmonization averages available reference arrays after removing six previously verified empty slots. These are different estimands; a preferred-reference sensitivity remains unperformed. A fixed 2012 cutoff includes post-designation years for Channel Islands MPAs, so it must not be described as a uniformly pre-designation selection rule.

## What reconstructs and what does not

`Rscript code/tools/audit_tom_delivery_2026_10_20.R` applies the documented missingness rules to the saved masks and cached 2021 NetCDF. All 24 MPA arrays reconstruct exactly, including missingness. Fifty of 72 same-labelled MPA/reference arrays match exactly. Eight additional nonempty reference masks match the opposite reference label at South Point, Farnsworth, Naples and Campus Point. Label exchange alone does not change an equal average of both arrays.

Five nonempty reference groups remain unmatched to either same-site reference array: both San Miguel groups, both Harris Point groups, and Skunk Point Ref_1. Three empty second-reference masks (Point Conception, Abalone Cove and Point Vicente) have nonzero supplied arrays; these are not assumed to be genuine zero observations or silently removed. Six empty reference groups at Crystal Cove, Laguna Beach and Dana Point are invalid in the delivery; they are outside the 20-MPA delivered primary-input labels. Empty-array MATLAB sums can return zero, so absence of a mask is not proof of zero canopy.

The 72 comparisons comprise 24 MPAs, 55 nonempty groups overall and 17 empty groups. The five unresolved nonempty groups and three nonzero empty slots require contributor/version confirmation. The script records every mismatch rather than identifying coincidental matching values or selecting a reconstruction that improves the scientific result. Outputs: `reconstruction_check.csv`, `mask_array_crossmatches.csv`, `current_input_matches.csv` and `pre2012_reference_correlations.csv`.

The current folder contains 26 accessible files. Its README also names files absent from the folder listing, including `MPA_Final_2022.mat`, shoreline workspace, `MPA_Runs_new.xlsx`, BACI/map ZIPs and the January 2023 patch workbook. A full fresh geometry/workspace rerun was therefore not claimed. The recovered mask workbook is byte-identical to the earlier email attachment.

## Scope and remaining decisions

The site list identifies San Miguel Special Closure separately from Judith Rock SMR; it does not justify excluding Judith Rock or assigning an SMR implementation date to the special closure. Tom found no dive-survey-polygon extraction. The 2026 Q2 workbook is a separate patch-level product and was not substituted for the 1984–2021 primary input. Current raw, harmonized and primary numerical outputs remain unchanged.

For KFM, Scotty confirms visual-fish data, juvenile/adult categories with a <10 cm juvenile definition, and no explicit YOY identification. Counts are separate from PISCO, but KFM invertebrate size-frequency data may be reused. The source processor demonstrably reads combined KFM/PISCO/LTER size records and supplies the same matched size pool to KFM and other programmes. Separate count collection therefore does not establish independent derived biomass estimates. Original age-category aggregation in the MBON handoff still requires checking; the reply is not permission to manufacture lengths or infer YOY from a juvenile category.
