# No-take classifications and satellite source audit

Follow-up, 17 October: [source reconciliation](SOURCE_RECONCILIATION_2026-10-17.md)
verifies Naples against the published technical report, checks delivered
satellite missingness and reference-centroid distances, and identifies dated
KFM/PISCO overlap candidates. The original audit below is retained as history.

## Outcome

Several missing details are recovered, but this is not a clean close-out. The
source checks expose incorrect regulatory descriptions and a likely satellite
input problem that should be resolved before another coauthor circulation.
No primary input, inclusion rule, model, figure or manuscript file was changed
in this audit. Gull Island KFM remains primary with its omission sensitivity.

## No-take rules versus the inherited inclusion list

Jenn's manuscript comment 24, returned on 15 October, says monitoring groups
have a list of SMCAs treated as de facto rocky-reef SMRs and that she believes
she gave it to Emily. That list was not recovered from the emails or files
searched. Current CDFW rules do resolve several individual cases:

| Area | Verified current rule and implication |
| --- | --- |
| Blue Cavern Onshore SMCA | No-take, not sheephead-only. |
| Long Point, Point Dume and Santa Barbara Island SMRs | No-take reserves. Programme-specific omissions cannot be explained by permission to fish the other focal taxa. |
| Farnsworth Onshore SMCA | Selected pelagic fish and squid exceptions, not ordinary permission to harvest sheephead, lobster, urchins or kelp. It is not a ban on all finfish. |
| Point Dume SMCA | Selected pelagic fisheries are allowed; those exceptions do not establish ordinary focal-taxon harvest. Tribal exemptions require separate treatment. |
| Swami's SMCA | Shore hook-and-line fishing is allowed. The claim that all finfish take is prohibited is false. A functional rocky-reef classification would need an explicit gear/depth rationale. |
| Painted Cave and Anacapa Island SMCAs | Lobster take is allowed; they are not fully protected for the predator-urchin-kelp question. Sheephead is not their only protected focal taxon. |
| Cat Harbor and Dana Point SMCAs | Fishing for sheephead and harvesting lobster/urchins are allowed under the specified ordinary methods. They are not sheephead refuges. |
| Arrow Point to Lion Head Point SMCA | Recreational finfish/plants and broad commercial take are allowed under applicable regulations; not full focal-taxon protection. |

Sources, checked 16 October 2026: CDFW [Blue Cavern](https://wildlife.ca.gov/Conservation/Marine/MPAs/Blue-Cavern),
[Long Point](https://wildlife.ca.gov/Conservation/Marine/MPAs/Long-Point),
[Point Dume](https://wildlife.ca.gov/Conservation/Marine/MPAs/Point-Dume),
[Santa Barbara Island](https://wildlife.ca.gov/Conservation/Marine/MPAs/Santa-Barbara-Island),
[Farnsworth](https://wildlife.ca.gov/Conservation/Marine/MPAs/Farnsworth),
[Swami's](https://wildlife.ca.gov/Conservation/Marine/MPAs/Swamis-San-Elijo-Lagoon),
[Painted Cave](https://wildlife.ca.gov/Conservation/Marine/MPAs/Painted-Cave),
[Anacapa](https://wildlife.ca.gov/Conservation/Marine/MPAs/Anacapa-Island),
[Cat Harbor](https://wildlife.ca.gov/Conservation/Marine/MPAs/Cat-Harbor),
[Dana Point](https://wildlife.ca.gov/Conservation/Marine/MPAs/Dana-Point), and
[Arrow Point](https://wildlife.ca.gov/Conservation/Marine/MPAs/Arrow-Point-to-Lion-Head-Point).
`data/external/source_metadata_2026-10-16/no_take_review.csv` records each case.
These are current ordinary-fishing rules, not an audit of every historical
amendment, special permit or tribal exemption over the monitoring period.

The operative eleven-name `SHEEPHEAD_ONLY_MPAS` list is inherited analysis
curation, not a legal class. The archived Donham `pBACIPS_PISCO_V10.R` removes
these names before re-including selected sheephead records. The data-processing
repo repeats that operation in `code/R/07_combine_and_export.R`; its
`data/mpa_take_regulations.csv` still has VERIFY cells and inaccurate blanket
finfish language. Neither that file nor its code comments independently prove
contributor approval. No targeted regulatory confirmation was found in the
March/April emails checked.

The PISCO master-site CSV contains designation/status fields, not a recovered
habitat-specific de facto no-take list. Emily's archived code/notes and
`Old/Sites_Years.xlsx` did not supply that list either. The archived exclusion
code documents what was done, but does not justify each exclusion.

### What can be resolved without another email

Stop describing the eleven names as areas that protect only sheephead or all
permit take of another focal taxon. Stop saying the four re-included areas all
prohibit finfish fishing. Legal protection, programme/site suitability and
analysis inclusion must be separate columns. The individual current rules
above are verified; Jenn's exact monitoring classification and the site-specific
reason for excluding eligible non-sheephead histories remain open.

## Satellite files recovered

| Source | What it establishes |
| --- | --- |
| Tom, 2 July 2020, message `1731158abac79364`, `MPA_Sites_Controls.xlsx` | Original MPA/control centroids and pre-designation correlation-based control selection. Email describes controls generally 4-7 km from MPA centroids and at least 500 m beyond the edge, with exceptions. Annual analyses used maximum annual biomass. Workbook ends in 2018 and is not the current 1984-2021 table. |
| Tom, 19 December 2020, `1767de09ce2383a5`, `BACI.zip` | Separate quarterly/annual and biomass/proportion products were delivered. Not evidence that every later table used the same extraction. |
| Tom, 12 August 2022, `182945f103c3bafb` | Links the canopy/environment NetCDF; identifies a July 2022 Landsat update and full years 1984-2021. |
| Emily, 15 August 2022, `182a22aeadfdedf4` | Distinguishes the aggregate MPA values she was already using from the new pixel-level product; asks how those site means were calculated. |
| Tom, 21 October 2022, `183fa8c993c82b45`, `MPA_Ref_Locations.xlsx` | Actual MPA/reference pixel coordinates. Email says 1-3 km reference distances and removal of pixels inside any MPA, but the attached coordinate sets do not straightforwardly agree with that stated distance rule. |
| Tom, 2 April 2024, `18ea01364b19f9c7` | MPA/control canopy and environmental products were still described as intended releases accompanying this paper; not an exact archived public release identifier. |

The linked Drive file is `1353POLzjcAZ5lrw9oHd8D9hLpSLe_AlW`,
`CAkelpCanopyEnv_2021.nc`, 207,796,876 bytes. Retrieved metadata report a
16 February 2023 modification; this is not proof of byte identity to its July
2022 version. Original attachments/NetCDF are cached in ignored
`local_archive/metadata_audit_2026-10-16/`; metadata and comparisons are retained
under `data/external/source_metadata_2026-10-16/`.

### Resolved: source quantity and most spatial identifiers

The actual NetCDF has 332,640 stations, 152 quarters and years 1984-2021.
Its case-sensitive `Units` attributes distinguish:

- `biomass`: `wet Kg/900m^2 pixel`, described as giant-kelp canopy biomass.
- `area`: `m^2 canopy/900m^2 pixel`, a different variable.
- `_FillValue = -1`, distinguished from observed zero.

All pixel coordinates in the workbook's 72 MPA/reference slots match the
NetCDF at six decimal degrees. Twenty workbook MPA centroids map to all twenty
names in the delivered CSV. The CSV repeats each MPA's same centroid in all
three status/replicate rows; those lat/lon columns are not reference locations.

Summing pixel biomass within each quarter and taking annual maxima reproduces
38 of the 60 CSV rows exactly on all jointly finite years (at least twenty
years required). Only six also reproduce the complete missingness pattern.
Other differences include missing-value treatment, nonzero value differences,
and reference labels swapped between the workbook and CSV. Thus the source
quantity is strongly supported, but a full byte/numeric reconstruction is not
complete. Keep the analysis derivative dimensionless and do not equate wet
canopy biomass with diver-derived dry foliar biomass.

### New input concern: empty reference slots encoded as zero

Six CSV second-reference rows contain zero for all 38 years and correspond to
empty Ref_2 slots in the recovered workbook:

- Skunk Point, Santa Barbara Island, South La Jolla and Cabrillo SMRs.
- Point Dume SMCA and Point Dume SMR.

All six contribute Landsat histories to the current primary synthesis.
The inherited and current preprocessing average reference replicates, so
including an empty all-zero slot halves every positive reference mean when
the actual reference is observed. In ten MPA-years, it also supplies a zero
where removing the candidate empty slot leaves no observed reference at all.
This is strong evidence of a structural-zero problem, not proof that the
reference kelp was absent. The 2022 workbook has three further empty Ref_2 slots
whose delivered rows are not all zero; those cannot be removed by the same rule.

The descriptive arithmetic check is
`satellite_empty_slot_arithmetic.csv`; no primary inputs were changed. Do not
claim that fixing this doubles a fitted MPA effect. Before/after differencing
can cancel constant scaling, whereas zero offsets, missing years, routing and
post-only estimates can change the fitted result. Refits are needed.

### Resolved: some satellite references overlap

Naples Ref_2 and Campus Point Ref_2 share 320 pixels: 80.8% of the former's
396 pixels and 63.9% of the latter's 501 pixels. Their finite values match the
CSV's reference replicate 1 (labels were swapped). Four pairwise overlaps
occur across the recovered reference sets overall; the others involve San
Miguel/Judith Rock/Harris Point. The Naples/Campus overlap directly affects
retained histories. It is no longer correct to call all satellite reference
overlap unverified. Exact cross-history covariance is not estimated here.

### Still unresolved: extraction rule and masks

The October 2022 email describes reference pixels 1-3 km from an MPA centroid.
Measured distances from the supplied centroids exceed 3.1 km in 30 of 31
nonempty reference sets; Carrington's Ref_1 spans about 4.50-6.13 km.
Rounding a centroid to four decimals cannot account for that discrepancy.
The email and workbook must therefore not be treated as a single verified
geometry recipe. The workbook identifies actual candidate pixels, not their
selection code or the exact rule used for every delivered series.

Remaining details: exact extraction script/version, valid-pixel thresholds,
landmask/reprocessing version, divergent reference slots, and the explanation
for the residual annual-value/missingness differences. The January 2023 and
September 2026 patch/metapopulation emails concern another workflow; their
patch thresholds and pixel shifts must not be transferred to this MPA table.

## Follow-Up Correction

Following Adrian's authorization, the six verified empty second-reference rows
were removed before averaging. Observed zeros were retained; missing reference
observations were not replaced with zeros. The versioned handoff contains 718
paired MPA-years, and diver inputs are unchanged. See
`docs/LANDSAT_REFERENCE_REFIT_2026-10-16.md` and
`data/harmonized/landsat_reference_correction.json` for the correction receipt.
The notes below describe the state when the source audit was first completed;
the corrected manuscript export is V10, with V9 preserved.

## Original Next Corrections

1. Reconcile and remove the six confirmed candidate empty reference slots in
   an explicit, versioned input correction; refit the affected Landsat models
   and downstream kelp/historical summaries before changing result text.
2. Replace incorrect legal descriptions and build a taxon/gear/site-specific
   inclusion ledger. Restore or exclude histories only after checking coverage
   and contributor rationale, not from the old variable name.
3. Use recovered pixel identities to update the shared-reference sensitivity;
   preserve the distinction between matched pixels and unresolved geometry.
4. Ask Jenn/Emily for the original de facto list and site-specific inclusion
   reasons. Ask Tom which geometry/mask/coverage rule produced this exact CSV,
   with the mismatches and empty slots attached as a short diagnostic table.

No email sent. The manuscript and V9 exports remain unchanged and should not
be described as cleared of these newly identified input/classification issues.
Reproduction commands and source receipts are in the evidence README. Four
focused source-evidence regression tests pass; no new model fits were run.
