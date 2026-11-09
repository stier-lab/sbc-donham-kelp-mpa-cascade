# Programme and predator-size checks

Run `Rscript code/R/38_coauthor_programme_checks.R`, then
`Rscript code/tests/test_coauthor_programme_checks.R`. Module 38 is also sourced
in isolation by `code/R/run_all.R`, after the primary analyses.

These checks do not change the 145 retained primary contrasts, model selection,
zero offsets or inclusion rules. The retained-primary refit must reproduce
Table 2 within 1e-7. Each omission repeats the two response-specific joint REML
models and the nine-coefficient BH family. PISCO-only models omit the source
random intercept because there is only one source. They also change geography,
calendar coverage and precision, so they are not a controlled method experiment.

## What changed

- Removing KFM sheephead, Skunk Point or Gull Island KFM leaves all nine primary
  directions and adjusted detection classifications unchanged.
- PISCO-only estimates differ: lobster density and purple-urchin biomass are
  detected, while kelp is not. Do not claim that every programme independently
  reproduces the pooled pattern.
- KFM and PISCO sheephead density have opposite protected/reference mean log
  ratios at Anacapa and Harris Point over identical observed post-designation
  years. Both are positive at Scorpion and South Point. These are post-only
  descriptive levels, not the primary mixed change-and-level contrasts.
- All eight exploratory predator-biomass/urchin-biomass and predator/kelp
  regressions have nominal intervals spanning zero. Their eight-test BH family
  is separate from the primary tests. Two links also occur in Table 3.
- PISCO sheephead size checks retain surveyed zero-count transects, sum counts
  over non-canopy levels and average transects, zones and sites. All-size and
  >30 cm TL comparisons use identical paired years. Positive-count site-years
  with missing sizes are omitted, not converted to zero large fish. The check
  uses primary sheephead-density MPA membership, so biomass coverage can exceed
  the primary biomass cohort. It is not a pBACIPS refit.
- Lobster size checks combine PISCO density with count-weighted MPA-status-year
  size pools. The >=82.55 mm CL cutoff is the 3.25-inch legal-size reference,
  not a universal urchin-feeding threshold. Expected mass/fraction calculations
  do not propagate size-pool uncertainty or establish unbiased visual detection.
  Primary biomass is unchanged; a new harmonised release is required before
  substituting these reconstructed values into the primary analysis.

## Inputs and outputs

Standard checks consume vendored primary summaries and raw-response exports.
Size checks additionally require contributing-programme raw files under
`data/PISCO/` and `data/ALL_sizefreq_2024.csv`; these are not silently downloaded
or redistributed. Input MD5 tables record the exact local files. In a checkout
without raw fish inputs the size-status table reports that limitation; do not
mistake an older size output for a fresh calculation.

Output families are `table_s_programme_*`, `table_s_sheephead_programme_*`,
`table_s_pisco_island_sheephead`, `table_s_extended_pathway*`,
`table_s_sheephead_size_*`, `table_s_predator_size_*`, and `table_s_lobster_size_*`.
Each pathway has one point per MPA; sources are averaged using the Table 3
convention, including its independence assumption for response variance.
Source-specific primary contrasts are supplied before averaging. Exact annual
membership is exported for every size and shared-year check.

KFM raw processing exports timed-roving ratios separately as response `RD`.
The primary density cohort uses visual-fish transects (`Den`), not those timed
counts. No KFM fish sizes are reconstructed.

## Map correction

Long Point's inherited coordinate (-118.21, 33.24) is outside its supplied CDFW
reserve polygon. Figure 1 now uses a projected interior polygon point for its
display location; the original statistical metadata are not rewritten.
`table_s_map_coordinate_correction.csv` records the display coordinates and route.
Ocean-depth shading and its legend are omitted to make protection polygons
legible. Figure count remains five; no titles/subtitles were added.
