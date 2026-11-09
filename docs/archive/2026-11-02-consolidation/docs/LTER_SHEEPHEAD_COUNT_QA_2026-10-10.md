# LTER Sheephead COUNT QA: 2026-10-10

**CURRENT STATUS (2026-10-10): COUNT CORRECTED; FULL REFIT COMPLETE.**

Exact-release COUNT semantics were verified, the focused correction was applied to both canonical export copies, and preservation tests passed.
The parent reports the full refit completed in 3.8 minutes with exit 0; primary membership remains 144. Table 1 was refreshed (LTER sheephead biomass mean 9.734, SE 2.087).
Numerical refit outcomes are documented separately in the forthcoming [LTER COUNT refit report](/Users/adrianstier/Donham-Stier-CA-MPA-2026/docs/LTER_COUNT_REFIT_2026-10-10.md), being written by Anscombe. This QA note records input provenance and correction tests, not the numerical refit audit.

## Historical Read-Only Audit (2026-10-10; Before Correction)

The following original audit sections describe the pre-correction snapshot. Their references to an unresolved issue or unchanged data/estimates are historical, superseded by the current-status banner and authorized correction follow-up below.

### Scope And Disposition

Read-only audit of the existing LTER sheephead size-to-biomass calculation and its exported histories. This note preserves completed checks; it does not authorize a harmonization change. No data-processing scripts, input data, primary estimates or manuscript files were edited, and no effect-size or pooled models were refitted. The parent is disclosing the unresolved issue during a readability revision, not changing data or estimates.

**Confirmed:** the current calculation sums one allometric individual-mass estimate per sized observation without weighting by COUNT. Some observations have COUNT greater than one, and the existing harmonized biomass and annual log contrasts reproduce that unweighted calculation.

**Unresolved:** data-owner confirmation of the exact local-export COUNT convention, whether the present treatment was intentional, and the consequences of any approved correction for model selection, primary effect sizes, uncertainty and pooled results. This audit makes no claim that pooled results would be unchanged.

## Inputs And Source Locations

- Raw fish: [Annual_fish_comb_20240307.csv](/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data/LTER/Annual_fish_comb_20240307.csv).
- Site mapping: [SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv](/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data/MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv).
- Upstream import, filtering, biomass and aggregation: [06_lter_processing.R](/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/code/R/06_lter_processing.R:753).
- LTER allometric function: [01_utils.R](/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/code/R/01_utils.R:179); coefficients: [00c_analysis_constants.R](/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/code/R/00c_analysis_constants.R:246).
- Normalization and log contrast helpers: [01_utils.R](/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/code/R/01_utils.R:383).
- Upstream exports: `/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/output/harmonized/harmonized_raw_responses.csv` and `/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/output/harmonized/harmonized_response_ratios.csv`.
- Canonical analysis copies: [harmonized_raw_responses.csv](/Users/adrianstier/Donham-Stier-CA-MPA-2026/data/harmonized/harmonized_raw_responses.csv) and [harmonized_response_ratios.csv](/Users/adrianstier/Donham-Stier-CA-MPA-2026/data/harmonized/harmonized_response_ratios.csv).
- Final primary membership: [sumstats_final.csv](/Users/adrianstier/Donham-Stier-CA-MPA-2026/data/sumstats_final.csv:29).
- Retained history index: [table_s_pair_history_atlas_index.csv](/Users/adrianstier/Donham-Stier-CA-MPA-2026/tables/table_s_pair_history_atlas_index.csv:29).
- Explicit Campus Point source selection: [08_effect_sizes.R](/Users/adrianstier/Donham-Stier-CA-MPA-2026/code/R/08_effect_sizes.R:1054); linear-estimate input handling: [add_linear_effect_size](/Users/adrianstier/Donham-Stier-CA-MPA-2026/code/R/08_effect_sizes.R:864).

CSV line locators below count the header as line 1.

### Input Fingerprints At Audit

| Input | MD5 |
|---|---|
| Raw fish CSV | `ddf0cb97132a35aede43a8057709ffee` |
| Site mapping CSV | `3b664e8e393e72e59ffe45e854cce9ae` |
| PISCO species attributes CSV | `6ccfacd18fff07208f510384e914a3b2` |
| Upstream 06_lter_processing.R | `5c2ecf1a9f85066a6d2067c25b4c62e7` |
| Upstream 01_utils.R | `896cc42a328364ec7206e98d965f1760` |
| Upstream 00c_analysis_constants.R | `5a435683ffd751841c054386f4b7a455` |
| Canonical harmonized raw responses | `91bf103c6b64bd4beba305c3ad80e7e0` |
| Canonical harmonized response ratios | `93406d70cca3af23b3ba1963b2599dbc` |
| Final primary effect-size CSV | `0dc1df3c099234762484c99b35a1037d` |
| Pair-history atlas index | `0038ddec697e4bc9b94fd95d3bf192f3` |
| Analysis 08_effect_sizes.R | `3fb51114cc53f8df013a629d5a58a2c2` |

These identify the inspected snapshot, not a requirement that future inputs remain identical.

## Raw Row Meaning And Counts

The raw file contains **93,811 data rows**, including **1,499 Semicossyphus pulcher rows**. COUNT reads as integer and SIZE as numeric. COUNT has no NAs or fractional values. All sheephead observations have AREA = 80 m2. Counts greater than one occur with a single recorded size, so these are not uniformly one-fish rows.

| COUNT | All raw sheephead rows | Retained mapped rows |
|---:|---:|---:|
| -99999 | 10 | 0 |
| 0 | 627 | 350 |
| 1 | 709 | 456 |
| 2 | 63 | 44 |
| 3 | 47 | 36 |
| 4 | 20 | 12 |
| 5 | 9 | 7 |
| 6 | 9 | 3 |
| 7 | 1 | 1 |
| 8 | 2 | 1 |
| 10 | 2 | 2 |
| Total | 1,499 | 912 |

The site join preserves the 1,499 raw sheephead rows: no join fan-out occurs for this subset. Applying the current nonempty-MPA mapping, AHND exclusion, Campus Point/Naples restriction and COUNT != -99999 filter retains 912 rows. Of these, 562 have positive counts and all have positive sizes; there are no positive-count missing/nonpositive sizes. The 350 zero-count rows have SIZE = -99999; their invalid allometric mass becomes NaN and is zero-filled by the current code.

| Retained MPA/status | Rows | Positive-count rows | COUNT > 1 rows | Sum of COUNT |
|---|---:|---:|---:|---:|
| Campus Point / MPA (IVEE) | 110 | 10 | 0 | 10 |
| Campus Point / reference | 142 | 28 | 2 | 30 |
| Naples / MPA (NAPL) | 492 | 461 | 93 | 659 |
| Naples / reference | 168 | 63 | 11 | 89 |

The **106 retained multi-fish rows** are NAPL 93, AQUE 11 and MOHK 2. The earlier count of 88 multi-fish NAPL/IVEE rows excluded year 2000; 93 includes five NAPL rows in 2000. These are different scopes, not conflicting counts.

### Upstream Semantics Caveat

The [official annual reef-fish protocol](https://sbclter.msi.ucsb.edu/external/Reef/Protocols/Kelp_Forest_Community_Dynamics/SBC_LTER_protocol_Reed_Kelp_forest_community_Density_fish_20181128.pdf) describes recording fish numbers, species and total lengths estimated to the nearest cm within 40 m x 2 m transects. This and the local COUNT/SIZE patterns support interpreting COUNT as fish multiplicity associated with a sized record. However, the precise column definition and data-entry convention for this local combined export were not independently recovered: the EDI metadata endpoint for the package identified in upstream DATA_SOURCES.md (`knb-lter-sbc.17.29`) was inaccessible, including an HTTP 403 response. Confirm the local export's convention with Emily/Miller before changing harmonization. Do not assume that repeated same-size rows are accidental duplicates or silently deduplicate them.

## Calculation And Aggregation Trace

The sheephead section reads the raw CSV directly, filters species, joins Sites2, restricts mapped MPAs/sites, and removes COUNT = -99999. There is **no prior aggregation or expansion into individual fish** before line 803.

At line 803 the implemented row contribution is `W(SIZE) = 0.0144 * SIZE^3.04`, not `COUNT * W(SIZE)`. NAs/NaNs are then set to zero. Subsequent operations are:

1. Sum COUNT and biomass separately within YEAR, SITE, TRANSECT, status, CA_MPA_Name_Short, ChannelIsland and MPA_Start.
2. Average transect totals within YEAR, SITE, status, CA_MPA_Name_Short, ChannelIsland and MPA_Start.
3. Average site means within status, CA_MPA_Name_Short, MPA_Start and YEAR.
4. Divide biomass by 80 for the raw response export.
5. For annual log contrasts, scale biomass by the common maximum within MPA/taxon across statuses and years, add half the minimum positive proportion, and calculate log(MPA/reference).

None of these operations later multiplies biomass by COUNT. The in-memory alternative below weights each sized row by its count while keeping every other grouping, averaging and normalization rule unchanged. It is a diagnostic comparison, not a replacement dataset or a model refit.

## Confirmed Numerical Examples

For NAPL, 2000-10-04, transect 8:

| Raw CSV line | SIZE (cm) | COUNT | Current contribution (g) | Count-weighted contribution (g) |
|---:|---:|---:|---:|---:|
| 53974 | 15 | 1 | 54.16016 | 54.16016 |
| 53975 | 20 | 1 | 129.86547 | 129.86547 |
| 53976 | 25 | 2 | 255.91758 | 511.83516 |
| 53977 | 30 | 4 | 445.46246 | 1781.84982 |
| 53978 | 35 | 1 | 711.75307 | 711.75307 |

These five rows represent COUNT = 9. Current transect biomass is 1,597.159 g, versus 3,189.464 g with count weighting. This 2000 Naples example is **not** an annual paired input: the reference observation is absent for that year. It is not part of the retained Campus Point history.

The implemented, unweighted reconstruction matches all **90 exported LTER sheephead biomass values** to a maximum absolute difference of `4.796163e-14`, and all **45 annual log contrasts** to `4.662937e-15`. Both the upstream exports and their canonical analysis copies match. There are 92 unpaired status/year aggregates before requiring both statuses; Naples MPA 2000 and Campus Point reference 2004 have no counterpart and do not enter these paired exports.

| Series/year | Current raw biomass (g/m2) | Count-weighted diagnostic (g/m2) |
|---|---:|---:|
| Campus Point reference, 2020 | 2.95515610 | 3.75887078 |
| Naples MPA, 2019 | 75.40261601 | 86.28064147 |
| Naples reference, 2019 | 23.32125539 | 45.23874480 |

Campus Point's 2020 annual log contrast changes from -2.4382246 to -2.6767381 in the diagnostic calculation (difference -0.2385135). Naples' 2019 annual log contrast changes from 1.1732296 to 0.6455636 (difference -0.5276660). The factors vary across status and year, so this discrepancy does not cancel as a constant common scaling.

Across paired years, count weighting changes 1/22 Campus Point annual log contrasts and 18/23 Naples contrasts. Among exported raw biomass values it changes 0/22 Campus Point MPA values, 1/22 Campus Point reference values, 18/23 Naples MPA values and 4/23 Naples reference values. These are descriptive input comparisons, not revised primary results.

### Retained Campus Point History

The final primary key `Campus Point SMCA | S. pulcher | LTER | Bio` is retained, with Model = Linear, AnalysisType = pBACIPS and N = 22. Its supplied annual history is **2001-2003 and 2005-2023**, missing 2004; designation is 2012. The atlas index records 22 pairs, 11 before/designation and 11 after.

- **2019 enters this retained history:** BA = After, time = 7; its Campus Point log contrast is 3.1582132 and does not change under the count-weighted diagnostic.
- **2020 enters this retained history:** BA = After, time = 8; this is the affected Campus Point observation. Its current raw reference value is at canonical raw-response CSV line 445 and its annual log contrast at response-ratio CSV line 2692.
- The raw NAPL 2000 and paired Naples 2019 examples are not Campus Point observations. The Naples/LTER/sheephead/Bio key is not in the current final primary set.

The explicit Campus Point branch supplies this full history to add_linear_effect_size; the function fits the supplied data without a post-only year filter. **No corrected Campus Point estimate, SE, selected-model comparison or pooled result was computed.** Consequences require a separately authorized data-owner review and downstream refit.

## Programme-Specific Sheephead Allometries

| Programme | Equation implemented for individual mass | Count handling | Source |
|---|---|---|---|
| PISCO, UCSB and VRG | `W = 0.0238 * TL^2.8847` | Row biomass explicitly multiplied by count | Species attributes, CSV line 150; 04_pisco_processing.R TL branch |
| LTER | `W = 0.0144 * TL^3.04` | Currently one mass estimate per sized row, without COUNT weighting | Constants, helper and 06_lter_processing.R line 803 |
| KFM | No sheephead biomass conversion | Count-derived Den/RD responses | 05_kfm_processing.R sheephead section |

The PISCO lookup specifies TL in cm, mass in g, source_LW_LC = VRG and reference notes = Williams et al. (2013). [PISCO species attributes](/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data/PISCO/spp_attribute_table_downloaded_9-13-22_SHSPUL.csv:150) are joined to both UCSB and VRG fish rows; [04_pisco_processing.R](/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/code/R/04_pisco_processing.R:778) explicitly multiplies the TL conversion by count. Original-source verification of the PISCO coefficient attribution was not undertaken here.

LTER constants and the helper label TL in cm and mass in g. KFM's [fish section](/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/code/R/05_kfm_processing.R:839) supplies count-derived responses, not sheephead biomass. Consequently, one universal sheephead conversion is not an accurate description of this analysis.

### Hamilton And Caselle (2015): Confirmed Support

The local PDF identity was verified from its title page: Hamilton, S.L. and Caselle, J.E. (2015), *Exploitation and recovery of a sea urchin predator has implications for the resilience of southern California kelp forests*, Proceedings of the Royal Society B 282:20141817, DOI `10.1098/rspb.2014.1817`.

Its Methods section 2(b), PDF page 2, explicitly reports `W = 0.0144 * TL^3.04`, TL in centimetres, n = 499 and r2 = 0.98. The continuation on page 3 says the estimated weights and density information were used to calculate biomass in metric tons per hectare. Thus it supports **both the body-size plus density biomass workflow and the exact coefficients used by the current LTER helper**. It does not establish that current PISCO uses those coefficients, nor does it validate omitting COUNT from the local LTER export. This resolves the earlier uncertainty about Hamilton (2015) supporting the LTER coefficients; the upstream code's Cowen attribution is not treated as proof of the coefficient's original provenance.

- Local PDF: [Hamilton_2015_ExploitationRecoverySeaUrchinPredator.pdf](/Users/adrianstier/Donham-Stier-CA-MPA-MS-2026/literature/pdfs/Hamilton_2015_ExploitationRecoverySeaUrchinPredator.pdf).
- Existing extracted text: [Hamilton_2015_ExploitationRecoverySeaUrchinPredator.txt](/Users/adrianstier/Donham-Stier-CA-MPA-MS-2026/literature/source_text/Hamilton_2015_ExploitationRecoverySeaUrchinPredator.txt:118), with the page-3 continuation around line 303.
- PDF SHA256: `05a978625c1160c9a7200416f46c867fb65b1faf591fcce31ea14cdb9cf7dd20`.
- Existing text SHA256: `70b6daeb0214447902543512dcca16550b73591ea85158604e7a910ea1560557`.

## Completed Read-Only Verification Approach

The numerical checks were run with `Rscript --vanilla -e` in a local R scope using read.csv, dplyr and tidyr. No complete upstream module was sourced. No output file was written. To reproduce the completed comparison:

1. Read the raw fish and site CSVs listed above, retaining the original CSV line as row index plus one. Select SCIENTIFIC_NAME == Semicossyphus pulcher, left-join SITE to site_id, and apply exactly the mapped-MPA, AHND and COUNT filters described above. Check that the sheephead join has 1,499 rows and the retained subset has 912.
2. Parse, but do not source, upstream 00c_analysis_constants.R and 01_utils.R. Evaluate only the two SHEEPHEAD_BIOMASS_ALLOMETRIC constants and the definitions of bio_sheephead, calculate_proportions, calculate_log_response_ratio and safe_divide in a new local environment. This uses the current source equations and normalization helpers without executing pipeline work.
3. Compute the implemented biomass with bio_sheephead(SIZE), replacing its NAs/NaNs with zero exactly as the source does. In a second in-memory column only, calculate COUNT times that biomass. Do not overwrite any input column/file used by the production pipeline.
4. Carry both columns through the exact transect sum, site mean and MPA/status/year mean grouping keys listed in the aggregation trace. Divide both by 80 for comparison to raw biomass exports.
5. For each column separately, call calculate_proportions with y = SPUL, pivot status to mpa/reference, and call calculate_log_response_ratio. Keep missing pairs absent. Do not fit an lm, pBACIPS or meta-analysis model.
6. Compare the current reconstruction to both copies of the harmonized raw and ratio CSVs, filtering source = LTER, species = Semicossyphus pulcher and resp = Bio. Join raw values on MPA/year/status and ratios on MPA/year. A tolerance of 1e-10 covers the observed CSV rounding; expect 90 matching raw values and 45 matching annual log contrasts.
7. Compare the count-weighted diagnostic to the reconstructed current values, and tabulate changed values by MPA/status/year. Recover the row/transect and annual examples above. This comparison is not an assertion about fitted outcomes.
8. Read sumstats_final.csv and the atlas index to verify the retained Campus Point key and its 22-year history; confirm 2020 is present. Preserve the input fingerprints and disclose any differences from this audit snapshot rather than silently reusing its numerical findings.

The following existing read-only command was used to verify the Hamilton Methods passage directly from the PDF without creating another text file:

```bash
pdftotext -f 2 -l 3 -layout \
  /Users/adrianstier/Donham-Stier-CA-MPA-MS-2026/literature/pdfs/Hamilton_2015_ExploitationRecoverySeaUrchinPredator.pdf - \
  | rg -n -C 4 '0\.0144|length.*weight|values and|density information|estimate biomass'
```

## QA Record Disposition

LTER sheephead biomass currently sums one allometric individual-mass estimate per sized observation without weighting by COUNT. This reproduces existing harmonized exports and differs from count-weighted totals for 106 retained multi-fish rows. The affected Campus Point 2020 observation enters the retained primary history. Local-export row semantics require data-owner confirmation; implications for the primary Campus Point estimate and pooled results remain unresolved pending separately authorized harmonization review and refitting. The readability revision discloses the issue and the programme-specific equations without changing data or estimates.

# Authorized Correction Follow-Up (2026-10-10)

This follow-up supersedes the earlier unresolved COUNT-semantics status below.
The original audit and its before-correction fingerprints are retained as a
historical record. User authorization did not itself establish row semantics.

## Exact Primary Metadata Verification

The EML for **knb-lter-sbc.17.39**, published 2024-03-07, was retrieved via the
DataONE resolve service after direct PASTA access returned an authorization
error:

- [Primary EML via DataONE](https://cn.dataone.org/cn/v2/resolve/https%3A%2F%2Fpasta.lternet.edu%2Fpackage%2Fmetadata%2Feml%2Fknb-lter-sbc%2F17%2F39).
- [Canonical EML identifier](https://pasta.lternet.edu/package/metadata/eml/knb-lter-sbc/17/39).
- Exact-release DOI: [10.6073/pasta/6ea061f660a7dc939686f51613b5ff0e](https://doi.org/10.6073/pasta/6ea061f660a7dc939686f51613b5ff0e), explicitly present as the EML DOI alternate identifier. The [EDI citation service](https://cite.edirepository.org/cite/knb-lter-sbc.17.39?style=RAW) independently returns version 39, Reed and Miller, and publication date 2024-03-07. Verified 2026-10-10; this matches the manuscript's Reed and Miller (2024b) DOI, so no COUNT citation-route replacement is needed.
- `dataset/dataTable/physical/objectName`: `Annual_fish_comb_20240307.csv`.
- `physical/authentication[@method='MD5']`: `ddf0cb97132a35aede43a8057709ffee`.
- `dataTable/numberOfRecords`: **93,811**, matching the local source file.
- `attribute[attributeName='COUNT']/attributeDefinition`: "The number of individuals of the same size counted".
- SIZE is estimated individual size, in centimetres; AREA is surveyed area, in square metres. COUNT and SIZE use `-99999` for unavailable/not-recorded values.

The local raw file's MD5 matches the exact-release EML authentication value.
Thus COUNT is a row multiplicity, not one unit of biomass per CSV row. This is
documentary verification, not a claim of a new communication with the provider.
The protocol linked by the EML describes fish total length and 40 x 2 m surveys:
[official reef-fish protocol](https://sbclter.msi.ucsb.edu/external/Reef/Protocols/Kelp_Forest_Community_Dynamics/KFCS_Protocol_Density_of_reef_fish.pdf).
The currently served protocol is dated 2025-10-15; it corroborates the workflow
but does not replace the exact 2024-release EML. The older source inventory's
`.17.29` label does not identify this checksum-matched local CSV.

## Approved Route And Applied Change

The parent preserved the analysis inputs, code and outputs under
`local_archive/lter_count_refit_2026-10-10/before/` before refitting. Original
upstream scripts are also preserved there in `upstream_raw_scripts.tar` (MD5
`f69dfe1328dc59ecaa88fac8d79c1731`). The upstream module-06 pre-edit diff was clean.
The original upstream summary is preserved as `before/upstream_average_responses.csv`.

The approved canonical route remains the local sibling data-processing repo,
then its copied tracked analysis inputs. The public archived remote was not
unarchived or pushed. The focused upstream change in
`code/R/06_lter_processing.R` is:

```r
lter.fish.sub.site$biomass <- lter.fish.sub.site$COUNT * bio_sheephead(lter.fish.sub.site$SIZE)
```

No allometry, missing-size policy, exclusion, site mapping, hierarchy, density
calculation or non-LTER processing changed. The new sibling-repo driver,
`code/R/export_lter_sheephead_count.R`, rebuilds old and corrected sheephead
histories with the existing helpers and grouping, proves old-history
concordance, then replaces only LTER / Semicossyphus pulcher / Bio records.
Default execution is read-only; production writes require `--apply`.

Only `harmonized_raw_responses.csv` and `harmonized_response_ratios.csv` were
updated, in both the sibling `output/harmonized/` and analysis `data/harmonized/`:

| Export | Target rows | Changed rows | Non-target rows preserved | Corrected MD5 |
|---|---:|---:|---:|---|
| Raw responses | 90 | 23 | 6,190 | `206e530c49963ef870d981e14fee24fb` |
| Response ratios | 45 | 45 | 3,296 | `ed8e5f3413e1ffba8eea448743c9068d` |

All 45 ratio records change in at least one normalized field; **19** annual
log contrasts change beyond `1e-10` (Campus Point: 2020 only; Naples: 18 years).
Every non-target CSV line is byte-preserved, as are target records' other
parsed columns. Landsat and site-metadata inputs are unchanged. The two export
locations have identical corrected checksums. The upstream checksum manifest
was refreshed; the parent owns the analysis manifest after registry refresh.

The upstream `tables/average_responses.csv` was independently recomputed using
module 07's `Rmisc::summarySE` specification. All 28 group keys and N values
are unchanged; the 27 non-target groups differ by at most `1.136868e-13` from
recomputation, and their original CSV lines were retained. Only sheephead /
LTER / Bio was replaced: N=90, mean=9.733950158468197,
SD=19.799286460837898, SE=2.087028042079398, CI half-width=4.146880264884062.
The parent separately added everyday analysis-side summary regeneration from
tracked raw responses; this correction did not edit that helper.

## Verification And Reproduction

From the sibling data-processing repo, using the preserved analysis baseline:

```sh
Rscript code/tests/test_lter_sheephead_count.R "$BASELINE_HARMONIZED_DIR"
Rscript code/R/export_lter_sheephead_count.R "$BASELINE_HARMONIZED_DIR" output/harmonized "$ANALYSIS_HARMONIZED_DIR"
# After reviewing the dry run and preserving all baselines:
Rscript code/R/export_lter_sheephead_count.R "$BASELINE_HARMONIZED_DIR" output/harmonized "$ANALYSIS_HARMONIZED_DIR" --apply
```

The tests write temporary copies only. Passed checks include exact raw-release
fingerprint, 912 retained rows / 106 COUNT>1 rows, old-history concordance
(maximum discrepancy: raw `4.796163e-14`, ratio fields `4.547474e-12`), actual
module-06 count assignment, 19 changed annual contrasts, unchanged Campus MPA
raw values, both-target parity, non-target line/column preservation, immutable
baseline, missing-key refusal, baseline-as-target refusal and repeat/double-
correction refusal. A separate post-export read-only check confirmed all
non-target rows and lines, both production-copy checksums and untouched Landsat
and site metadata.

The descriptive summary can be reproduced without processing raw surveys or
fitting models:

```r
raw <- read.csv("output/harmonized/harmonized_raw_responses.csv")
fresh <- Rmisc::summarySE(raw, measurevar = "value",
                        groupvars = c("taxon_name", "source", "resp"))
```

The retained primary branch is Campus Point / S. pulcher / LTER / Bio in
analysis module 08: its full 22-year linear history includes the affected 2020
observation. Naples annual Bio also feeds temporal and sensitivity analyses,
although Naples LTER sheephead Bio is not a retained primary key in the
pre-correction final set. Downstream models, aggregates, tables, cache and
figures therefore require the parent's refit; no unchanged pooled result is
claimed. The parent's run must regenerate the reference registry after module
08 because module 27 checks current harmonized-input and primary-result hashes.
Dan's atlas (29) and six-panel primary pathway display (31) are already in the
main orchestrator and must be regenerated from the refitted inputs/results.
No model fit, full harmonization/pipeline run, manuscript edit, email, commit or
push was performed by this focused correction task.
