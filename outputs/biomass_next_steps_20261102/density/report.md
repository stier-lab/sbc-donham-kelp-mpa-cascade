# Density baseline rebuild — diagnostic implementation

Parent session `9901638869734200a6633d0f88ec31ef`. All writes belong in the chosen scratch output directory; repository and raw source inputs are read only. No own continuity start, shared status edits, source changes, selection, model fits, messages, pushes or further agents.

## Implemented result

Six density histories use their own **558 assigned station-year uses**, preserving researcher assignments, all available original years and all **63 post-protection comparison-years**. Recovered **33 prior station-year uses / six comparison-years**: Abalone lobster 2011 (4); Point Vicente red and purple 2008 (6 each); Vicente lobster 2008–2010 (17). No additional station or year was recruited from raw availability.

All six recovered comparison-years have finite full-available-zone densities, but **only three have complete INNER/MID support**. ROCKY_POINT_S in 2008 has DEEP only, affecting all three Vicente density responses. Its station rows and all three 2008 annual rows remain present; shared-zone values are explicitly unknown, never an average of the remaining references. No zone spans the full Vicente density roster. Even the other three recovered years have different full available zone sets among stations: availability does not establish habitat or seasonal comparability.

| Comparison | Primary prior | Old all-month prior | Rebuilt full prior | Shared INNER/MID prior | Recovered station uses |
|---|---:|---:|---:|---:|---:|
| Abalone Cove SMCA / M. franciscanus | 1 | 1 | 1 | 1 | 0 |
| Abalone Cove SMCA / P. interruptus | 1 | 0 | 1 | 1 | 4 |
| Abalone Cove SMCA / S. purpuratus | 1 | 1 | 1 | 1 | 0 |
| Point Vicente SMCA / M. franciscanus | 4 | 3 | 4 | 3 | 6 |
| Point Vicente SMCA / P. interruptus | 4 | 1 | 4 | 3 | 17 |
| Point Vicente SMCA / S. purpuratus | 4 | 3 | 4 | 3 | 6 |

Primary and rebuilt full records contain 2011 at Abalone and 2008–2011 at Vicente. INNER/MID support reduces Vicente to 2009–2011. Thus all three Vicente responses meet three and four prior-year **count-only** rules in the full-zone version; they meet three but not four in the shared-zone version. Abalone meets neither. All six histories have changing station membership through time. These are response-specific diagnostics, **not** adopted fixed-station/shared-benthic-response eligibility, evidence of parallel trends, or permission to select pBACIPS. The adopted shared-response benthic limitation (including lobster biomass) remains unchanged. Prior means year < 2012; designation year 2012 is not prior. The 63 legacy recovery years use t > 0 to preserve that review's comparison exactly.

## Files and interpretation

- `raw_count_unit_registry.csv`: 4218 species-unit rows on 1406 distinct observed physical units. Stable SHA256 IDs are built from site/date/zone/method/transect; species IDs add classcode. IDs exclude comparison and species at the physical level, allowing common weights across taxa/uses. They are **dated sampling-unit identities**, not proof of permanent transects. Raw row locators are 1-based data records excluding the CSV header; source bytes are hashed. Retain the same IDs when sharing observations between density and biomass.
- `density_station_year_roster.csv`, `density_count_unit_uses.csv`, `per_visit_density.csv`: authorized membership and its count/visit joins; no implicit joins by generated row number.
- `station_zone_eligibility.csv`: exact available, missing shared and extra zones, explicit eligibility, original value and biomass-roster membership.
- `annual_station_density_variants.csv`, `annual_comparison_density_variants.csv`: original verified station values (inherited comparator, not claimed newly re-derived), rebuilt all-month/full available zones, and all-month/INNER+MID when both strata are supported. Full-role arithmetic means use every assigned station; any unavailable assigned value makes that variant unknown. Frozen offsets are retained. The original verified station comparator is distinct from the manuscript-facing primary values in `primary_annual_comparison.csv`.
- `recovered_station_year_uses.csv`, `recovered_comparison_years.csv`, `baseline_coverage_diagnostic.csv`: exact omissions, all three variants of the six recovered years (18 rows), primary before-year coverage and three/four-prior count checks. No outcome-based screening.
- `biomass_roster_count_units.csv`: **2214 current INNER/MID count-unit uses**, with original biomass series ID, role, exact site/year/zone/date/method/transect, count, area, unit density, visit density/total count/transect count, physical IDs, species and provenance. These reproduce all **1086** existing shared-zone visit uses. Other zones are available in the full registry, but are not part of this current biomass bootstrap export.
- Identity CSVs, `verification.json`, `input_hash_manifest.json`, `test_density_baseline.py`, `tests.json`, `tests.log`: numerical checks and reproducibility evidence. Tests are a separate command; a rebuild alone does not refresh the test log.

## Count and aggregation conventions

The established VRG/UCSB campus filter is preserved (UCSB DEEP/SPECIALO excluded; VRG DEEP retained). All retained sites here use VRG swath. Every focal species is joined to the physical-unit roster formed from **all observed taxa**, assigning a surveyed zero only when the focal classcode is absent from an observed unit. An explicitly missing count remains unknown. Multiple focal records within a unit sum as in the source reconstruction; exact duplicate source records are checked and cause failure. Current inputs have no exact duplicates or missing focal counts. No site/date/zone/unit with no survey is invented.

Counts divide by the established **60 m²** swath area. Average transects within site/date/zone/method; average all valid-month visits equally within zone/year; average zones equally within station/year; average assigned stations equally within each role/year. Repeat visits are neither pooled with unequal transect weights nor represented as independent annual samples. Method is retained in physical keys; the current saved visit schema omits it, so ambiguous multi-method joins fail explicitly. Date IDs are not inferred survey campaigns.

## Verification and limits

All 78 original shared-zone all-month ratios (63 post years plus baseline/designation years) reproduce with maximum absolute error 4.44e-16. All 1,086 existing biomass-roster visit densities reproduce to 1.78e-15. Independent tests reconstruct counts from the raw CSV using Python's standard CSV reader, then visits and annual ratios with dictionary grouping; they also check the 78 original full-zone all-month ratios. Synthetic edge tests cover missing strata, missing values, equal visit/zone weighting and ID collisions. Input hashes are checked before/after the rebuild.

This restores density evidence, not biomass size support, causal comparability or primary eligibility. Shared labels alone do not establish matching depths or survey seasons; station membership varies. No date-gap, size-donor or new exclusion rule is imposed. Raw/primary/paper files remain untouched. Parent integration and scientific policy decisions remain separate.

## Reproduction

Requires Python 3.9+ (tested here with the installed Python), pandas and numpy. The script imports no repository code and performs no work on import.

```sh
python rebuild_density_baseline.py --repo /path/to/sbc-donham-kelp-mpa-cascade --raw /path/to/MLPA_kelpforest_swath_2024.csv --output /path/to/scratch/density
python test_density_baseline.py --repo /path/to/sbc-donham-kelp-mpa-cascade --output /path/to/scratch/density > /path/to/scratch/density/tests.log
```

`DONHAM_ANALYSIS_REPO`, `PISCO_SWATH_CSV`, and `DENSITY_OUTPUT_DIR` are optional environment equivalents. If raw is omitted, the sibling data-repository convention is used. Output is mandatory; repository writes are restricted to outputs/biomass_next_steps_20261102. Raw input paths are never valid output locations. Portable input paths and full SHA256 fingerprints are recorded in the manifest.
