# Comparison with original KFM individual measurements

November 2, 2026. This check covers every red- and purple-urchin site/year/species group in the combined size file, comparing exact millimeter frequencies and animal counts with the older individual-measurement CSV and Scott's original workbook. These are the KFM animal sizes used in the current biomass comparisons. It does not certify the other 24 taxa in the combined KFM file or count/density datasets.

| Independent source | Groups in combined file | Exact matches before corrections | Exact matches in verified pipeline input |
|---|---:|---:|---:|
| Older individual measurements, through 2011 | 904 | 873 | 903 |
| Scott's 2012–2022 workbook | 558 | 518 | 558 |
| 2023, without independent source | 20 | Unverified | Unverified |

A group means one site, year and species, not one animal. The corrected input matches all 1,461 groups marked verified. No new mismatch was found in that verified set. This confirms measurement identity, not total biomass uncertainty or ecological comparability.

## What disagreed in the combined file

- Historical records: 15 groups combined Wyckoff Ledge with excluded Miracle Mile; 14 substituted excluded Southeast Reef for SE Sea Lion; Hare Rock purple 1987 and Pelican Bay red 1995 differed from the older individual measurements.
- Recent records: 20 groups combined Wyckoff with Miracle Mile and 20 substituted Southeast Reef. All 40 corrected groups now match Scott's intended physical-site measurements exactly.
- Hare Rock now uses the 102 originals corroborated by the published NPS table. Pelican Bay remains unverified: 135 animals in the combined file versus 140 in the legacy source. The prior same-model sensitivity left the Scorpion red-urchin finding unchanged; see the linked report below.
- There are 20 unverified 2023 groups across the entire raw file, not just the four discussed previously. Existing site exclusions remove 16; the remaining four have no uses in the current review size pools or corresponding primary biomass histories. No 2023 groups were added.

## What exists only in the originals

The reverse comparison found 288 site/year/species groups without a separately labeled counterpart in the combined file: 188 historical and 100 from Scott. This is not equivalent to 288 lost analytical comparisons.

- Miracle Mile and Southeast Reef account for 69 groups. They remain excluded; the absence of separate labels partly reflects the known mixing/substitution problem.
- SE Sea Lion has 39 historical groups from 1985–2004 absent from the combined file. They are recorded as source-only candidates, not silently inserted into the existing analysis. Their use would require checking corresponding count coverage and existing eligibility rules.
- Rodes Reef, Fry's Harbor and Webster's Arch account for the remaining 180 groups. Their source availability is documented separately from the decision to include a comparison. Existing researcher exclusions are unchanged.

## How the check works

The script reads Scott's original workbook directly and independently confirms that all 55,142 rows retain site, species, exact size, count and start/end date in the normalized copy. It uses exact sizes, not Excel-corrupted size-bin labels. Legacy sizes are converted from centimeters to millimeters. Repeated identical measurements remain separate animals; no deduplication is performed.

Comparisons use physical site, year, species and the full size-frequency distribution. They do not claim one-to-one specimen identity where the combined file lacks specimen IDs or dates. Source-only groups are retained in a separate table. The original combined file is compared separately from the corrected input, preserving a visible record of every discrepancy.

Reproduce: `python3 code/tools/compare_kfm_original_measurements.py`. This run changed audit outputs only. Raw files, pipeline inputs, primary analysis and manuscript were not modified; all 487 protected hashes remain unchanged. No emails or push.

## Tables

- [All 1,482 groups](all_group_comparisons.csv)
- [Exact size-frequency disagreements](size_frequency_differences.csv)
- [Source-only groups](source_only_groups.csv)
- [Source-only summary](source_only_summary.csv)
- [Taxa checked and outside scope](taxon_scope.csv)
- [Input hashes and verification](verification.json)
- [Previously tested analysis effects](../kfm_remaining_sources_20261102/README.md)
