# KFM donor-alias audit — November 2, 2026

## Outcome and scope

**No additional post-2011 year can be safely reconstructed from the local evidence examined.** Keep the six biomass histories provisional and their density counterparts retained. No new rule or correction was adopted. This was a read-only child audit under parent session `393ab715dc104db49d989bacdd9f2f3e`; all new files are in this scratch directory. No repository, status, primary output, raw data, or consolidated table was changed. No messages, downloads, models, pushes, or additional agents were launched.

The six histories are four KFM comparisons and two PISCO comparisons using KFM size donors—not six KFM-program comparisons. Independent checks reproduce 65/65 historical alias histograms (51 Wyckoff/composite and 14 Southeast Reef), four 2015 grouped-distribution identities, and 14/14 Santa Barbara Island physical count-site checks. The existing guard blocks exactly the six affected histories when applied to the actual pool rows. Six targeted guard cases pass, including allowance for original permitted-site sizes and unrelated kelp measurements. These checks support the current hold; they do not certify ecological representativeness or total uncertainty.

## Exact years and count units

Here, a comparison-year is one MPA × species × program × response × post-protection year, not an independent survey, animal, or unique donor-year. Across the six existing annual histories there are **100 post-protection comparison-years: 45 with existing historical diagnostic replacements, 54 requiring exact original size evidence, and one genuine-zero reference year needing no alias size replacement**. Shared donor measurements are reused across KFM/PISCO; do not call these 54 independent missing samples.

| History | Alias role | Existing diagnostic replacements | Exact size measurements still unresolved | Additional detail |
|---|---|---|---|---|
| Harris Point red, KFM | Reference a-k-01 | 2004–2011 (8) | 2012–2019, 2021–2022 (10) | 18 post years |
| Harris Point purple, KFM | Reference a-k-01 | 2004, 2006–2011 (7) | 2012–2019, 2021–2022 (10) | 2005 reference is a verified zero; 18 post years |
| Harris Point red, PISCO | Reference a-k-01 | 2004–2011 (8) | 2012, 2014, 2016, 2018, 2019, 2021, 2022 (7) | 15 post years; no retained 2015 alias use |
| Harris Point purple, PISCO | Reference a-k-01 | 2004–2011 (8) | 2012, 2014, 2016, 2018, 2019, 2021, 2022 (7) | 15 post years; no retained 2015 alias use |
| Santa Barbara Island red, KFM | Protected a-k-14 | 2005–2011 (7) | 2012–2019, 2021–2022 (10) | 17 post years |
| Santa Barbara Island purple, KFM | Protected a-k-14 | 2005–2011 (7) | 2012–2019, 2021–2022 (10) | 17 post years |

`history_years.csv` supplies exact semicolon-separated year sets, including the longer pre-protection donor-use histories and historical full-histogram verification years. The source request's 2012–2022 range is appropriate; absence of 2020 from these retained histories does not justify adding it. The 2015 grouped summaries identify the alias composition but **do not resolve exact sizes for biomass**. The four retained KFM 2015 comparison-years remain in the unresolved measurement column. Historical replacements remain diagnostic, not adopted full-history estimates.

The apparent 100 minus 45 minus 54 discrepancy is explained by Harris Point KFM purple in 2005: its reference biomass is zero in the donor inventory; all **12 raw MBON count rows** for the relevant site/species/year are zero. There is no reference size pool used in that year. Do not label it a missing-size year or invent a donor.

## Physical identity and guard

- Original site 1 is Wyckoff Ledge; 21 is Miracle Mile. All 51 overlapping delivered a-k-01 urchin species-years reproduce original 1+21 distributions. This is a distribution identity, not an individual-record crosswalk; earlier years can have no contribution from 21.
- Original site 14 is Southeast Sea Lion; 37 is Southeast Reef. All 14 a-k-14 urchin species-years in 2005–2011 reproduce site 37 sizes. MBON a-k-14 count densities independently reproduce site 14 original densities for the same 14 species-years. A common delivered label is therefore not proof of physical count–size matching.
- The 2015 NPS spreadsheet independently reproduces the four alias distributions after grouping into its published bins. No individual sizes, nonlinear masses, or dates were inferred from those summaries.
- Original `pBACIPS_PISCO_V10.R` and current data-repo `03_data_import.R`/`01_utils.R` explicitly exclude sites 21 and 37. The original script reads the combined CSV; it does not reveal its construction. The current review guard is consistent with those exclusions and source-specific: it does not automatically reject original permitted-site 14 sizes or Wyckoff kelp.
- `audit_biomass_joint.py` passes canonical `Taxa` values to the guard, and excludes whole affected histories before drawing. Direct evaluation over the actual pool finds exactly the six expected histories. No present-path guard defect was found. This remains a review guard, not a production fix. Reusing it with other species spellings/source identifiers would require explicit normalization and tests.

## Existing holdings and correspondence

All **44 paths** in the prior holdings inventory still exist and match their recorded SHA-256 hashes (42 resolved targets/34 hashes was the prior inventory's deduplication scope, not 44 independent deliveries). The available original individual size table ends in 2011; the 2015 NHSF workbook is binned. The combined CSV has later records but unresolved original physical identity and missing KFM visit dates. Neither cached biomass estimates nor the public processing script supply missing original measurements.

A fresh filename search of the analysis repo, data repo, and Emily's local Kelp MPA project located the legacy table, the 2015 NHSF workbook, and two combined-size CSV paths; no newer original natural-habitat size table or Access database was found by those filename patterns. `candidate_paths.txt` records exact matches. This search used `rg --files --hidden` and respects ignore rules; it is not proof that no unrelatedly named, ignored, remote-only, or inaccessible file exists. The hash recheck separately covers all 44 previously inventoried paths. No new NAS-wide or online search is claimed.

A cached correspondence text search covered **554 files**, with 75 keyword-hit files; locators are in `correspondence_search.json`, without copied message bodies. The Oct 19 Scott source message was decoded and checked: his recorded delivery covered count methods, while invertebrate size frequencies were used by PISCO but not on that recorded shared-file list. The later cached Oct 29 incoming record addresses permissions, not a new size delivery; Oct 31 selection records do not supply original size measurements. The Nov 2 source request is already **sent**, per the saved receipt/status and approved body. Do not repeat the earlier holdings README's draft/not-sent status or send a duplicate request. No later original-size answer was located in these local caches; no live inbox freshness claim is made.

The saved NPS package record lists internal tables with `CanDownload: false`; this is dated acquisition evidence, not a fresh access check. Do not rerun its Access-based processing script without the original database. Keep raw records private pending existing permissions confirmation.

## Recommendations before guide updates — not adopted rules

1. Keep current density eligibility and the six biomass holds. The current guide's main conclusion is correct; no new numerical correction or exclusion is supported by this audit.
2. In `measurement_support_current.csv`, attach the exact affected donor role/site, verified historical replacement years, unresolved exact-measurement years, and a distinct identity-only-2015 flag. Avoid overloading `unsupported_years` if its existing PISCO meaning is missing biomass support; use clearly named alias fields or link a separate alias-year table. Preserve the 2005 zero exception.
3. In `unresolved_issues.csv` U08 and the guide README, distinguish four KFM and two PISCO histories. Add a link to the already-sent Nov 2 request and describe the next step as checking for/receiving the original records, not drafting/sending another request. Keep later identity unknown except where independently checked.
4. Minimal next action: receive the requested original red/purple sizes, original site IDs, survey dates, units and site key if needed. For these six histories prioritize permitted sites **1 and 14**; source records for **21 and 37**, or a documented original-to-combined transformation, will establish how the aliases arose. This prioritization is a recommendation, not a narrower replacement for the already-sent request. Ask for the combined-file construction code only if it is not supplied with the response.
5. Once originals arrive, validate the original site/event keys and reconstruct permitted distributions, retaining the current PISCO threshold and annual coverage unless separately authorized. Recompute physical donor support and apply the guard in the coordinated production rebuild before promotion. Do not split aggregate rows, use binned means as exact masses, zero unknowns, drop years silently, or select corrections by slope.

## Reproduce and evidence

Run `python3 /tmp/donham_swarm_20261102/kfm/check_kfm.py` (pandas, numpy, openpyxl required). It reads sources and writes only beside itself. No production script is executed or imported with write side effects; guard definitions are evaluated without bytecode writes. Results: `checks.json`, `legacy_identity_checks.csv`, `count_identity_checks.csv`, `bins2015_checks.csv`, `guard_checks.csv`, `zero_year_check.json`, `history_years.csv`, and `holdings_recheck.csv`. The initial run reached all assertions but failed serializing a NumPy integer in its final summary; fixed only in scratch, then rerun successfully. No model runs were needed.

Exact numeric input paths are in `script_examined_files.txt`; additional read/search paths and limitations are recorded in `examined_files.md`. `findings.csv` separates verified checks, unresolved evidence, and recommendations. Original source data and private message bodies have not been copied into this report.
