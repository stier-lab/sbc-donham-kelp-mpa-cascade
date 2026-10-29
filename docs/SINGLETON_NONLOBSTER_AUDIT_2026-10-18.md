# Singleton Nonlobster Biomass Audit, 2026-10-18

## Finding

**Further primary corrections ARE required beyond VRG lobster.** Two currently retained primary red-urchin biomass histories contain positive-count observations drawing from exactly one expanded individual:

- `Campus Point SMCA|M. franciscanus|PISCO|Bio`: MPA side, 2016, 2017, 2018, 2022, 2023; 28 positive-count zone records.
- `Naples SMCA|M. franciscanus|LTER|Bio`: reference side, 2010; four positive-count transect records.

**KFM is a clean negative:** no KFM cache record matches a singleton pool, including all retained primary histories. No other primary PISCO/LTER urchin biomass key has a positive-count singleton match. Neither affected red-urchin key appears in the current `table_s_primary_pathway_membership.csv`; this does not exempt them from correction in the broader primary results.

This audit was read-only except this report. No processing script, cache, input, harmonized output, membership, site classification, or missingness designation was changed. Calculations below are in-memory checks, not a rebuilt analysis or revised model estimates.

## Scope and Provenance

Paths below use:

- `D = /Users/adrianstier/Donham-Stier-CA-MPA-Data-2026`
- `A = /Users/adrianstier/Donham-Stier-CA-MPA-2026`

Inspected original processing in `D/code/R/03_data_import.R`, `04_pisco_processing.R`, `05_kfm_processing.R`, `06_lter_processing.R`, and the helper/constants in `01_utils.R` and `00c_analysis_constants.R`. Read `D/data/ALL_sizefreq_2024.csv`, the PISCO master site table, MBON geolocation, all three urchin caches in **both repositories**, `A/data/harmonized/reference_primary_keys.csv`, `harmonized_raw_responses.csv`, `A/data/sumstats_final.csv`, and the primary pathway membership table.

The raw size-frequency file has 335,498 rows. Target-taxon import filtering selects 140,917 rows including lobster; existing site exclusions leave 110,613. The reconstructed merged urchin pool has 107,020 rows. Its counts contain no NAs, negative values, or fractional values. No downloads or authentication were needed.

The reconstruction was independently checked by evaluating the original import filter/recode expressions from `03_data_import.R:82-153` and initial urchin subset/merge expressions from `04_pisco_processing.R:520-530` in memory, followed by the original metadata backfill. `all.equal()` confirmed identical unfiltered and >=25 mm pool summaries to the audit reconstruction. No full pipeline scripts were sourced.

## Exact Pool Definition

`bootstrap_biomass()` at `D/code/R/01_utils.R:355` expands each selected row with `rep(size, count)`, concatenates these vectors, and removes NA sizes. The vulnerable call is `sample(a, n, replace = TRUE)`. The relevant predicate is therefore:

```r
a <- unlist(lapply(t2, function(j) rep(size_freq_table$size[j],
                                     size_freq_table$count[j])),
            use.names = FALSE)
a <- as.numeric(na.omit(a))
affected <- count > 0 && length(a) == 1L && a[1] >= 1
```

Counting one source row or one distinct size is insufficient. A single size row with count seven yields seven entries and does not trigger numeric-scalar sampling semantics. Conversely, many rows with NA sizes can leave exactly one usable individual after expansion and NA removal.

The original import includes VRG/UCSB `STRPUR`/`MESFRA`, KFM full-name urchin codes, and LTER `SFL`/`SFS`/`SPS`/`SPL`. It applies the existing `EXCLUDED_KFM_SITES` and `EXCLUDED_REFERENCE_SITES` lists, recodes taxa, merges PISCO metadata by `(site, site_new, campus)`, and backfills only missing MPA/status metadata using MBON `site_id`. The site table is deduplicated by `(site, CA_MPA_Name_Short)` exactly as in script 04.

All three bootstrap lookups use **MPA name + status + year + species**, not size-source campus or the individual count site's identity. PISCO uses the combined pool filtered to size >=25 mm (`04:574-576`); KFM uses the unfiltered combined pool (`05:295`); LTER uses that unfiltered pool after removing NA size rows (`06:183-184`). The helper's later NA removal makes KFM/LTER usable expanded counts identical for these cells.

Consequently, PISCO counts at IV_REEF_E/W legitimately match size records from LTER site `a-l-07` under the existing rules. This audit preserves those rules and makes no site reclassification.

## Cache Census

Counts below apply to both repositories' cache keys. The older A KFM/LTER caches each contain one entirely NA placeholder row, omitted from the census. The numeric biomass contents differ between repositories, as discussed below.

| Program | Valid cache rows | Positive-count rows | Rows matching singleton pools, including zeros | Positive-count singleton rows | Positive singleton rows in retained primary years |
| --- | ---: | ---: | ---: | ---: | ---: |
| PISCO | 8,940 | 7,385 | 96 | 50 | 28 |
| KFM | 1,484 | 1,465 | 0 | 0 | 0 |
| LTER | 1,156 | 729 | 72 | 26 | 4 |

The PISCO >=25 mm pool contains 29 singleton cells, of which 25 are requested by cache records and 19 have positive-count cache records. The unfiltered pool contains 34 singleton cells; KFM requests none, and LTER requests 12, all with at least one positive-count record. These are pool cells, distinct from the multiple zone/transect records drawing from each cell.

Primary membership was checked by exact `MPA|Taxa|Source|Bio` match to `reference_primary_keys.csv`, then exact membership of the cache year in `observed_years`. Both affected keys also appear with `Primary == "Y"` in current `sumstats_final.csv`, with `AnalysisType == "CI"` and `Model == "Mean"`.

## Affected Primary Keys

### PISCO Campus Point Red Urchin

Every row below is the pool key `(Campus Point SMCA, mpa, year, MESFRAAD)`, with exactly one expanded individual. Source row numbers are **1-based data rows excluding the CSV header** in `D/data/ALL_sizefreq_2024.csv`. All six source records in this section and the Naples section have count exactly one.

| Year | Singleton size (mm) | Size-source campus/site | CSV data row | Positive cache records / all matching records |
| --- | ---: | --- | ---: | ---: |
| 2016 | 40 | LTER / a-l-07 | 316245 | 6 / 6 |
| 2017 | 70 | LTER / a-l-07 | 317394 | 5 / 6 |
| 2018 | 90 | LTER / a-l-07 | 318301 | 5 / 6 |
| 2022 | 100 | LTER / a-l-07 | 321411 | 6 / 6 |
| 2023 | 65 | LTER / a-l-07 | 322313 | 6 / 6 |

Exact cache keys are the Cartesian product of the listed years, sites `IV_REEF_E` and `IV_REEF_W`, species `MESFRAAD`, and zones `INNER`, `MID`, `OUTER`, with two zero-count exceptions. Counts below fully identify the affected positive records; zero entries return before sampling and are unaffected.

| Year | IV_REEF_E counts INNER / MID / OUTER | IV_REEF_W counts INNER / MID / OUTER |
| --- | --- | --- |
| 2016 | 7 / 21 / 99 | 12 / 60 / 5 |
| 2017 | **0** / 28 / 55 | 11 / 17 / 50 |
| 2018 | 5 / 18 / 16 | 12 / 57 / **0** |
| 2022 | 26 / 18 / 53 | 15 / 117 / 59 |
| 2023 | 2 / 20 / 58 | 18 / 16 / 8 |

For a singleton of diameter `d`, the correct unnormalized total is deterministic: `count * 0.00059 * d^2.917`, with bootstrap SE zero under this empirical pool. The D cache has totals only 0.2501-0.2713 times those deterministic values. The A PISCO cache also contains finite, depressed values for all 28 positive records, rather than corrected singleton totals.

Applying the existing normalization (divide by two transects, average zones/sites, divide by 60 m2) reproduces the **current harmonized MPA annual biomass values exactly** from D cache contents:

| Year | Current retained MPA biomass | Deterministic singleton correction |
| --- | ---: | ---: |
| 2016 | 2.108953533 | 7.876960676 |
| 2017 | 8.332160793 | 31.805071334 |
| 2018 | 11.567886405 | 44.408762448 |
| 2022 | 41.956718780 | 161.031931817 |
| 2023 | 5.087271850 | 19.415456544 |

These are annual biomass checks, not newly fitted response ratios. All five years occur in the retained primary history, on the After side. The reference observations in those years were not marked affected merely because they share an MPA/year.

### LTER Naples Red Urchin

Exact pool key: `(Naples SMCA, reference, 2010, MESFRAAD)`. The sole usable individual is 24.9 mm, from LTER site `a-l-03`, original class `SFS`, CSV data row **311865** (2010-07-22). Before LTER's explicit NA-size filter, the same cell has 11 source rows; after expansion/NA removal, exactly one individual remains. The other rows do not make it a nonsingleton pool.

Exact positive cache keys are `(site = a-l-03, year = 2010, y = MESFRAAD, transect)` as follows. All have `quad = 12`, `area = 1`, MPA `Naples SMCA`, status `reference`.

| Transect | Count | D cached total biomass | Correct deterministic total biomass |
| --- | ---: | ---: | ---: |
| c-l-00010 | 26 | 45.159032996 | 181.358892525 |
| c-l-00009 | 4 | 7.021426008 | 27.901368081 |
| c-l-00008 | 3 | 5.258118516 | 20.926026061 |
| c-l-00005 | 2 | 3.487127421 | 13.950684040 |

The same pool is matched by zero-count transects `c-l-00006` and `c-l-00007`, which are unaffected. D cached positive totals are 0.2490-0.2517 times their deterministic values. Averaging the six transect totals divided by 12 reproduces current retained reference biomass **0.84619034641177** exactly. The deterministic corrected value is **3.39079125981959**. Both status observations are currently present in the harmonized 2010 history; the MPA-side value is 10.5356728 and is not assigned this reference singleton pool. Year 2010 is explicitly retained in the primary registry, on the Before side.

Important cache distinction: the older A LTER cache stores NA biomass for these four positive-count keys, whereas D stores finite bug-affected values. Across all 26 positive singleton LTER records, D has 26 finite values; A has seven finite values and 19 NAs. This is an observed cache-content difference, **not** an inference that any currently retained observation is missing. Current harmonized Naples 2010 reference biomass matches D, not the A cache.

## KFM Negative and Constant-Size Controls

No KFM cache key matches an expanded pool of length one, whether primary or otherwise. The smallest nonempty pool for positive KFM counts has four expanded individuals (Gull Island SMR MPA purple urchin 2005; Anacapa Island SMR 2003 MPA red urchin 2022). Restricting to exact retained primary years, the minimum is 11. No KFM primary correction is required for this singleton mechanism.

Three PISCO cache pool cells have only one unique size but multiple expanded individuals and are **not** affected:

| MPA | Status | Year | Species | Size (mm) | Source rows | Expanded individuals |
| --- | --- | ---: | --- | ---: | ---: | ---: |
| Campus Point SMCA | mpa | 2014 | MESFRAAD | 70 | 1 | 7 |
| Campus Point SMCA | reference | 2018 | MESFRAAD | 80 | 1 | 4 |
| Naples SMCA | reference | 2013 | STRPURAD | 40 | 2 | 49 |

LTER also has numerous constant-24.9-mm pools with multiple individuals. For example, Naples MPA red urchin 2010 has 22 expanded individuals across six source rows, and Naples reference purple urchin 2010 has 34 across 15 source rows. Neither is a numeric-scalar singleton. The meaning of the stored 24.9-mm class value was not changed or reassessed here.

An in-memory base-R demonstration with `set.seed(42)` gave range 1-24 for `sample(24.9, 1000, replace = TRUE)`, versus exactly 24.9 for every draw from the actual singleton. This verifies that noninteger 24.9 is vulnerable too.

## Other Positive Singleton Cache Keys

These keys are affected cache calculations but **not** currently retained primary urchin biomass histories. Sizes/counts come from the reconstructed program-specific pools. Listing them does not recommend adding, removing, or reclassifying histories. Species codes are R = `MESFRAAD`, P = `STRPURAD`.

PISCO positive singleton keys outside the two primary findings:

| MPA | Status | Year | Species | Size (mm) | Positive cache rows | Count-site IDs |
| --- | --- | ---: | --- | ---: | ---: | --- |
| Arrow Point to Lion Head Point SMCA | mpa | 2008 | R | 126 | 1 | SCAI_LION_HEAD |
| Arrow Point to Lion Head Point SMCA | mpa | 2012 | R | 106 | 2 | SCAI_INDIAN_ROCK |
| Cat Harbor SMCA | reference | 2012 | R | 112 | 3 | SCAI_IRON_BOUND_COVE |
| Cat Harbor SMCA | reference | 2020 | P | 58 | 3 | SCAI_IRON_BOUND_COVE |
| Crystal Cove SMCA | mpa | 2019 | P | 52 | 2 | CRYSTAL_COVE |
| Dana Point SMCA | mpa | 2011 | P | 62 | 1 | DANA_POINT |
| Dana Point SMCA | reference | 2022 | R | 82 | 1 | SAN_MATEO_KELP |
| Dana Point SMCA | reference | 2022 | P | 36 | 1 | SAN_MATEO_KELP |
| Farnsworth Onshore SMCA | mpa | 2023 | R | 111 | 2 | SCAI_CHINA_POINT |
| Farnsworth Onshore SMCA | reference | 2023 | R | 113 | 1 | SCAI_SALTA_VERDE |
| Matlahuayl SMR | mpa | 2019 | R | 69 | 1 | MATLAHUAYL |
| Matlahuayl SMR | reference | 2012 | R | 135 | 1 | CHILDRENS_POOL |
| Matlahuayl SMR | reference | 2012 | P | 69 | 1 | CHILDRENS_POOL |
| Swamis SMCA | reference | 2011 | P | 77 | 2 | LEUCADIA |

The other six PISCO cache-matched singleton cells have only zero counts and never reach sampling: Arrow Point to Lion Head Point SMCA/mpa/2023/R; Blue Cavern Onshore SMCA/mpa/2023/P; Cat Harbor SMCA/mpa/2022/R; Point Dume SMR/reference/2022/R; South La Jolla SMR/empty stored status/2022/P; Swamis SMCA/reference/2022/R. Empty stored status was preserved literally.

LTER positive singleton keys outside the Naples primary finding:

| MPA | Status | Year | Species | Size (mm) | Positive cache rows | Positive count-site IDs |
| --- | --- | ---: | --- | ---: | ---: | --- |
| Campus Point SMCA | mpa | 2008 | P | 24.9 | 1 | a-l-07 |
| Campus Point SMCA | mpa | 2016 | R | 40 | 1 | a-l-07 |
| Campus Point SMCA | mpa | 2017 | R | 70 | 1 | a-l-07 |
| Campus Point SMCA | mpa | 2018 | R | 90 | 1 | a-l-07 |
| Campus Point SMCA | mpa | 2022 | R | 100 | 1 | a-l-07 |
| Campus Point SMCA | mpa | 2023 | R | 65 | 1 | a-l-07 |
| Campus Point SMCA | reference | 2002 | P | 24.9 | 4 | a-l-06, a-l-08 |
| Campus Point SMCA | reference | 2004 | R | 24.9 | 3 | a-l-06, a-l-08 |
| Campus Point SMCA | reference | 2009 | R | 84 | 2 | a-l-08 |
| Campus Point SMCA | reference | 2010 | R | 24.9 | 3 | a-l-06, a-l-08 |
| Naples SMCA | reference | 2001 | P | 24.9 | 4 | a-l-03 |

The MBON table explicitly assigns `a-l-03` to Naples/reference; `a-l-07` to Campus Point/mpa; and `a-l-01`, `a-l-06`, `a-l-08` to Campus Point/reference. No alternate geographic assignment was inferred.

## Required Disposition

Correct singleton handling in the shared helper and refresh affected **PISCO and LTER urchin** biomass caches and downstream harmonized histories/results in addition to the separate lobster work. A helper edit alone will not correct preexisting caches. Refit the two affected primary red-urchin histories and regenerate dependent primary summaries. Current primary estimates are Campus Point/PISCO red-urchin biomass mean -0.2601963 (SE 0.3809033) and Naples/LTER red-urchin biomass mean 1.2305186 (SE 0.2048905); these are identification values, not corrected estimates.

No KFM singleton correction is indicated. No retained primary purple-urchin history has a positive singleton match. The present pathway-membership table contains neither affected red-urchin key; any broader rebuild should verify that membership and resulting products rather than assume all outputs change. This audit does not quantify refitted effects, alter uncertainty policy, or generalize the negative KFM result to other biomass bugs. Changing sampling behavior can also change later RNG draws unless the main repair preserves the previous RNG consumption; incidental Monte Carlo differences are separate from direct singleton exposure.

## Snapshot Hashes

SHA256 values identify the inspected snapshot; D and A cache files are not interchangeable. Hashes were read after the numerical checks.

| File | SHA256 |
| --- | --- |
| D/data/ALL_sizefreq_2024.csv | `711e2ff55973a3592ac7aa6a00f37cfbe33c482d9684eaeecc41519fc2552741` |
| D/data/MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv | `c5a21d53b5f630348c5894f209627b25b3a190b0c28d3eef791bd2982dc1c059` |
| D/data/cache/pisco_urchin_bootstrap.rds | `294afb7986274aef44cd737efcc71fa21e62e2362c0bc05c6860ae40f98ce365` |
| D/data/cache/kfm_urchin_bootstrap.rds | `2ea64da70501b275c7a2bfc62caa3dce0ccf731d6fdae6fdd629ccbef66c1630` |
| D/data/cache/lter_urchin_bootstrap.rds | `7cc1a983b529d7ec87962273002b72511775ed9b34d5be03fa105d66e50582df` |
| A/data/cache/pisco_urchin_bootstrap.rds | `247db2934d39ac416062a7d31b13551a7358f0c97ef9e0e83140a2baeaff8779` |
| A/data/cache/kfm_urchin_bootstrap.rds | `3608cba3a56972240f0ee470bf66f07ddeea770133291afc63b02a6febb9b7a7` |
| A/data/cache/lter_urchin_bootstrap.rds | `55db7aa8e8b5999004c795a882de52ff86b8ccb5bedf035d47373c09eb00fc3d` |
| A/data/harmonized/reference_primary_keys.csv | `bbac189ec70df814b0c97142ea259c00254394ff72f062e859d00b965e24cd12` |
| A/data/harmonized/harmonized_raw_responses.csv | `7f5e3b552bc93173534927ecb6ab6eb5908ad92a6af3a47053644e2a34d1adfd` |
