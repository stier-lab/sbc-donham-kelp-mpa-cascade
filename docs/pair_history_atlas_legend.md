# Pair-history atlas legend

Annual paired histories for the 145 primary MPA/source/taxon/response series in `data/sumstats_final.csv`, covering all 21 final MPAs. Each of 37 pages contains four facets (A-D, top to bottom). Facet strips identify the full MPA name, source, abbreviated taxon and biomass or density response.

Blue (#0072B2) solid lines with circles denote MPA observations; grey (#707070) dashed lines with triangles denote paired reference observations. Direct labels identify the last observed values; thin dotted label leaders are annotations, not additional observations. Lines connect consecutive observed annual values only. Missing years and missing values break the lines; isolated observations remain visible as markers.

Vertical grey dot-dash lines mark each MPA's actual `MPA_Start` from `data/harmonized/harmonized_site_metadata.csv`. Subtle gold shading marks the complete 1997-1998 and 2014-2016 calendar-year windows (1996.5-1998.5 and 2013.5-2016.5 on the annual axis). These windows provide temporal context; they are not estimates of local disturbance exposure.

The vertical axis shows the actual `mpa` and `reference` proportions supplied in the harmonized paired-response tables, on their shared upstream normalization. Values are plotted unchanged, with the same 0-1.25 axis across all facets and pages; no MPA/reference rescaling, source averaging, physical-unit reconstruction, smoothing or model fit is applied. Values above 1 in the harmonized files are retained. These are normalized responses, not absolute biomass or density, community fractions, or estimates of among-program absolute abundance. The supplied tables have no annual uncertainty estimates, so no confidence ribbons are inferred.

## Coverage and index

All 145 final primary keys have matched annual paired observations; there are no unmatched final series. The selected histories have 0 missing MPA values and 0 missing reference values. The atlas displays 2827 paired annual rows (5654 individual values), spanning 1982-2023. Annual paired counts match the final effect-size `N` for every series. 97 series contain internal missing calendar years.

Exactly 50 harmonized series whose full key is absent from the final primary set are omitted, including rejected series and other responses. No raw-response or external data files are read.

52 selected histories have no Before-period observations (year <= `MPA_Start`) in the supplied paired inputs. Their facets therefore show only the observed later history; no baseline is imputed.

| Source | Response | Series |
| --- | --- | ---: |
| KFM | Biomass | 12 |
| Landsat | Biomass | 16 |
| LTER | Biomass | 7 |
| PISCO | Biomass | 44 |
| KFM | Density | 17 |
| LTER | Density | 4 |
| PISCO | Density | 45 |

`tables/table_s_pair_history_atlas_index.csv` has one row per facet, keyed by `page` and `panel`. `png` and `pdf` are repository-relative paths to individual atlas pages. Include each distinct `png` once, sorted by `page`; QA PNGs in `plots/pair_history_atlas/qa/` are not atlas pages. Pages follow MPA name, source (KFM, LTER, PISCO, Landsat), response (Bio, Den) and taxon (predators through kelp) order.

`dataset_year_min/max` describe the full input-file range; `year_min/max` describe each selected paired history. `n_before` counts finite annual pairs with year <= `MPA_Start`; `n_after` counts pairs with year > `MPA_Start`, matching the harmonized BA convention. `missing_years` lists absent internal calendar years separated by semicolons. The index also records the input dataset, full taxon name, final `N`, analysis type, selected model, missing-value counts, common scale and export specifications.

## Outputs and reproduction

- `plots/fig_s_pair_history_atlas.pdf`: combined multipage vector PDF.
- `plots/pair_history_atlas/page_001.png` onward: 170 x 200 mm pages at 300 dpi.
- `plots/pair_history_atlas/page_001.pdf` onward: matching single-page vector PDFs.
- `plots/pair_history_atlas/qa/contact_sheet_01.png` onward: six-page QA contact sheets, ordered left to right, then top to bottom. Sheet 01 contains pages 001-006, sheet 02 pages 007-012, and so on.
- `plots/pair_history_atlas/qa/page_001_{deutan,protan,desaturate}.png`: accessibility checks.

Run `Rscript code/R/29_pair_history_atlas.R` or `source("code/R/29_pair_history_atlas.R")`.
The script resolves the repository root from its own path, runs in local scope, and does not attach packages or load pipeline objects. Required packages are ggplot2, ggrepel and ragg. magick and colorspace add reproducible QA assets. Figures contain no title or subtitle.

Input MD5 checksums for this render:

| Input | MD5 |
| --- | --- |
| `data/harmonized/harmonized_response_ratios.csv` | `e98d1e3eede110d30e1cfb4c5ec208af` |
| `data/harmonized/harmonized_landsat_rr.csv` | `aa1080e0127e5cacfb357ed073fe53ee` |
| `data/harmonized/harmonized_site_metadata.csv` | `2474b071b2f623c277ba46e94a03b7a1` |
| `data/sumstats_final.csv` | `153842a6bab8d85d1152a4be5a05b29c` |
