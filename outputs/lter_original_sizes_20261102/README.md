# Original SBC LTER size source located

The NAS is mounted and the original LTER benthic survey tables are present. The previously reported LTER source-availability gap is resolved. A source-backed production adapter and explicit treatment of group-average sizes/missing sizes still need integration; the analysis has not been changed by this check.

## Verified source files

- SBC LTER package `knb-lter-sbc.19.32`: `Annual_Quad_Swath_All_Years_20250903.csv`, 459,200 rows.
- SBC LTER package `knb-lter-sbc.19.33`: `Annual_Quad_Swath_All_Years_20260923.csv`, 478,384 rows.
- Both NAS files match their official EML MD5 checksums. A compressed copy of version 33 in the SSWD/kelp source archive also matches after decompression.
- NAS locations and verification: [source table](nas_sources.csv).

The shared catalogue initially provided metadata while the NAS was unmounted. After mounting, the actual files were checked. Public unauthenticated PASTA access returned 403 and the local EDI token was expired; no expired credential was used. The public DataONE route also returned version 33 successfully. No new data-access permission is implied by these checks.

## Comparison with our dataset

The original source reconstructs all **2,706** focal LTER size rows from the combined file on physical site, year, species, count and size, preserving repeated rows. This uses the existing six-site scope, positive counts, years through 2023, original `SFL/SFS/SPL/SPS` codes, centimeters × 10, and `-99999` as missing. Sizes are compared after rounding to six decimal places in millimeters to accommodate decimal serialization; this is not a claim of byte identity.

**2,702 rows also agree on date.** Four Arroyo Quemado purple-urchin rows have July 19, 2022 in the combined file but July 27 in both verified 2025 and 2026 source releases. Counts and sizes are unchanged: 8 animals at 30 mm, 16 at 40 mm, 9 at 35 mm and one at 24 mm. We have not established which earlier version introduced the date difference.

Annual site/species size frequencies and count-weighted mean masses therefore do not change from this comparison. Visit-level matching could be affected by the eight-day difference and should use a documented source-version decision. No model refit or primary-date correction was applied.

[Species comparison](comparison.csv) · [Combined unmatched rows](combined_unmatched.csv) · [Original unmatched rows](original_unmatched.csv) · [Site crosswalk](site_crosswalk.csv) · [Verification](verification.json)

## What the measurements mean

Official EML defines `SIZE` as the **average size of individuals counted**, in centimeters. `COUNT` is the number of individuals in the survey area. These are survey-group size observations, not necessarily separately measured individual animals. Preserving the original transect, quadrat and side fields will make this distinction visible in the pipeline.

Of the 2,706 rows, **1,809 have size values and 897 have missing sizes** in both sources. Missing size does not mean zero count or zero biomass. Small-class values of 2.49 cm already occur in the original source; they were not introduced by a new 24.9-mm substitution in this check. Their measurement/class-representative interpretation still matters for uncertainty.

The catalogue also contains separate urchin-size surveys in packages 26 and 175. Package 26 includes experimental treatments. Their existence does not authorize mixing them into the monitoring size pool; they are unnecessary for explaining the 2,706 existing rows.

## Additional source lead

A standalone PISCO file exists at `SBC-Data-Catalogue/pisco-data-catalogue/data/raw/PISCO_kelpforest.1.11/PISCO_kelpforest_sizefreq.1.11.csv`. It has 13,370 UCSB rows (2003–2024) and 5,770 UCSC rows (2003–2011), with no VRG rows. This is a source candidate, not yet a verified replacement for the PISCO portion of the combined file. Earlier file/email searches did not cover this shared-cache location adequately.

Reproduce the LTER comparison with `python3 code/tools/check_lter_original_sizes.py`. Raw files, numerical inputs, manuscript and supplement remain unchanged. No emails or push.
