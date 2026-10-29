# Sheephead Coefficient Source Check

Checked 12 October 2026. This note preserves the source mismatch found before correction. The subsequent approved correction adopts Hamilton & Caselle (2015) in both PISCO and LTER, without rewriting the original lookup. See `SHEEPHEAD_ALLOMETRY_REFIT_2026-10-12.md` for the rebuild and numerical comparison.

The original PISCO processing used the SPUL row of `spp_attribute_table_downloaded_9-13-22_SHSPUL.csv` in the data-processing repository. Its coefficients are a = 0.0238 and b = 2.8847, with wet mass in g and total length in cm. The row names VRG as the lookup source and Williams et al. 2013 as the length-weight reference. Hamilton et al. 2011 is the separate length-conversion reference; the TL biomass branch does not use that conversion.

Hamilton & Caselle (2015), Methods page 2, report W = 0.0144 TL^3.04 with TL in cm. The LTER route uses that equation. It is not the current PISCO equation. Primary source: https://pmc.ncbi.nlm.nih.gov/articles/PMC4286036/.

Williams et al. (2013), Table 2, printed page 224, give the sheephead row as n = 194, standard length in mm, a = 8.45e-05, b = 2.80 and R-squared = 0.92. That does not directly support the lookup's total-length equation. DOI: https://doi.org/10.3160/0038-3872-112.3.217. The primary PDF was read and the table visually checked through its author-hosted copy, linked from https://www.claisselab.com/publications.

## Upstream Table Trace

The upstream PISCO tables were located through the public [PISCO/CRANE processing code](https://github.com/zoekitchel/PISCO-CRANE-processing/blob/39abf7e97e9b2720e4cb5110700312c39e22e611/scripts/from_PISCO/MLPA_data_merge_fish_2023.R). Its comments link the current and old master sheets. Both were inspected read-only after checking tab metadata and headers, then locating SPUL in the species-code column.

| Source | SPUL location | a | b | Length | Weight | Stated length-weight reference |
| --- | --- | --- | --- | --- | --- | --- |
| Project `..._SHSPUL.csv` | CSV line 150 | 0.0238 | 2.8847 | TL, cm | g | Williams et al. 2013 |
| [Current PISCO master](https://docs.google.com/spreadsheets/d/1h-gAvQWT3sw1PjEau4dI0r31a1w7qUnAx6PaqC_Yk4Q/edit#gid=258316809) | `species_attribute_table`, row 156 | 0.0000845 | 2.8047 | SL, mm | g | Williams et al. 2013 |
| [Old PISCO master](https://docs.google.com/spreadsheets/d/1OsGt4m4n247TTvRzAa1TQB0FrMsBR9Vaa7TdtW_vXQs/edit#gid=258316809) | `MLPA_kelpforest_fishes_life_history_table`, row 154 | 0.0000845 | 2.8047 | SL, mm | g | Williams et al. 2013 |
| [Pinned January 2025 PISCO export](https://github.com/zoekitchel/PISCO-CRANE-processing/blob/39abf7e97e9b2720e4cb5110700312c39e22e611/scripts/from_PISCO/spp_attribute_table_20250113.csv) | CSV line 153 | 0.0000845 | 2.8047 | SL, mm | g | Williams et al. 2013 |
| Hamilton & Caselle 2015 | Methods, page 2 | 0.0144 | 3.04 | TL, cm | g | Published equation |

The upstream exponent 2.8047 is consistent with Williams's printed 2.80, but the paper does not establish the extra digits. The project copy differs in the two coefficients, length units and length type while retaining the Williams citation and the separate Hamilton 2011 conversion fields. These comparisons identify a source mismatch; they do not establish who made the change or recover the master sheet as it stood on the filename's September 2022 download date. A length-unit conversion alone cannot change the exponent from 2.8047 to 2.8847.

The two original project CSV copies are byte-identical (SHA-256 `d8793c327a508ee0c812aef9262308afef71385bc3d96bf494fde5895bd34289`). Their MD5 `6ccfacd18fff07208f510384e914a3b2` agrees with Drive metadata, which gives a last-modified date of 6 September 2023. The legacy `pBACIPS_PISCO_V10.R` reads the sheephead-specific table and describes it as using Hamilton's parameters; its separate LTER section uses 0.0144 and 3.04. Thus the supplied-file discrepancy predates the current rebuild.

Additional checks did not identify a source for the exact 0.0238/2.8847 pair. Miller et al. 2008 does not include sheephead in its length-weight table. Hamilton's two 2011 papers were checked; the separate SL/TL conversion is not evidence for the project's length-weight coefficients. FishBase snapshots 19.04, 21.06 and 23.01 give sheephead a = 0.0289, b = 3, TL (reference 40637, two observations); snapshot 26.06 also has an SL row with a = 0.02903, b = 3 (reference 137577, three observations). Neither is the project equation. Snapshot data were read through the [documented rfishbase data source](https://github.com/ropensci/rfishbase), not inferred from species summaries. Targeted public-code and exact-number searches found no verified source; absence from these searches is not proof that no unpublished fit exists.

## Approved Correction

Hamilton & Caselle (2015) supplies the verified total-length equation already used by LTER: **W(g) = 0.0144 x TL(cm)^3.04**. The upstream PISCO/Williams route is a distinct standard-length equation. The exact published source of the project's 0.0238/2.8847 pair remains unverified.

The original spreadsheet was not edited. Following approval, a sourced in-memory SPUL override adopts Hamilton's equation. A complete staged harmonization, scoped promotion, explicit registry refresh and full refit were carried out; manuscript numbers and figures were refreshed from those results. This is an analytical correction, not a citation-only substitution. The exact origin of the obsolete pair remains unverified, but is no longer the source of the current PISCO biomass analysis.
