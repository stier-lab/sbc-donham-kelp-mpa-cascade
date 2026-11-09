# Painted Cave size records and Emily’s biomass workflow

November 2, 2026. Original files retained; manuscript, supplement and frozen results unchanged.

## What happened

The original PISCO file identifies 99 size-frequency records as Painted Cave East in 2003. They represent 69 red urchins and 70 purple urchins. Emily’s archived combined size file already contains the conflicting South Point grouping label; her archived site table correctly assigns Painted Cave East to Painted Cave. The current combined file is byte-identical to that archived file. The code that originally assembled this combined file was not found among the archived GitHub scripts, so the exact editing operation that introduced the label remains unknown.

This did **not** send these animals into South Point biomass. Emily’s code joins sizes to site metadata using physical site, grouping label and program together. The inconsistent combination matches neither site. Reproducing that join leaves all 99 records without a comparison assignment. Correcting the alias assigns all 99 to Painted Cave and none to South Point.

## Correction and exclusions

The source adapter now reads the checksum-verified PISCO 1.11 roving file, converts centimeters to millimeters, and obtains the grouping label from the physical site and program in the checked site table. It does not copy the combined file’s erroneous grouping label. All 3,970 historical UCSB measurement records match the original source with multiplicity.

Painted Cave remains excluded under the existing ecological scope because lobster take is permitted. The selection ledger preserves all 332 Painted Cave size rows in the historical source scope, including the 99 corrected rows, but keeps them out of candidate size donors. The other 3,638 historical UCSB records enter the source-derived candidate table; existing physical reference exclusions remain applied by the production importer. Additional 2022–2024 source observations remain outside the approved replacement years. VRG originals remain unresolved.

A new guard rejects a grouping label that conflicts with the physical site before comparison membership is assigned. Tests cover this exact Painted Cave/South Point failure, duplicate records, unknown mappings, source checksums and year/exclusion scope.

## What changes in results

The retained donor-membership table is identical before and after correcting the 99 labels. The correction only restores membership in an excluded comparison. It therefore changes no retained annual biomass inputs or fitted results through this pathway. No model rerun is needed for this identity-only correction; no new fit or total-uncertainty calculation is claimed. The independent-source adapter is also checked against the entire historical measurement multiset, rather than only these 99 rows.

## How Emily estimated biomass

The archived `Github/MPA_Kelp/code/pBACIPS_PISCO_V10.R` and `MPA_code_notes.docx` establish the following:

- The size table combines PISCO, VRG, KFM and LTER measurements. Researcher reference exclusions are applied before selecting donors (script lines 7–38).
- PISCO urchin donors are assigned to comparisons, restricted to diameters of at least 25 mm, and selected by comparison, protected/reference role, year and species (lines 294–336). Despite a nearby comment mentioning site/year/zone, the implemented donor selection does not require the same physical station, survey date, depth zone or monitoring program.
- Each recorded size is repeated according to its count. The code draws the number of animals recorded by the count survey, with replacement, 1,000 times; it converts each sampled animal’s size to mass, sums those masses and averages the totals (lines 338–373). This is the source of cross-program borrowing, not independent replication of each borrowed size sample.
- Red and purple urchins use separate allometric exponents. UCSB lobsters are calculated directly from swath lengths; VRG lobsters use roving size samples (lines 117–212). These measurement pathways must remain distinguishable.
- Her notes explicitly say that many size measurements used for PISCO biomass came from KFM’s separate size-frequency surveys. This corroborates the code; it does not by itself establish that every donor pool is seasonally or spatially representative.

The source-first pipeline preserves these distinctions and logs donor identities. Matching of survey dates, depths, physical stations and shared-donor uncertainty remains a separate review. No legacy pooling or missing-data convention is silently endorsed by verifying the original files.

## Evidence and reproduction

- [Corrected records, including retained original labels](corrected_records.csv)
- [Counts by species](correction_summary.csv)
- [Before/after comparison memberships](membership_changes.csv)
- [Source selection, including excluded and outside-scope rows](pisco_roving_selection_ledger.csv)
- [Verification and archived-source hashes](verification.json)

Run `python3 code/tools/audit_painted_cave_sizes.py` in the analysis repository. Run `make check` in the data repository to rebuild source inputs and test guards. Original files in Emily’s Drive folder are read-only evidence and were not edited.
