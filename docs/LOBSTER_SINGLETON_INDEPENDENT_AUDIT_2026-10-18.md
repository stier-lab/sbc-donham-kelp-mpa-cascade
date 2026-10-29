# Independent adversarial lobster singleton audit

Date: 2026-10-18; amended after scoped promotion on the same date. Read-only investigation; this report is the only file written by this audit. No pipeline, cache-correction script, or export script was executed by this auditor. R computations were performed in memory from source inputs and caches. The source working tree was changing during the investigation; the observations below distinguish the initial pre-correction snapshot, the subsequently observed corrected source cache, and the analysis cache. The initial report incorrectly attached the pre-correction hash to the corrected cache; the amendment below verifies that distinction directly against the archived backup.

## Finding and disposition

**The narrow deterministic correction survives independent attack. Exactly six pre-promotion primary PISCO lobster biomass annual/status records required replacement.** Fourteen positive site/year/zone cache cells are affected overall, but several belong to MPAs/statuses absent from the primary lobster export. Do not promote all corrected cache cells to new primary membership. Adrian reports that scoped promotion is now complete: six raw records and 31 response-ratio rows changed, with all nontarget bytes retained. This amendment verifies cache lineage and the module 43 path change; it does not independently re-audit that promotion or rerun models. The main full pipeline is running separately.

The original source cache from Git HEAD reconstructed all 62 overlapping pre-promotion harmonized primary VRG biomass rows, with maximum absolute error **4.263256414560601e-14 g/m2**. The source working-tree cache observed later in the investigation contained the deterministic singleton corrections. Against that pre-promotion export, it reconstructed 56 unchanged rows and the six corrected values below. The separate analysis cache was not the primary baseline: it matched only six of these 62 pre-promotion annual/status records within 1e-10.

The scientific statement supported is: the inherited singleton sampling behavior underestimated six exported annual/status lobster biomass observations; replacing their underlying singleton zone totals by count times the measured animal's allometric mass corrects that software error while preserving the original aggregation and membership.

## Evidence and lineage

Paths below use `SOURCE = /Users/adrianstier/Donham-Stier-CA-MPA-Data-2026` and `ANALYSIS = /Users/adrianstier/Donham-Stier-CA-MPA-2026`.

| Observed file | MD5 |
| --- | --- |
| ANALYSIS/data/harmonized/harmonized_raw_responses.csv, pre-promotion snapshot | c817b10fac60d714be510fc8547a2685 |
| SOURCE/data/cache/vrg_panint_bootstrap.rds, initial pre-correction snapshot | dfaf16358d517518e16c64db9316f156 |
| ANALYSIS/local_archive/lobster_reconciliation_2026-10-18/before_vrg_panint_bootstrap.rds, verified pre-correction backup | dfaf16358d517518e16c64db9316f156 |
| SOURCE/data/cache/vrg_panint_bootstrap.rds, corrected cache verified at amendment | dccdf45b4dc7fccea4b36461a5cc50dc |
| ANALYSIS/data/cache/vrg_panint_bootstrap.rds | 073c6f13a8a93e902a50aba15021a1dc |
| SOURCE/data/ALL_sizefreq_2024.csv | 97013fa6283a5534cf37207d5486b2f9 |
| SOURCE/data/PISCO/master_site_table_Emilyedit.csv | 53ffa708c3e3cf6899ba4a20e3620fce |
| ANALYSIS/data/PISCO/master_site_table_Emilyedit.csv | 53ffa708c3e3cf6899ba4a20e3620fce |

Both repositories' ALL_sizefreq_2024.csv paths are symlinks to the same Google Drive source file. The site tables have identical content. Both caches contain 1,086 rows and six columns (`site`, `year`, `zone`, `transect`, `biomass`, `count`). The initial source-versus-analysis cache comparison on all five non-biomass fields produced exactly 1,086 rows; 480 finite biomass values differed, and each cache contained 31 missing biomass values. That initial comparison preceded the targeted source-cache correction; its hash must not be assigned to later corrected observations.

At amendment, the corrected source cache is 6,910 bytes with modification time Oct 18 22:01:21, 2026 (local time); the ignored pre-correction backup is 6,960 bytes with modification time Oct 18 22:01:20. The backup MD5 equals the initial source hash. Independently comparing the backup with the current source R object confirms identical non-biomass columns and biomass missingness, with exactly the same 14 biomass cells changed at tolerance 1e-10. These file times and hashes identify the observed snapshots; they are not immutable upstream provenance.

Crucially, I independently read the committed source cache using `git show HEAD:data/cache/vrg_panint_bootstrap.rds`, decompressed it in memory, and compared its R object with the observed working-tree source cache. All five non-biomass columns, including their order and types, are identical. Exactly 14 biomass values differ at tolerance 1e-10, all listed below; biomass missingness is unchanged. This is strong evidence of a targeted correction rather than a fresh bootstrap. It also establishes why using the differing analysis cache would introduce unrelated changes.

## Exact primary replacements

All keys additionally require `source = PISCO`, `taxon_name = Panulirus interruptus`, and `resp = Bio`. Values are g/m2. CSV line numbers include the header and identify the observed pre-promotion analysis input; use keys for implementation.

| MPA | Year | Status | CSV line | Existing biomass | Correct biomass |
| --- | --- | --- | --- | --- | --- |
| Abalone Cove SMCA | 2016 | mpa | 4507 | 0.156735525187098 | 0.5706962100368140 |
| Abalone Cove SMCA | 2017 | mpa | 4512 | 0.259976460638117 | 0.8943045986153320 |
| Abalone Cove SMCA | 2023 | mpa | 4542 | 0.226070423521575 | 0.8943045986153320 |
| Matlahuayl SMR | 2023 | mpa | 4793 | 0.369642664035510 | 1.3414568979229979 |
| Point Vicente SMCA | 2012 | mpa | 4869 | 0.249736529818614 | 0.9481899092784327 |
| Abalone Cove SMCA | 2013 | reference | 5075 | 0.112173155031274 | 0.4398965982624868 |

The existing harmonized density values need no replacement for this correction. Recomputed annual ratios and fitted results can change; six changed raw values do not imply only six changed downstream ratio rows, because normalization and zero-offset rules can depend on the full history.

## Exact singleton pools

I reconstructed the primary size join using the original first occurrence of each `(site, CA_MPA_Name_Short)` in the source site table, excluded the source reference/KFM site lists before joining, and joined sizes on **site + site_new + campus**, retaining all months. Pool identity is **MPA + year + site_status**. The helper expands size-frequency counts and drops missing sizes; singleton means one animal in that expanded vector, not one source record or one distinct size.

These are all 16 VRG singleton pools from 2011 onward. Each source record has count 1. Source row is the one-based data-row index in ALL_sizefreq_2024.csv, excluding its header.

| MPA | Year | Status | Measured site | mm CL | Month | Source row | Cache/export scope |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Abalone Cove SMCA | 2013 | empty | 120_REEF | 90 | 7 | 349 | Existing cached counts zero |
| Abalone Cove SMCA | 2017 | empty | 120_REEF | 50 | 11 | 568 | No matching cached zone |
| Abalone Cove SMCA | 2021 | empty | 120_REEF | 50 | 6 | 1035 | Positive cached zone, empty status excluded from export |
| Abalone Cove SMCA | 2016 | mpa | ABALONE_COVE_KELP_W | 60 | 6 | 2812 | Primary replacement |
| Abalone Cove SMCA | 2017 | mpa | ABALONE_COVE_KELP_W | 70 | 6 | 2980 | Primary replacement |
| Abalone Cove SMCA | 2023 | mpa | ABALONE_COVE_KELP_W | 70 | 6 | 3722 | Primary replacement |
| Abalone Cove SMCA | 2013 | reference | BUNKER_POINT | 80 | 10 | 7309 | Primary replacement |
| Crystal Cove SMCA | 2012 | mpa | CRYSTAL_COVE | 70 | 11 | 8739 | No matching cached zone; excluded MPA |
| Point Dume SMR | 2020 | reference | ESCONDIDO_W | 70 | 12 | 10053 | No matching cached zone; lobster excluded at combination/export |
| Point Dume SMR | 2023 | reference | MALIBU_BLUFFS | 70 | 9 | 20231 | Positive cache cells; lobster excluded at combination/export |
| Matlahuayl SMR | 2023 | mpa | MATLAHUAYL | 70 | 10 | 20399 | Primary replacement |
| Point Vicente SMCA | 2012 | mpa | POINT_VICENTE_W | 100 | 8 | 25149 | Primary replacement |
| Blue Cavern Onshore SMCA | 2011 | reference | SCAI_RIPPERS_COVE | 60 | 9 | 39445 | Positive cache cells; lobster excluded at combination/export |
| Blue Cavern Onshore SMCA | 2020 | reference | SCAI_RIPPERS_COVE | 50 | 10 | 39657 | Positive cache cell; lobster excluded at combination/export |
| Swamis SMCA | 2020 | mpa | SWAMIS | 80 | 10 | 48173 | Positive cache cell; lobster excluded at combination/export |
| Swamis SMCA | 2023 | mpa | SWAMIS | 60 | 10 | 48252 | Positive cache cell; lobster excluded at combination/export |

The broad size table also has UCSB singleton pools. Those are not VRG-bootstrap corrections: UCSB primary lobster biomass follows the separate direct transect-size conversion. A search for singleton pools without intersecting the actual VRG cache keys overstates scope.

## Exact positive cache cells

Mass is `0.001352821 * CL_mm^2.913963113` g, using full primary precision. Single-animal masses are 120.7747810886406 g at 50 mm, 205.4506356132531 g at 60 mm, 321.9496555015195 g at 70 mm, 475.0883261234858 g at 80 mm, and 910.2623129072954 g at 100 mm. Do not multiply these size-frequency lengths by 10: that conversion belongs to UCSB swath lengths recorded in cm.

| Site | Year | Zone | Transect denominator | Count | Original source HEAD total g | Deterministic total g |
| --- | --- | --- | --- | --- | --- | --- |
| 120_REEF | 2021 | MID | 2 | 1 | 32.45571848907334 | 120.7747810886406 |
| ABALONE_COVE_KELP_W | 2016 | INNER | 2 | 1 | 56.42478906735537 | 205.4506356132531 |
| ABALONE_COVE_KELP_W | 2017 | OUTER | 2 | 1 | 93.59152582972220 | 321.9496555015195 |
| ABALONE_COVE_KELP_W | 2023 | INNER | 2 | 1 | 81.38535246776684 | 321.9496555015195 |
| BUNKER_POINT | 2013 | MID | 2 | 1 | 121.14700743377620 | 475.0883261234858 |
| ESCONDIDO_W | 2023 | MID | 2 | 2 | 168.89949411610215 | 643.8993110030390 |
| MALIBU_BLUFFS | 2023 | MID | 2 | 1 | 83.11430022716273 | 321.9496555015195 |
| MATLAHUAYL | 2023 | MID | 2 | 1 | 88.71423936852243 | 321.9496555015195 |
| POINT_VICENTE_W | 2012 | INNER | 2 | 1 | 239.74706862586916 | 910.2623129072954 |
| SCAI_RIPPERS_COVE | 2011 | INNER | 1 | 2 | 108.82575392876561 | 410.9012712265061 |
| SCAI_RIPPERS_COVE | 2011 | MID | 2 | 15 | 811.98502441434346 | 3081.759534198796 |
| SCAI_RIPPERS_COVE | 2020 | INNER | 2 | 1 | 32.52490036319318 | 120.7747810886406 |
| SWAMIS | 2020 | INNER | 2 | 1 | 130.13036606112104 | 475.0883261234858 |
| SWAMIS | 2023 | INNER | 2 | 1 | 49.36369693575755 | 205.4506356132531 |

Every positive singleton working-tree cache total agrees with `count * mass` to floating-point tolerance. SCAI_RIPPERS_COVE MID differs from direct multiplication by only 4.55e-13 g, consistent with equivalent floating-point summation.

## Aggregation attacks and clean negatives

- **Aggregation denominator:** the cache's `transect` is a denominator, not a transect ID to average directly. Primary processing groups by site/year/MPA/status/zone/denominator, sums cache totals, divides by that denominator, averages zones equally within site, then sites equally within MPA/status, and divides biomass by 60 m2. I reproduced this hierarchy. Weighting by pooled animal counts, pooled transects, or the number of sized animals would change the estimand.
- **Zero zones/sites:** keep surveyed zero-count zones and sites. Abalone Cove mpa replacements are mass/(2 transects * 3 zones * 60). Matlahuayl is mass/(2 * 2 * 60). Point Vicente is mass/(2 * 4 zones * 2 sites * 60). Abalone Cove reference 2013 is mass/(2 * 3 zones * 3 sites * 60). These independently recover the six replacement values.
- **Missing-size exclusion:** primary changes missing numeric biomass to zero, then removes entire site-years with positive counts and zero biomass. None of the positive singleton target rows is missing or zero before correction. No new site-year should become eligible through this correction; the observed cache comparison confirms unchanged missingness.
- **Mapping ambiguity:** the correction script checks that candidate cache-row IDs are unique after its many-to-many metadata join. The observed 14-cell result contains no conflicting duplicated cache assignment. Preserve the inherited first site/MPA row rule. Recategorizing sites is a separate scientific change.
- **Pool support:** PANINT rows examined after primary joining had no missing/nonpositive sizes or nonpositive counts. Consequently weighted moments and exact count expansion agree for these inputs. The correction script omits EXCLUDED_KFM_SITES, whereas import applies it, but those IDs are KFM sites and do not affect the UCSB/VRG PANINT pools here. This is a clean negative for the present dataset, not a general exemption from reproducing import filters.
- **Taxon/campus:** pool construction in primary 04 joins all cleaned PANINT size records; VRG restriction applies to the receiving count records, not an extra filter on the source size pool. Retain this distinction even though the identified VRG singleton pools each have one VRG record.

## Remaining traps and limits

1. **Code fix versus cached correction.** Source 01_utils.R now indexes `a[sample.int(length(a), n, replace=TRUE)]`, eliminating numeric singleton expansion. Source 04 still loads the existing cache whenever present unless forced. Fixing only the helper would leave old cached results in place. The observed correction script addresses the cache separately and preserves non-singleton draws; that approach passes this audit.
2. **Random-number stream changes on a forced rebuild.** The index-based helper changes RNG consumption for singleton pools relative to sampling integers 1:CL. A full forced bootstrap can therefore alter later non-singleton pools even with the same seed. Fourteen changed cells describes this targeted cache repair, not a guarantee about full rebuild output. No full rebuild was performed here.
3. **Constant size is not necessarily singleton.** A pool containing several animals all of the same length is already sampled as a vector of length greater than one and does not trigger the R scalar shortcut. Conversely one source row with count greater than one is not singleton. Use the expanded vector length.
4. **Conditional zero variance does not validate sparse biological support.** A singleton has mathematically deterministic mass conditional on its sole observed size, but that does not imply population mean size is known without uncertainty. The correction addresses an implementation error, not sparse-sampling uncertainty or on-/off-transect representativeness.
5. **Secondary differences persist for legitimate reasons.** Analysis 38_coauthor_programme_checks.R restricts size pools to May-October and joins by site; primary 04 uses all cleaned months and joins sizes on site/site_new/campus. Secondary biomass multiplies harmonized density by pooled expected individual mass, while primary VRG biomass retains cache count support, missing-size site-year exclusion, and the hierarchical site/zone means. UCSB primary biomass is direct individual-size conversion. Replacing all primary values with secondary biomass would mix these differences into the singleton fix. The current shared lobster mass helper uses full coefficients, so rounded allometry is not an explanation for its present exact-mass reconstruction.
6. **Earlier reconciliation cache finding, now superseded.** At the initial inspection, analysis 43_lobster_biomass_reconciliation.R read ANALYSIS/data/cache/vrg_panint_bootstrap.rds, and its then-exported status reported max cached-primary discrepancy 0.528959780179221, matching this audit's analysis-cache reconstruction. At amendment, I verified that module 43 instead reads `Sys.getenv("DONHAM_PRIMARY_LOBSTER_CACHE", unset = file.path(dirname(root), "Donham-Stier-CA-MPA-Data-2026/data/cache/vrg_panint_bootstrap.rds"))`. The stale analysis-cache path finding is therefore resolved in the inspected code; any environment override must still identify the intended authoritative cache. I did not rerun the module or inspect results from the running full pipeline. The diagnostic limitation remains: zero-variance singleton MC errors yield NA rather than a useful standardized error, so a quiet MC-z check cannot exonerate the scalar bug. The earlier 51 overlap rows reflected secondary seasonal coverage, not the 62 primary VRG rows checked here.
7. **Export membership remains binding.** Empty statuses are removed in primary 04. Crystal Cove is in EXCLUDED_MPAS. Source 07 removes lobster from the inherited SHEEPHEAD_ONLY_MPAS cohort, including Point Dume SMR, Blue Cavern Onshore, and Swamis. Correcting their cached totals does not authorize adding lobster histories. No claim is made here that the inherited membership is scientifically optimal.

## Reproducibility ceiling

I inspected source 03_data_import.R, 04_pisco_processing.R, 01_utils.R, 07_combine_and_export.R, the proposed correct_lobster_singleton_cache.R, analysis modules 38, 39, 43 and lobster_biomass_utils.R, plus the raw harmonized export, source sizes/site mapping and both caches. Numerical checks independently reconstructed pools and annual values in R without sourcing executable pipeline modules. Git HEAD was used only to read the original cache.

This audit establishes the six raw replacements, the 14-cell cache scope, and the preservation of unrelated cached results. The amendment independently verifies the original backup versus corrected cache and the updated module 43 cache selection. Scoped promotion of six raw records and 31 response-ratio rows is reported complete by Adrian; this audit does not independently validate that promotion, revised response ratios, fitted models, manuscript conclusions, or the running full pipeline. No models were rerun. Snapshot hashes above should be retained with their explicit pre-/post-correction and pre-promotion labels when assessing later files.

## Independent post-promotion numerical check

Appended 2026-10-18 after the combined lobster and PISCO/LTER singleton-urchin promotions. This section supersedes the preceding statement that promotion had not been independently compared: it checks the combined promotion against ANALYSIS/local_archive/lobster_reconciliation_2026-10-18/before_harmonized, before_sumstats_final.csv, and before_tables. Earlier tables and hashes remain historical snapshots. Parent reports core modules 08/09/36/37 completed, with module 15, module 25, and secondary work still running. I did not execute any full module, write inputs, or await those jobs. The only fits run by this check were the focal four-point regression and its 24 permutations in memory; the HAC check inspects exported coefficients and independently recomputes BH adjustment, rather than reconstructing HAC estimates from annual data.

### Promoted data and unchanged support

The archived and current raw tables both contain 6,280 data rows, with identical headers, schema, order, and all identifiers/metadata. Exactly **23 raw rows** change, exclusively in `value`, all for biomass:

| Source | Taxon | Changed raw rows | Changed RR rows |
| --- | --- | --- | --- |
| PISCO | Panulirus interruptus | 6 | 31 |
| PISCO | Mesocentrotus franciscanus | 7 | 17 |
| LTER | Mesocentrotus franciscanus | 8 | 29 |
| LTER | Strongylocentrotus purpuratus | 2 | 2 |
| Total | | 23 | 79 |

The archived and current response-ratio tables both contain 3,341 rows. Exactly **79 RR rows** change, only within `lnDiff`, `mpa`, `reference`, and `Diff`; identifiers, membership, ordering, time variables, and remaining fields are identical. Each of those four numeric columns changes in 78 rows, not 79: 76 rows change all four fields, two LTER purple-biomass rows each retain one normalized component, and Abalone Cove PISCO lobster biomass 2013 changes `mpa` and `reference` while its `lnDiff` and `Diff` stay identical. There are **no missingness transitions** in these four fields. A changed normalized component therefore does not necessarily imply a changed annual lnRR.

In both CSVs, each data line outside the identified changed rows is byte-identical between archived and current text, as is the header; physical line counts are unchanged. Thus the changes are confined to the identified biomass rows, with all density data and other taxon/source rows retained.

The full 14-file archived harmonized directory has nine files unchanged by MD5 and five changed files. The unchanged files are harmonized_landsat_annual_index.csv, harmonized_landsat_rr.csv, harmonized_site_metadata.csv, kelp_monthly_calibration.json, LANDSAT_ANNUAL_INDEX.md, landsat_reference_correction.json, reference_primary_keys.csv, REFERENCE_REGISTRY.md, and sheephead_allometry.json. The changed files are the two promoted CSVs, checksums.md5, mapping_provenance.csv, and reference_sites.csv.

The two provenance-file changes do not change reference membership. All 861 rows and all 22 columns of reference_sites.csv retain their values except `exclusion_definitions_md5`, updated in all 861 rows. Mapping_provenance.csv retains all 14 rows and their roles, kinds, and paths, with four updated MD5 values: harmonized RR input, primary results, exclusion-definition source code, and reference-sites registry. Reference IDs, names, included/excluded flags, coverage, mapping hashes, footprint hashes, and reference-year counts are identical. A blanket claim that every other file is byte-identical would be too broad; the data and mapping content are unchanged, with the stated provenance/checksum refreshes.

### Primary histories and pooled classifications

Both before_sumstats_final.csv and current data/sumstats_final.csv have **145 unique `(MPA, Taxa, Source, Resp)` keys**, with exactly the same key set and all `Primary = Y`. Exactly **five histories change**, and **140 histories retain every parsed field**. The five changed histories are:

| Primary history | Previous mean | Current mean | Previous SE | Current SE |
| --- | --- | --- | --- | --- |
| Abalone Cove SMCA / lobster / PISCO / Bio | -1.669070768372160 | -1.034131783707469 | 0.601004643753685 | 0.444890297756517 |
| Campus Point SMCA / red urchin / PISCO / Bio | -0.260196321225985 | 0.247671033864686 | 0.380903349914306 | 0.265747903080593 |
| Matlahuayl SMR / lobster / PISCO / Bio | 0.231940106973038 | 0.394677877530655 | 0.645402518188181 | 0.543253436513325 |
| Naples SMCA / red urchin / LTER / Bio | 1.230518612026627 | 1.093007666967526 | 0.204890516788950 | 0.180685658645262 |
| Point Vicente SMCA / lobster / PISCO / Bio | 0.401370906575031 | 0.560812084752340 | 0.473734057482547 | 0.378339910803440 |

Only Mean, SE, SD, and CI change in all five; DW_stat and DW_pval change in four. Models, AnalysisType, sample sizes N, flags, locations, and all other fields retain their values. The Campus Point red-urchin history changes sign: unchanged pooled direction must not be misreported as unchanged direction for every individual history.

Comparing archived and current Table 2 confirms **9/9 pooled coefficient signs, 9/9 nominal CI exclusion classifications, 9/9 BH q < 0.05 classifications, all k values, and all MPA counts are unchanged**. The four BH-supported coefficients remain red-urchin biomass, kelp biomass, sheephead biomass, and purple-urchin density. Purple-urchin biomass, lobster biomass, red-urchin density, lobster density, and sheephead density remain unsupported at BH 0.05. Recomputing BH across all nine current raw P values agrees with Table 2 q values to 4.03e-16. Biomass coefficients and precision change; density coefficients are identical, but some density q values change because the nine-test BH family also contains the changed biomass tests. These are preserved classifications, not preserved estimates.

### Focal Table 3 regression and exact permutation

I rebuilt the focal biomass relationship directly from current primary histories, using the primary Table 3 aggregation: equal means across source-specific effects within MPA/taxon/response, response variance `sum(SE^2)/n_sources^2`, then REML `metafor::rma` with Knapp-Hartung inference. The same **four MPAs** occur before and after promotion: **Abalone Cove SMCA, Campus Point SMCA, Naples SMCA, and Point Vicente SMCA**. The primary reference registry and its primary-key file retain their content as described above; no reference reassignment is indicated.

The independently fitted lobster-biomass slope is **-0.5616019887499334**, SE **0.1152898433697170**, nominal 95% CI **[-1.0576541480370345, -0.0655498294628323]**, P **0.0396531878371408**, reproducing current Table 3 (-0.5616). The nominal interval excludes zero.

I enumerated every one of the **24 = 4!** permutations of the four response effect labels, holding predictor values and response-variance positions fixed, and refitted the same model. All 24 fits succeeded. **Three** absolute slopes are at least as large as the observed absolute slope (tolerance 1e-12), giving the two-sided exact-enumeration **P = 3/24 = 0.125**, without a Monte Carlo plus-one adjustment. This matches the inspected permutation export. The exact enumeration is conditional on this particular label-exchangeability scheme; it neither resolves shared-reference dependence nor propagates predictor estimation uncertainty. It does not support a conventional P < 0.05 claim despite the nominal parametric interval.

Qualification: the inspected table_s_primary_pathway_membership.csv still contained old `verified_slope = -0.421545450712644`, old table3_slope -0.4215, and old primary/Table 3 hashes. It is a stale downstream snapshot at this check, not evidence against the current regression reproduced above. Parent should refresh it with the remaining outputs before manuscript integration. I did not use it as the current numerical input.

### Purple-urchin biomass HAC qualification

The current module 36 input manifest matches all six inspected current input/code MD5 values, including current primary histories and harmonized RR data. Recomputing BH separately within each exported nine-coefficient scenario agrees with its q values to at most 6.67e-16.

| Scenario | Estimate | Nominal 95% CI | Raw P | Nine-test BH q |
| --- | --- | --- | --- | --- |
| Retained primary | -0.403981266476213 | [-0.825655823539694, 0.017693290587268] | 0.060147890483349 | 0.093997533527630 |
| Calendar HAC automatic | -0.485209416554644 | [-0.909454069045695, -0.060964764063593] | 0.025560018387222 | 0.046008033097000 |
| Calendar HAC lag 3 | -0.500303010521991 | [-0.921729378464972, -0.078876642579010] | 0.020628173121859 | 0.037130711619346 |
| Calendar HAC lag 1 | -0.464146729584557 | [-0.888504207762653, -0.039789251406461] | 0.032486874106056 | 0.058476373390902 |

Thus automatic and lag-3 scenarios cross BH 0.05 while the retained primary and lag-1 scenario do not. **Lag 1 is nonsignificant after BH correction, but its nominal CI excludes zero and raw P is below 0.05.** Describe the change as sensitivity-dependent adjusted support, not a uniformly robust primary detection. Source module 36 floors history-level synthesis SEs at their primary values; changing synthesis weights and fitted heterogeneity can still alter pooled coefficients and support. The HAC values above are independently checked table arithmetic and input currency, not a fresh annual-level HAC reconstruction. Temporal module 37 and the remaining running secondary analyses were not independently audited here.

### Post-promotion snapshots

| Current file | MD5 at this check |
| --- | --- |
| data/harmonized/harmonized_raw_responses.csv | 95db0b520f097497cc6f189a62329dec |
| data/harmonized/harmonized_response_ratios.csv | 4ca8444a3bf4ab0fc720594cde63307e |
| data/sumstats_final.csv | fe918c4ae52e64f5382a82a599c2dacb |
| tables/table_02_meta_analysis.csv | d4d363034b713d3f11d28f5cba3c715c |
| tables/table_03_cross_taxa_meta_regression.csv | 31f0f895e75fabc8811af86fc3f70839 |
| tables/table_s_primary_hac_pooled.csv | aba5f5cca58dec26d94008889215fe55 |

These are observed files during an active parent pipeline, not certification that every downstream artifact has been refreshed. This appended check changes only this report; manuscript integration remains with the parent.
