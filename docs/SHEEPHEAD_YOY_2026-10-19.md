# Sheephead Size-Proxy Diagnostic, 2026-10-19

## Current Result

**Baseline guards passed for all 816 retained annual raw cells, all 30 retained PISCO sheephead primary contrasts and their SEs, and all nine published joint estimates/SEs/BH values.** The Rscript entry-point run completed on the current input fingerprints. All 45 protected files remained byte-identical, and subsequent output tests and current-hash checks passed. No primary model or input was changed.

Unique physical SPUL records contributing to retained PISCO histories contain 20,877 fish across 17,762 rows. The <10 cm proxy comprises 946 fish (4.5313% of count) and 0.04770% of calculable raw biomass. The <=10 cm proxy comprises 1,356 fish (6.4952%) and 0.10446% of biomass. Exactly 10 cm contributes 410 fish; it is why the scenarios differ. One retained positive-count record, representing three fish, lacks a valid length. Those fish are not allocated to either size proxy.

The broader UCSB/VRG campus-filtered file contains 29,259 SPUL rows representing 35,079 fish: 1,417 fish <10 cm (4.0395% of count; 0.04528% of biomass) and 2,132 fish <=10 cm (6.0777%; 0.10298% of biomass). Its 715 fish at exactly 10 cm and eight invalid-length positive records (26 fish) differ from the retained-history subset. The original full-file 29,792 SPUL rows also include 533 UCSC rows that the pipeline excludes. There are no missing SPUL counts or zero-count SPUL rows in these two raw-share populations; the source's complete species roster subsequently adds survey zeros.

The fixed-route joint sheephead results are:

| Response | Scenario | Estimate | 95% CI | Nine-family BH q |
| --- | --- | ---: | --- | ---: |
| Biomass | Baseline | 1.168170 | 0.770225 to 1.566114 | 0.000001127 |
| Biomass | Remove <10 cm | 1.168944 | 0.771023 to 1.566865 | 0.000001108 |
| Biomass | Remove <=10 cm | 1.171156 | 0.772332 to 1.569981 | 0.000001118 |
| Density | Baseline | -0.319458 | -1.445061 to 0.806146 | 0.572539 |
| Density | Remove <10 cm | -0.316787 | -1.433935 to 0.800361 | 0.572867 |
| Density | Remove <=10 cm | -0.300816 | -1.416940 to 0.815308 | 0.591982 |

**Neither cutoff changes any of the nine BH detection decisions.** The largest pooled sheephead change is +0.002987 for biomass and +0.018642 for density (both <=10 cm). This supports insensitivity of these retained joint results to these two PISCO size exclusions, conditional on the unchanged sampling and model routes; it does not independently verify fish age.

Small pooled changes do not imply every local contribution is small. The maximum absolute history-contrast change is 0.010730 for biomass and 0.164422 for density under <10 cm; under <=10 cm it is 0.092711 for biomass and 0.226384 for density. The latter two occur at Long Point SMR. No retained history changes sign. Individual annual cells lose as much as 18.64%/42.67% of biomass and 81.48%/88.89% of density under the strict/inclusive scenarios, respectively. The strict biomass maximum and both density maxima occur at Abalone Cove SMCA, 2016, MPA side; the inclusive biomass maximum occurs at Long Point SMR, 2023, reference side. The annual table retains these cells and their original zeros; the raw pooled biomass share should not be used to characterize all local histories. Differences at floating-point precision can occur when reconstructing rounded CSV values and are not substantive contributions.

## Scope

This is an optional PISCO raw-input diagnostic and fixed-route sensitivity. It does not edit primary data, caches, models, runners, or manuscripts. The two scenarios remove SPUL observations with finite positive total length **<10 cm** or **<=10 cm** by setting their counts to zero. Both preserve the complete original survey/species roster. These are metadata-derived size proxies, not independently aged young-of-year fish. Neither retained >=10 cm nor >10 cm fish should be called adult-only.

The authoritative local lookup is `data/PISCO/spp_attribute_table_downloaded_9-13-22_SHSPUL.csv`: the SPUL row has `yoy_size = 10`, `juv_size_cutoff = 10`, and `yoy_code = SPULYOY`. The strict/inclusive boundary convention is unverified, so both scenarios are reported. The raw file contains 29,792 SPUL rows and no SPULYOY rows, with observed nonmissing lengths from 1 to 100 cm. These full-file counts include UCSC; the source pipeline subsequently retains UCSB and VRG only. Eight positive-count SPUL records (26 fish) have invalid/missing lengths in the full file; no age or size is imputed for them.

KFM is not filtered: the available count inputs do not provide individual lengths or age classes. The parent's LTER inventory (39 positive-count fish <10 cm, 28 at exactly 10 cm, total positive count 1,198, all positive counts sized) does not establish that PISCO's metadata cutoff classifies LTER fish ages. No LTER age-class sensitivity is substituted for an independent classification.

## Reproduction and Guards

Run `Rscript --vanilla code/R/44_sheephead_yoy_check.R` from the analysis repository. An optional argument supplies the sibling data-processing repository; the default resolves it relative to the analysis root. The module executes under Rscript and can be sourced to load definitions without executing. Missing raw files or sibling source produce an explicit `skipped_raw_unavailable` status instead of new model claims. A complete status is required before treating any sensitivity table as current.

The module parses the source `04_pisco_processing.R` fish block from `fish_att_path` through the Sheephead-records message and evaluates it in an isolated environment with the source constants, helper definitions, and the original first-unique `(site, MPA)` site table. It retains the verified count-weighted sheephead allometry. Cutoffs act after the original campus filter and before merging/aggregation. Sensitivities reuse the baseline's per-row allometric results, zeroing the selected finite-length SPUL biomasses and counts; all aggregation statements and survey rosters remain the original source's. The unchanged merged-row identity is asserted before reuse.

Aggregation sums non-canopy vertical levels within transects, then averages transects within zones, zones within sites, and sites within MPA/status/year. The source SOUTH-region filter, reference and MPA exclusions, and minimum-five-year site coverage are retained. Coverage counts the union of original fish site-years and source-filtered summer swath site-years, matching the combined source roster; it is not recalculated after zeroing small-fish counts. Density and raw biomass use the same 60 m2 denominator as the pipeline. Annual contributions are aligned to the existing harmonized raw keys, without adding or dropping keys.

Before any sensitivity output is accepted, the reconstructed baseline must reproduce retained annual raw values, retained annual adaptive-corrected log ratios, and every retained PISCO sheephead primary estimate **and SE**. Exact primary model routes, observed years, and supplied time are held fixed; Mean contrasts and 11-year Linear contrasts use `primary_uncertainty_utils.R` reconstruction and gradients. No model selection is repeated. Adaptive zero corrections are recalculated for each scenario using the original full MPA/taxon source histories, before subsetting to the identical retained contrast years.

The joint sensitivity uses `uncertainty_joint()`, updating only PISCO sheephead primary estimates/SEs and leaving all other primary entries unchanged. Its baseline must reproduce published joint estimates, SEs, and BH-adjusted P values. BH adjustment remains across the full nine-test family in each scenario, not just the two sheephead tests.

Raw shares are count-weighted physical-record summaries, separate from the nested primary aggregation. Shared sites can contribute to more than one MPA in the primary rules, but a raw record is counted only once in the retained-record share denominator. Raw biomass shares use the verified allometry and omit undefined per-record masses; they are not a substitute for the source's zero/missingness and nested-averaging conventions.

## Outputs and Tests

Current output filenames are `tables/table_s_sheephead_yoy_` followed by:

- `status.csv`: completion/skip state and inventory counts.
- `annual.csv`: baseline/filtered annual MPA/status values, removed amounts and fractions, by response.
- `raw_shares.csv`: study-campus and unique retained-record count/biomass shares, the 10 cm boundary count, and missing-size/count/zero inventories.
- `contrasts.csv`: fixed-route estimates/SEs, differences, observed years, and input times.
- `pooled.csv`: all nine joint tests per scenario, confidence intervals, BH P values, and detection-change flags.
- `inputs.csv`: fingerprints of consumed raw, source-code, primary, and comparison inputs.
- `protected.csv`: before/after fingerprints of primary harmonized files, cache files, primary statistics, and the joint comparison table.

Run `Rscript --vanilla code/tests/test_sheephead_yoy_check.R` for cutoff index/equality, roster, monotonicity, missingness, fixed-clock, optional-skip, and existing-output checks. The optional `--integration` argument repeats the full raw-input computation and all baseline/protection guards with output writing disabled. The ordinary Rscript module run exercises the actual entry point; no full `run_all` invocation is required.

Interpretation is conditional on the retained model routes, data years, original sampling weights and zero conventions. A 10 cm size threshold does not verify chronological age, maturity, or comparable age classification across programs. This diagnostic does not assess independent aged YOY cohorts or introduce new primary model claims.
