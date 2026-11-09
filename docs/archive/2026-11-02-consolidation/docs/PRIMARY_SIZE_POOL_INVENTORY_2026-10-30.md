# Primary lobster and purple-urchin size inventory — October 30, 2026

Inventory of all 15 retained biomass histories for these two taxa: 204 paired years / 408 protected-or-reference year records. The inventory preserves sample membership. The subsequent missing-length correction is documented below. Seasonal reanalysis remains deferred.

## What matters biologically

The main concern is lobster size representation at Abalone Cove, followed by Matlahuayl and Point Vicente. Their primary resampled biomass histories include singleton measurements. Purple-urchin pools are generally larger, with isolated sparse years and demonstrable reuse of the same source measurements. Larger pools do not remove uncertainty from shared samples or habitat representativeness.

Five LTER lobsters counted with missing lengths formerly contributed zero biomass: one at Campus Point in 2017 (IVEE) and four at Naples in 2014 (NAPL). The source marks those lengths -99999. The corrected pipeline now estimates their mass from measured lobsters at the same site and survey date, retaining the original missing-length flags. Both histories enter the four-MPA lobster–purple-urchin biomass relationship. The focused treatment comparison leaves the pooled and pathway interpretations unchanged; see `LOBSTER_MISSING_LENGTHS_2026-10-30.md`.

## Which comparisons have small samples

Each denominator below is a positive-biomass protected or reference year, not an independent survey replicate. The <10 and <20 counts are descriptive screening thresholds, not exclusion rules or universal adequacy criteria. n counts animals represented by recorded frequencies, not independent dives. The direct LTER rows use measured animals in count surveys; they are not separate resampling pools. Zero-biomass years remain in the full inventory but do not require a size pool.

| Comparison | Method | Fewer than 10 / positive side-years | Singletons | Minimum n | Median n |
|---|---|---:|---:|---:|---:|
| Abalone Cove SMCA — lobster / PISCO | Size pool | 10/20 | 4 | 1 | 10 |
| Abalone Cove SMCA — purple urchin / PISCO | Size pool | 2/26 | 0 | 2 | 66 |
| Anacapa Island SMR 2003 — purple urchin / KFM | Size pool | 0/33 | 0 | 11 | 497 |
| Campus Point SMCA — lobster / LTER | Direct lengths | 0/24 | 0 | 15 | 174 |
| Campus Point SMCA — purple urchin / PISCO | Size pool | 2/26 | 0 | 4 | 64.5 |
| Gull Island SMR — purple urchin / PISCO | Size pool | 1/30 | 0 | 4 | 136.5 |
| Harris Point SMR — purple urchin / PISCO | Size pool | 0/30 | 0 | 16 | 92.5 |
| Matlahuayl SMR — lobster / PISCO | Size pool | 4/14 | 1 | 1 | 14 |
| Naples SMCA — lobster / LTER | Direct lengths | 1/24 | 0 | 6 | 58.5 |
| Naples SMCA — purple urchin / PISCO | Size pool | 0/26 | 0 | 37 | 142.5 |
| Point Vicente SMCA — lobster / PISCO | Size pool | 3/23 | 1 | 1 | 30 |
| Point Vicente SMCA — purple urchin / PISCO | Size pool | 1/27 | 0 | 8 | 161 |
| Scorpion SMR — purple urchin / PISCO | Size pool | 0/36 | 0 | 79 | 168.5 |
| South Point SMR — purple urchin / KFM | Size pool | 0/34 | 0 | 186 | 673.5 |
| South Point SMR — purple urchin / PISCO | Size pool | 0/28 | 0 | 168 | 416 |

Singleton lobster pools occur at Abalone Cove reference in 2013 and protected reef in 2016, 2017 and 2023; Matlahuayl protected reef in 2023; and Point Vicente protected reef in 2012. A corrected resampling algorithm now preserves those lengths, but cannot estimate unsampled size variation. No retained purple-urchin pool is a singleton.

Purple-urchin pools below10 occur at Abalone Cove protected reef in2022 (n=2) and2023 (4); Campus Point reference in2017 (9) and protected reef in2023 (4); Gull Island protected reef in2005 (4); and Point Vicente protected reef in2016 (8). Exact source sites and all years are in the inventory CSV.

## Shared measurements

- South Point KFM and PISCO purple-urchin histories reuse size records in all14 common paired years, on both sides (28side-years). Every PISCO size record in those pools is also used by KFM; KFM additionally includes smaller urchins. Thus these are not identical biomass estimates, but their size-derived errors are dependent. Shared animals represent28–88% of the KFM pools.
- Campus Point and Naples PISCO purple-urchin reference pools reuse129 measured animals from Arroyo Quemado West in2020 and107 in2021. These constitute roughly65% of the Campus Point reference pools and41–43% of the Naples pools.
- These are exact same-file/row matches across the 15 retained histories. This does not prove the absence of physical overlap among differently recorded samples, shared count-survey errors or common coefficient error. Within each comparison, the same MPA/status/year pool can also supply multiple count samples; repeated use does not add independent size information.

## Recommended response

1. Completed: source checks confirm antennae-only sightings rather than recoverable lengths. The pipeline now uses explicitly flagged same-site/date mean individual mass. The inventory validates those corrected annual values; raw missing lengths remain visible.
2. Prioritize Abalone Cove and Point Vicente for any subsequent biomass uncertainty analysis: both enter the four-MPA lobster–urchin relationship. Matlahuayl also has sparse pools and affects the pooled lobster result, but does not enter that particular pathway. Preserve current comparison membership while separating size uncertainty from changes in sampled sites/years.
3. Carry one realization of each shared size sample to all of its consumers, using record-level overlap for partly overlapping pools. Include independent count-survey units where recoverable. Singleton and constant-size pools need explicit limitations or justified external information; resampling them cannot create a meaningful confidence interval for population size variation.
4. Do not remove all years below an arbitrary n threshold. A targeted influence check can be informative, but would also change temporal coverage. Current evidence identifies where uncertainty arises; it does not establish the direction or magnitude of its effect on pooled or trophic conclusions.
5. Keep the present manuscript disclosure. Decide on numerical propagation after the missing-length treatment and sampling-unit structure are resolved. Do not add an extensive new manuscript section simply to narrate this inventory.

## Reproduction, validation and scope

Run `Rscript code/R/55_primary_size_inventory.R`, then `python3 code/tests/check_primary_size_inventory.py`. DONHAM_DATA_REPO can identify the source workspace. This is a standalone inventory module, not a new default pipeline stage. Source inputs and code fingerprints are in `tables/table_s_primary_size_inputs.csv`.

Outputs: `table_s_primary_size_inventory.csv` (all408side-years), `comparisons.csv` (15histories), `sharing.csv` (exact cross-history source reuse), and `direct_validation.csv` (48LTER side-year reconstructions), all with the `table_s_primary_size_` prefix. Restricted row-level membership and missing-length details remain in git-ignored `local_archive/primary_size_inventory_2026-10-30/`; they must not be published as part of this audit. The test requires those local audit ledgers and source files.

Checks passed: exact retained keys/years and both sides; every positive resampled biomass cell has a pool; source row frequency counts independently checked; summary thresholds independently tallied; exact overlap independently reconstructed; all47 overlapping pool mean masses match the prior lobster reconciliation; all48 direct LTER annual biomass values reproduce to <1e-10; all five missing lengths are source sentinel values; protected primary inputs/results retain their hashes. No coefficient-error draws, uncertainty refits, seasonal restrictions, source downloads or emails were performed.

The inventory traces the primary size-pool lookups, rather than estimating whether surveyed sizes represent all counted animals. It does not establish observer error, independent dive counts, effective sample size, allometric covariance or total confidence-interval coverage. Those require further evidence.
