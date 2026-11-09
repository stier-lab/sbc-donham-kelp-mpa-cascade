# Five missing LTER lobster lengths: source check and treatment comparison

## Finding

These are antennae-only sightings, not recoverable recorded lengths in the checked files. All five have COUNT=1, NUM_AO=1 and SIZE_MM=-99999. Official metadata defines NUM_AO as antenna-only observations and -99999 as not recorded/not available. Four are at Naples (NAPL), August26,2014; one is at Campus Point (IVEE), August22,2017. A sixth comparable record at Carpinteria is outside the two retained primary lobster histories.

The source supplied with the paper and its archived GitHub-folder copy are byte-identical. The newer official SBC77.13 release, dated June18,2026, retains all five missing lengths with the same dates, transects, replicates and antenna-only flags. Its lab-NAS copy has MD5 f99d82e37f961b87a8c8457733b6e867, matching the official metadata; 460159bytes. The attempted EDI download returned403, but the verified NAS copy resolved access. No credential was needed for the completed source check. No individual field sheets or unpublished measurements were available in the inspected project files, and this is not a claim that those cannot exist elsewhere.

Source: Reed and Miller, SBC LTER reef lobster abundance, size and effort, package knb-lter-sbc.77.13; DOI https://doi.org/10.6073/pasta/ecb6064d4483b250d6c2897d71caee00. The analysis still uses its supplied2023 source; this comparison did not replace the entire release or extend survey years.

## Primary treatment (applied October 30, 2026)

Keep each animal in the count and estimate its mass from measured lobsters at the same site and survey date, explicitly flagging the estimated contribution. Use mean individual mass after the nonlinear length-to-mass conversion, not the mass calculated from mean length. The available donor samples contain605 measured lobsters at IVEE on the same2017 date and159 at NAPL on the same2014 date. Their mean masses are380.36g and485.74g, respectively.

This requires the assumption that body-hidden animals have a comparable size distribution to visible animals at that site/date; the source cannot establish it. It is a transparent imputation, not recovery of true lengths. Retain alternative assignments as a focused check. For a subsequent full biomass-uncertainty analysis, retain imputation uncertainty and donor sampling dependence rather than treating estimated values as measured lengths. More repetitions or a zero conditional interval cannot resolve nonrandom visibility.

Assigning zero mass to an animal known to be present is unsuitable as the preferred estimate. Dropping whole affected replicate surveys also discards valid measurements and changes survey effort. It is useful as a comparison, not the preferred primary treatment. The same-site/date mean-mass treatment is now applied in source processing. Original lengths and counts remain intact; estimated contributions are explicitly flagged. The scoped promotion changes only two annual biomass values and their response ratios.

## What changes numerically

All treatments preserve the145 primary histories, the same paired years, existing reference assignments, retained contrast models/clocks, and the existing Figure4 biomass model structure. Only the two affected LTER annual biomass values are changed before regenerating response ratios, including the original adaptive correction. The replicate-exclusion treatment preserves years/sites but changes within-year survey effort; it is therefore not a pure size-assignment comparison.

| Treatment of unknown mass | Campus Point protected biomass2017 | Naples protected biomass2014 |
|---|---:|---:|
| Same-site/date mean individual mass | +0.165% | +2.516% |
| Same-transect/date mean individual mass | +0.172% | +2.225% |
| Smallest measured lobster at site/date | +0.014% | +1.064% |
| Largest measured lobster at site/date | +1.289% | +6.224% |
| Exclude affected replicate surveys | -0.649% | -9.407% |

Percentages are relative to the previous zero-contribution calculation. Observed length ranges are32–150mm at IVEE and60–110mm at NAPL. Assigning those endpoints is an observed-range sensitivity, not a guaranteed biological bound on unobserved lengths.

Across all five alternatives:

- All nine pooled response directions and all four FDR-supported classifications are unchanged.
- Figure4 species slope directions and nominal significance classifications are unchanged; every refitted fixed structure converged without singularity.
- The four-MPA lobster–purple-urchin biomass slope remains negative: -0.555412 with zero contributions and -0.554902 with the current site/date mean-mass treatment. Its nominal P changes from 0.037651 to 0.037699; exact 24-permutation P remains 0.0833 under every treatment. The small sample and shared-reference limitations remain.
- Campus Point’s retained contrast changes from 1.505652 to 1.505595, and Naples’ from0.848073 to0.850114 under site/date mean mass. These are log-scale contrast estimates, not percentage biomass increases.

This resolves the effect of these five records under the tested treatments. It does not establish that all biomass uncertainty, sparse PISCO lobster size pools or shared purple-urchin size samples are negligible. No new seasonal analysis or full joint uncertainty propagation was run. The standalone treatment comparison refits the primary Figure4 structure. The subsequent full pipeline rebuild also regenerates its existing AR1 uncertainty checks.

## Reproduction and checks

`Rscript code/R/56_lobster_missing_lengths.R` writes `tables/table_s_lobster_missing_{annual,contrasts,pooled,pathway,temporal,fit_checks,inputs}.csv`. This is a standalone sensitivity module; it does not modify source data, caches, primary results or manuscripts. `python3 code/tests/check_lobster_missing_lengths.py` independently checks source flags, same-date donor membership, hierarchical annual percentage changes, conclusion comparisons and fingerprints.

The primary-treatment reconstruction reproduces all48 LTER annual side values to<1e-10, both retained contrast estimates/SEs, the current nine pooled table entries and lobster–urchin pathway coefficient/P. Exact temporal-cache cohort, years and inputs match current files. Independent checks pass; git diff --check passes. Protected input/result hashes are unchanged. Original-row comparisons and donor details are private in ignored `local_archive/lobster_missing_lengths_2026-10-30/`.

Implementation: source `code/R/lobster_missing_length_utils.R` applies the flagged treatment before spatial joins. `code/promote_lobster_missing_lengths.py` preserves unrelated rows and records the release lineage. The five retained observations affect only Campus Point 2017 and Naples 2014; the sixth source record at Carpinteria does not enter primary comparisons. Reproduce the standalone inventory and treatment comparison after a primary rebuild with modules55 and56, then run their independent checks. These modules require the local source workspace and do not overwrite primary results.

After a scoped source promotion and primary-contrast rebuild, explicitly refresh the checksum-bound reference registry with `Rscript code/R/export_reference_registry.R /path/to/sbc-donham-kelp-mpa-data` before the design/geography modules. The registry checks must remain active; changed numerical inputs require a new fingerprint even when site membership is unchanged. The October30 analysis run stopped at this guard, then continued from its verified corrected snapshot after the registry refresh.
