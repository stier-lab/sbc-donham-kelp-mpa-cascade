# LTER COUNT Refit: Independent Numerical Audit

## Status And Scope

`run_all.R` completed at 12:04:31 on 2026-10-10, exit 0, 3.8 minutes (`logs/pipeline_20261010_120043.log`). **The COUNT-only numerical audit is closed and preserved separately from the subsequently discovered Step-direction correction.** This comparison uses the original preserved snapshot in `local_archive/lter_count_refit_2026-10-10/before`, not the preliminary module-08 outputs. Immediately before changing the Step helper, all **298 files** in the parent's `local_archive/lter_count_refit_2026-10-10/count_only` snapshot were independently checked against the audited live results; all were byte-identical. The parent subsequently confirmed that the first module-25 run had completed before this freeze. Its six frozen tables have now been independently checked; module 32 was not final-checked for this phase. No production data, models, manuscript or supplement files were changed by this audit. The subsequently authorized helper fix and combined-refit results are recorded in `docs/STEP_DIRECTION_FIX_2026-10-10.md` and `docs/STEP_DIRECTION_REFIT_2026-10-10.md`, not folded into these numbers.

The preserved COUNT-only primary CSV has SHA256 `7612d8d3b5ad47b4e529093908c73f80cfe75541b6b8c1b3c02ad5839192dc49`; Table 2 `c3a72bd0e35d05e0a8b11ee96c59d7f72b1dfa52a114be95447ae304e77a440d`; Table 3 `bd66e62c285f996753063e59d803decc53d223f76db32db8becaf1c515b6240e`.

The original post-`run_all` audit recorded **93 passing checks, zero failures and two review notices** (then-pending module 25; the pre-existing Naples routing caveat below). That dated result is retained. The expanded reusable checker, `code/audit_lter_count_refit.py`, now has **33 passing synthetic tests**. A later replay of the frozen COUNT-only results, including their module-25 tables, gives **99 passing checks, zero failures and one review notice** (Naples). The earlier `code/tests/test_response_summary.R` and upstream `code/tests/test_lter_sheephead_count.R` passes are unchanged; the upstream test covered the actual module-06 COUNT assignment, reconstruction, non-target preservation, immutable-baseline and double-correction tests, writing only temporary fixtures that it removed. No competing model refits were run.

## Input Preservation

- Raw export: 6,280 records before and after. All **6,190 non-target records are byte-identical**, with header, order and membership unchanged. Of 90 target LTER/sheephead/Bio records, 23 values change: Campus reference 1, Naples MPA 18 and Naples reference 4. No Campus MPA value changes.
- Response-ratio export: 3,341 records before and after. All **3,296 non-target records are byte-identical**, with header, order and membership unchanged. All 45 target records have changes to normalized values or ratios, but only **19 annual log contrasts** change: Campus 2020 alone and 18 Naples years. The distinction matters because within-series normalization changes `mpa` and `reference` even when their log ratio is unchanged.
- Campus 2020 reference biomass changes from 2.95515610 to 3.75887078 g/m2; its log contrast changes from -2.4382246 to -2.6767381 (difference -0.2385135).
- Landsat response ratios, the Landsat annual index and site metadata are whole-file byte-identical to the preserved baseline. These checks establish preservation of the harmonized exports, not an independent reprocessing of every original monitoring record.

## Primary Effects And Cohorts

All **144 primary keys** are unchanged. Only `Campus Point SMCA | S. pulcher | LTER | Bio` changes numerically in the shared full-precision primary columns:

| Quantity | Original | Corrected |
|---|---:|---:|
| Model-time-11 log effect | -2.235358737183282 | -2.283061430829335 |
| SE | 1.633285823501082 | 1.6349523994110413 |
| CI half-width, not an endpoint | 3.406974526760784 | 3.410450943191053 |
| Durbin-Watson statistic | 2.0359 | 2.0606 |
| Durbin-Watson P | 0.4385 | 0.4619 |

Campus remains **Linear / pBACIPS / N = 22**; Naples is still absent from the primary set. The effect-size module is byte-identical to its preserved source, so no routing or fitting rule changed in this refit. Membership and counts remain: 108 CI and 36 pBACIPS; 68 Mean, 67 Linear, 6 Step and 3 Sigmoid; **39 Mean/CI histories have before observations**. The atlas retains 144 keys across 36 original figure pages with identical panel addresses and annual footprints. Matched-design eligibility remains 59 histories. The refreshed registry contains 857 reference rows and all 144 primary series.

The archive lacks a standalone `sumstats_final.csv`, so the primary comparison uses all 12 full-precision columns in its `tables/table_model_selection.csv`. SD and the additional primary metadata columns are not independently covered by that numerical baseline comparison; atlas and design-coverage checks separately cover the annual footprints and classifications.

## Main Tables

Table 1 retains all 28 taxon/source/response rows. Only the LTER sheephead biomass summary changes: N = 90; mean 8.14155973468287 -> 9.7339501584682; SD 16.34285109254513 -> 19.7992864608379; SE 1.7226877637804616 -> 2.0870280420794; CI half-width 3.422943892531885 -> 4.14688026488406. These are descriptive pooled annual/status observations, not 90 independent MPA effects or a body-size effect.

Table 2 retains all nine tests. All five biomass coefficients shift slightly because they share the joint multilevel fit; changing one source effect does not imply that only the sheephead pooled coefficient can change. The four density fits **and their FDR q values are unchanged**. BH adjustment was independently recomputed over all nine raw P values for both snapshots. No coefficient changes direction or crosses the nominal or FDR 0.05 threshold.

| Biomass response | Original log coefficient | Corrected log coefficient | Corrected exp(coefficient), 95% CI | Corrected q |
|---|---:|---:|---|---:|
| Sheephead | 1.161057621 | 1.160842121 | 3.192620716 [2.029408683, 5.022560078] | 0.00002275240 |
| Lobster | 0.280381180 | 0.280534496 | 1.323837209 [0.719374088, 2.436208066] | 0.407606769 |
| Purple urchin | -0.389174294 | -0.389353530 | 0.677494713 [0.421611835, 1.088676951] | 0.159217493 |
| Red urchin | 0.651890038 | 0.651775391 | 1.918944684 [1.231406038, 2.990361088] | 0.009675693 |
| Giant kelp | 0.739255198 | 0.739116089 | 2.094083712 [1.324645672, 3.310460061] | 0.005805652 |

Biomass tau2 changes from 0.2455 to 0.2453 and QE from 575.2805 to 575.3786. The purple-urchin exponentiated percentage changes rounding from -32.2% to -32.3%. These exponentiated pooled coefficients combine model-time-11 change and post-only level estimands; they must not be described as a single common inside/reference biomass ratio, physical kelp biomass, causal protection effect or body-size response. The fixed endpoint is an observation-index endpoint for most diver histories, not necessarily eleven elapsed calendar years.

Table 3 retains all six pathways and their k values. Only sheephead biomass -> purple-urchin biomass changes: slope 0.2896 -> 0.2884; SE 0.2050 -> 0.2043; CI [-0.1951, 0.7743] -> [-0.1947, 0.7714]; P 0.2005924864 -> 0.2009271079. This remains an imprecise cross-MPA association, not evidence of no pathway or proof of a causal cascade. The other five pathway rows are numerically unchanged.

## Before-Trend Routing Caveat

The audit evaluates only the existing parsed `add_time_columns()` and `test_linear_before()` definitions, without sourcing the effect-size module. **Both baseline and current helpers sort by year**, then use sequential observation indices, not elapsed calendar time. Campus has 11 before/designation and 11 after observations; its before-slope P is **0.9420955674 in both snapshots**. Naples has 12 before/designation and 11 after observations; its P is **0.0667104504 -> 0.0988993259**.

The hardcoded statement that Naples biomass was linear in the before does not reproduce the documented P < 0.05 screen in either snapshot. This discrepancy **predates the COUNT correction**. The same-rule refit appropriately leaves membership fixed; it must not be represented as a new automatic selection decision, and this audit does not authorize admitting Naples. The earlier unsorted-row diagnostic giving Campus P = 0.671 was not the production helper and is superseded by the values above.

## Unaffected Evidence And Displays

All **28 archived disturbance, resilience and historical-heatwave tables are byte-identical**, including modern resistance/recovery, ecological memory, compound disturbance and the seven historical-heatwave tables. The historical primary contrasts remain resistance n = 14, 0.422062896 [-0.575136647, 1.419262439] and recovery n = 15, 0.505457067 [-0.387441892, 1.398356026], both q = 0.377169242. Failure to detect a difference is not an equivalence result.

The current canonical map, legacy case-study display, pooled-effects figure, recovery-curves figure, paired kelp-resilience figure and primary-pathway display each have a nonempty PNG at approximately 600 dpi and a one-page PDF. All 36 separate atlas PDFs have one page and the combined original atlas has 36 pages. The pathway membership and legend hashes match the refreshed primary and Table 3 inputs. PNGs for the map, legacy case-study display and paired kelp resilience are byte-identical to the original snapshot; pooled-effects, recovery and pathway PNGs change. Different PDF bytes alone are not evidence of a visual change. **These are metadata checks, not renewed visual QA; the new module-32 Figure 2 and companion display still need their final metadata check.**

## COUNT-Only Checker

The checker below tests the COUNT-only expectation (only Campus changes in the primary set). **It must not be treated as a final combined-refit gate after the separately authorized six-Step sign correction.** A combined refit must be compared to both preserved phases, allowing exactly the six listed Step means plus the Campus change relative to the original, and keeping SE/CI/N and membership fixed for the Step rows.

```bash
python3 -B code/audit_lter_count_refit.py --results-snapshot local_archive/lter_count_refit_2026-10-10/count_only --figures --pretrend --json
python3 -B -m unittest discover -s code/tests -p 'test_lter_count_refit_audit.py' -v
```

The CLI prints its report to stdout and does not write files. `--pending-25` explicitly excludes the in-progress prereview tables and is not a final sensitivity pass. All inspected files must remain byte-stable throughout the audit. This record retains the COUNT-only numbers; later combined-refit results belong in the separate Step-direction record.

## Frozen Module 25 Completion

All six published-slope/P references match the frozen COUNT-only Table 3; the six pathways retain five 10,000-permutation reconstructions and one exact 24-permutation test. Five grouped reconstructions do not reproduce the primary regression slopes, so their permutation P values are not replacements for Table 3. The nine separate grouped AR1 comparisons have median CI-width ratio **1.008**, with no sign or nominal-support flips. The DW distribution remains **47/137** tested histories, including six Step histories. The ten offset summaries contain **2,801** annual contrasts; independent adaptive-offset arithmetic reproduces **2,801/3,140 = 89.2038%** of matched raw/RR records. DW distribution and offset coverage are byte-identical to the original baseline. Copy timestamps do not establish the frozen runner's completion; that provenance is the parent's completion report plus the verified pre-Step snapshot. These are COUNT-only sensitivities, not the later combined-refit values.
