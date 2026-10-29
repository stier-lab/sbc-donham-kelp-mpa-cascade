# Step Direction Refit: Final Independent Numerical Audit

## Status And Phase Separation

The second full `run_all.R` completed at **12:22:26 on 2026-10-10**, exit 0, 4.18 minutes (`logs/pipeline_20261010_121815.log`). The parent subsequently completed the separate module-25 rerun, exit 0. This record audits the final **COUNT plus Step-direction** results against both the original `local_archive/lter_count_refit_2026-10-10/before` and frozen `count_only` phases. The COUNT-only results remain in `docs/LTER_COUNT_REFIT_2026-10-10.md`; they are not overwritten with the combined-refit numbers.

The only authorized executable change in module 08 is the label-matched `After - Before` Step contrast, documented in `docs/STEP_DIRECTION_FIX_2026-10-10.md`. Independent R parsing with `keep.source = FALSE` confirms that the entire module is identical to the original outside that helper. The analysis constants also have identical parsed expressions. Later clock clarifications changed comments, not clocks, formulas, inclusion rules or constants.

The final read-only checker reports **110 PASS, zero FAIL and one REVIEW** across 290 inspected files; the sole review is the pre-existing Naples routing discrepancy below. It validates input preservation, primary changes, Tables 1-3, nine-test BH adjustment, cohorts, module 25 and display metadata. All **33 synthetic audit tests** pass. The focused R tests independently confirm all six actual Step contrasts, both BA factor-level orders, character BA, reversed rows and preserved uncertainty. No production models, data, manuscript or supplement outputs were written by this audit. Only audit code, tests and these reports were updated.

## Primary Effects And Preserved Inputs

All **144 primary keys** remain: 108 CI and 36 pBACIPS; 68 Mean, 67 Linear, six Step and three Sigmoid. The original model-selection CSV covers 12 full-precision shared columns; the frozen COUNT-only primary CSV covers all **24 shared columns**. Relative to COUNT-only, exactly six `Mean` fields change and every other shared field is unchanged, including SE, SD, CI half-width, N, diagnostics, retained model and design.

| MPA / taxon / source / response | Original and COUNT-only mean | Corrected After - Before |
|---|---:|---:|
| Gull Island / red urchin / KFM / density | 0.430712882940 | -0.430712882940 |
| Naples / red urchin / LTER / density | 1.189470103710 | -1.189470103710 |
| Santa Barbara Island / purple urchin / KFM / density | 1.992445112450 | -1.992445112450 |
| Scorpion / sheephead / KFM / density | -1.543103140442 | 1.543103140442 |
| Scorpion / red urchin / KFM / density | -0.622103592406 | 0.622103592406 |
| Scorpion / red urchin / KFM / biomass | -1.253036353820 | 1.253036353820 |

Relative to the original, the only additional changed primary key is Campus Point sheephead/LTER/biomass from the COUNT correction: mean **-2.235358737183282 -> -2.283061430829335**, SE 1.633285823501082 -> 1.6349523994110413, CI half-width 3.406974526760784 -> 3.410450943191053, DW 2.0359 -> 2.0606 and DW P 0.4385 -> 0.4619. These final Campus values equal the frozen COUNT-only values; it remains Linear / pBACIPS / N = 22.

The Step correction adds **no raw/RR changes**. Independent byte checks against the original confirm all **6,190 non-target raw records** and **3,296 non-target RR records** are preserved, with headers, membership and order unchanged. The only input changes remain the authorized LTER sheephead biomass correction: 23 raw values and 19 annual log contrasts. Raw/RR exports retain 6,280/3,341 rows. Landsat RR, its annual index and site metadata remain whole-file byte-identical. Table 1 retains all 28 rows; only the LTER sheephead biomass summary changes, exactly as in phase 1.

Footprints remain **39 Mean/CI histories with before observations**, **59 matched-design-eligible histories**, **857 reference rows / 144 reference series**, and **144 atlas keys / 36 original figure pages**. No histories were admitted or dropped.

## Table 2 And FDR

All nine coefficient signs, k values, MPA counts and model labels are retained. BH q values were independently recomputed from all nine raw P values; exponentiated coefficients and interval bounds were independently checked. Both joint fits change because the corrected Step means change their shared variance components. Even lobster density, whose own primary means do not change, can therefore have a different pooled estimate and SE. Biomass tau2 is 0.2368 (original 0.2455; COUNT-only 0.2453); density tau2 is 1.358 (original/COUNT-only 0.839).

| Taxon / response | k / MPAs | Final exp(coefficient), 95% CI | Final raw P | Final BH q |
|---|---:|---|---:|---:|
| Giant kelp / biomass-index group | 33 / 19 | 2.230006 [1.401440, 3.548441] | 0.000961976 | 0.002885928 |
| Sheephead / biomass | 16 / 15 | 3.350078 [2.115960, 5.303985] | 0.000001473 | 0.000013258 |
| Red urchin / biomass | 14 / 9 | 2.109575 [1.345002, 3.308775] | 0.001472306 | 0.003312688 |
| Purple urchin / biomass | 10 / 9 | 0.692535 [0.428380, 1.119578] | 0.131741300 | 0.197611950 |
| Lobster / biomass | 5 / 5 | 1.419500 [0.767695, 2.624714] | 0.259743917 | 0.318256427 |
| Purple urchin / density | 13 / 11 | 0.136666 [0.043802, 0.426416] | 0.000877780 | 0.002885928 |
| Red urchin / density | 16 / 10 | 0.543057 [0.176012, 1.675520] | 0.282894602 | 0.318256427 |
| Lobster / density | 17 / 11 | 2.589584 [0.838168, 8.000711] | 0.096796508 | 0.174233714 |
| Sheephead / density | 20 / 15 | 0.726543 [0.235732, 2.239260] | 0.572539364 | 0.572539364 |

**Lobster density is the sole nominal/FDR support crossing.** Its original/COUNT-only exp(coefficient) was 2.955184 [1.394887, 6.260802], P = 0.005375385, q = 0.009675693. The final interval includes 1 and the final P/q exceed 0.05. It is now inconclusive, not evidence of no effect. Supported pooled contrasts decrease from five to four: **three positive biomass-group contrasts and one negative purple-urchin density contrast; no supported positive density contrast**.

These exponentiated pooled coefficients mix model-time-11 change contrasts with post-only level contrasts. They are not a universal inside/reference ratio, a common calendar-year-11 contrast, physical kelp biomass, a causal protection effect or a body-size effect. Most diver fits use sequential observation indices; the fixed endpoint does not impose a common elapsed-calendar duration.

## Table 3 Pathways

All six relationships and their k values are preserved. Five rows change relative to the original: four density pathways from the Step correction and the sheephead-biomass pathway from COUNT. Relative to COUNT-only, only the four density pathways change. All six slope intervals include zero and all six slope P values exceed 0.05; there is no newly supported primary pathway.

| Predictor -> response | k | Final slope, 95% CI | Final slope P |
|---|---:|---|---:|
| Purple-urchin density -> kelp index | 11 | -0.4860 [-1.0194, 0.0475] | 0.069379845 |
| Red-urchin density -> kelp index | 10 | -0.1764 [-0.4732, 0.1203] | 0.207592347 |
| Lobster density -> purple-urchin density | 11 | -0.2210 [-0.7998, 0.3577] | 0.410048560 |
| Sheephead density -> purple-urchin density | 11 | 0.5413 [-0.4950, 1.5776] | 0.267653145 |
| Lobster biomass -> purple-urchin biomass | 4 | -0.4215 [-0.8604, 0.0173] | 0.053853517 |
| Sheephead biomass -> purple-urchin biomass | 9 | 0.2884 [-0.1947, 0.7714] | 0.200927108 |

Three estimated slope directions change relative to the original. Neither those signs nor significant intercepts establish a cascade. The k = 4 biomass relationship remains unchanged and imprecise; failure to detect slopes is not evidence of equivalence or absence of trophic effects.

## Module 25 Final Checks

All **six dependent tables** are present and newer than the final primary CSV. Their published references match final Tables 2/3. No competing sensitivity models were run by this audit; checks use saved outputs and independent arithmetic.

- **Permutation:** six relationships, k = 11, 10, 11, 11, 4, 9; five Monte Carlo runs with 10,000 retained permutations and one exact run with 24 permutations. Only lobster biomass -> purple-urchin biomass reproduces the primary slope. The other five are distinct inverse-variance-grouped reconstructions, not the primary analyses.
- **AR1 inflation:** nine comparisons; median CI-width ratio **1.023**, sign flips **0/9**, nominal-support flips **0/9**. Independent DW-proxy arithmetic gives median rho **0.159125**, rounded 0.159, over 144 histories, with missing DW assigned rho = 0 by this helper. These summaries compare two separate legacy-grouped fits, not the primary joint fit.
- **DW:** 137 tested histories, **47 P < 0.05 (34.3%)**; Linear 24/67, Mean 19/64 and Step 4/6. Seven primary histories have no DW P in this table. Counts and medians independently match the final model-selection table and remain byte-identical to both earlier phases.
- **Offsets:** ten taxon/response summaries containing **2,801** annual contrasts. Independent reconstruction reproduces **2,801/3,140 = 89.2038%** of matched raw/RR records. Coverage is unchanged. The ten-row sensitivity table is byte-identical to COUNT-only; its original-to-final changes belong to COUNT, not the Step correction. These are annual-ratio shifts, not refitted effect-size or inference changes. The file's `rows_compared = 2801` is the retained summary count, not the 3,140-record fidelity denominator.
- **Point Conception:** one Landsat summary, 28 before / nine after observations; after mean lnRR **0.7094**, SE **0.1288**. Its saved main-fit reference refreshes from 0.7391 to **0.8020**, k = 33. This file does not refit the joint kelp model with Point Conception included; a caption must not call it an inclusion refit.

Two reporting pitfalls are load-bearing. Grouped permutation P values **0.0357** (purple-urchin density -> kelp) and **0.0469** (sheephead density -> purple-urchin density) newly cross 0.05, but their corresponding primary slopes have P = **0.0694** and **0.2677**. Do not promote these reconstructed tests to supported primary pathways. Similarly, the legacy AR1 table retains positive lobster density and negative purple-urchin biomass contrasts whose intervals exclude zero; these are not the final primary Table 2 results. Its zero printed P values are four-decimal rounding, not P = 0. The 0/9 flips apply only within its own paired legacy fits, not to before/after primary support.

## Unaffected Evidence And Metadata

All **28 archived disturbance, resilience and historical tables are byte-identical to the original**, including all seven historical-heatwave tables. Historical primary contrasts remain resistance n = 14, **0.422062896 [-0.575136647, 1.419262439]** and recovery n = 15, **0.505457067 [-0.387441892, 1.398356026]**, both q = **0.377169242**. The correction changes neither windows nor these results.

The map, legacy case-study, pooled-effects, recovery-curves, paired-resilience and primary-pathway displays have nonempty PNGs near 600 dpi and one-page PDFs. Map, legacy case-study and paired-resilience PNGs are byte-identical to the original. Pathway memberships and legend hashes match the final inputs. Module 32 retains **36 pBACIPS histories**, its nine-page companion and a one-page three-nonlinear-history comparison; all eight current code/input fingerprints match. An earlier read-only check caught stale hashes after comment updates; the current refreshed provenance and `--check-cache` validator now pass. These are metadata checks, not renewed rendered visual QA or verification of assembled manuscript figures.

No unexpected primary, raw/RR, cohort, clock or disturbance-table changes were found. The pre-existing Naples sheephead-biomass exclusion remains an explicit **VERIFY WITH EMILY** issue: the documented before-trend screen gives original P = 0.0667104504 and corrected P = 0.0988993259, neither below 0.05. Campus before-trend P remains 0.9420955674. This does not authorize automatically admitting Naples or altering routing.

## Recheck

```bash
python3 -B code/audit_lter_count_refit.py --step-direction-fixed --figures --shape-display --pretrend --json
python3 -B -m unittest discover -s code/tests -p 'test_lter_count_refit_audit.py'
Rscript --vanilla code/tests/test_step_effect_direction.R
Rscript code/R/32_pbacips_shape_display.R --check-cache
```

The CLI prints to stdout and never writes results. Snapshot fingerprints guard against changing inputs during inspection. Final SHA256: primary `89a20b19884bf00387ef9e6cf97d4ab76eefb50a09d45a33ac75153cca6afb6d`; Table 2 `4983678faf7f9641979fb2dd06043df83727204ac99bd3a9e7534925b729bf19`; Table 3 `09289d25b0ac5c45708a6b5a924ba02ac10e29fad9580533408227174878bb20`.
