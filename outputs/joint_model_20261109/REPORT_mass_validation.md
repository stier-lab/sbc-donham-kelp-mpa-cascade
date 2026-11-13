# Is the density, mass and biomass model sound? Validation report

Rules were fixed in `MASS_VALIDATION_PLAN.md` before the runs. Review only; nothing is adopted. Effects are log changes at ten years; intervals are 90%.


## Scorecard

| Check | Result |
|---|---|
| M1 convergence | FAIL |
| M2 seeds agree | PASS |
| M3 calibration | PASS |
| M4 robust to bimodal | PASS |
| M4 robust to visitnoise | FAIL |
| M4 robust to coupled | PASS |
| M5 time placebo | PASS |
| M6 space placebo (mass, biomass) | FAIL |
| M7 predictive check of mass | FAIL |
| M8 hold-out prediction | FAIL |
| M9 leave one MPA out | PASS |
| M10 independent arithmetic | PASS |
| M11 input integrity | PASS |

Known failures of the density part from 2026-11-10 that this battery does not re-test: sheephead time placebo; urchin count layer misfit; calibration narrowly below its rule; effect size depends on the curve prior.

## M1. Convergence

| Taxon | Random-walk prior, seed 1 | seed 2 | Independent-increment prior |
|---|---|---|---|
| Lobster | yes (1.002 / 2696 / 0) | yes (1.008 / 619 / 0) | NO (1.001 / 3124 / 1) |
| Sheephead | yes (1.002 / 1512 / 0) | yes (1.002 / 1859 / 0) | yes (1.002 / 2210 / 0) |
| Purple urchin | yes (1.004 / 1341 / 0) | yes (1.003 / 1601 / 0) | yes (1.002 / 1626 / 0) |
| Red urchin | yes (1.004 / 1448 / 0) | yes (1.002 / 1518 / 0) | yes (1.005 / 1171 / 0) |
| Kelp | yes (1.002 / 1474 / 0) | yes (1.004 / 1360 / 0) | yes (1.004 / 1197 / 0) |

(R-hat / effective draws / divergent transitions.)

**FAIL.**

## M2. Reproducibility

Largest difference between two independent fits: Lobster 0.003; Sheephead 0.002; Purple urchin 0.030; Red urchin 0.008; Kelp 0.018.

**PASS.**

## M3. Calibration when the model is right

Simulated on each taxon's real survey and measurement design, then refitted (2 chains, 1,500 + 1,500).

| Scenario | Taxon | Data sets | R-hat < 1.05 | Mass: truth | mean error | covered | excludes zero | Biomass: truth | mean error | covered | excludes zero |
|---|---|---|---|---|---|---|---|---|---|---|---|
| std | Lobster | 20 | 20 | +0.20 | -0.06 | 18/20 | 5/20 | +2.12 | -0.16 | 17/20 | 20/20 |
| std | Sheephead | 20 | 19 | +0.15 | -0.06 | 18/20 | 4/20 | +0.38 | +0.01 | 19/20 | 8/20 |
| std | Purple urchin | 20 | 6 | -0.10 | +0.01 | 18/20 | 4/20 | -0.75 | +0.10 | 18/20 | 7/20 |
| std | Red urchin | 20 | 9 | +0.12 | -0.01 | 18/20 | 6/20 | -0.13 | -0.01 | 18/20 | 2/20 |
| std | Kelp | 20 | 17 | +0.00 | -0.04 | 20/20 | 0/20 | +0.70 | -0.21 | 15/20 | 6/20 |
| null | Lobster | 20 | 19 | +0.00 | -0.01 | 19/20 | 1/20 | +0.00 | +0.03 | 18/20 | 2/20 |
| null | Sheephead | 20 | 18 | +0.00 | -0.03 | 19/20 | 1/20 | +0.00 | +0.03 | 19/20 | 1/20 |
| null | Purple urchin | 20 | 8 | +0.00 | -0.00 | 19/20 | 1/20 | +0.00 | -0.02 | 18/20 | 2/20 |
| null | Red urchin | 20 | 9 | +0.00 | +0.01 | 18/20 | 2/20 | +0.00 | +0.02 | 18/20 | 2/20 |
| null | Kelp | 20 | 16 | +0.00 | -0.04 | 20/20 | 0/20 | +0.00 | -0.08 | 14/20 | 6/20 |
| trend | Lobster | 20 | 19 | +0.00 | +0.01 | 20/20 | 0/20 | +1.92 | -0.09 | 17/20 | 20/20 |
| trend | Sheephead | 20 | 19 | +0.00 | -0.00 | 18/20 | 2/20 | +0.22 | +0.06 | 18/20 | 4/20 |
| trend | Purple urchin | 20 | 7 | +0.00 | -0.04 | 18/20 | 2/20 | -0.65 | +0.05 | 18/20 | 6/20 |
| trend | Red urchin | 20 | 11 | +0.00 | -0.00 | 20/20 | 0/20 | -0.25 | +0.00 | 18/20 | 2/20 |
| trend | Kelp | 20 | 16 | +0.00 | +0.01 | 20/20 | 0/20 | +0.70 | -0.17 | 16/20 | 6/20 |

- std, mass: coverage 0.92 of 100 (rule ≥ 0.84)
- std, biomass: coverage 0.87 of 100 (rule ≥ 0.84)
- null, mass: coverage 0.95 of 100 (rule ≥ 0.84)
- null, biomass: coverage 0.87 of 100 (rule ≥ 0.84)
- null: 5% of mass intervals exclude zero (rule ≤ 16%)
- trend, mass: coverage 0.96 of 100 (rule ≥ 0.84)
- trend, biomass: coverage 0.87 of 100 (rule ≥ 0.84)
- trend: 4% of mass intervals exclude zero (rule ≤ 16%)

**PASS.**

## M4. Robustness when an assumption is wrong

`bimodal`: two size cohorts instead of a Gamma. `visitnoise`: extra visit-to-visit noise in mean size. `coupled`: mean size falls with station density while the model assumes independence.

| Scenario | Taxon | Data sets | R-hat < 1.05 | Mass: truth | mean error | covered | excludes zero | Biomass: truth | mean error | covered | excludes zero |
|---|---|---|---|---|---|---|---|---|---|---|---|
| bimodal | Lobster | 20 | 19 | +0.20 | -0.04 | 19/20 | 3/20 | +2.12 | -0.15 | 18/20 | 20/20 |
| bimodal | Sheephead | 20 | 20 | +0.15 | -0.03 | 19/20 | 6/20 | +0.38 | +0.03 | 19/20 | 8/20 |
| bimodal | Purple urchin | 20 | 4 | -0.10 | +0.01 | 19/20 | 4/20 | -0.75 | +0.10 | 18/20 | 7/20 |
| bimodal | Red urchin | 20 | 7 | +0.12 | -0.02 | 18/20 | 5/20 | -0.13 | -0.02 | 18/20 | 2/20 |
| bimodal | Kelp | 20 | 10 | +0.00 | -0.03 | 20/20 | 0/20 | +0.70 | -0.21 | 15/20 | 6/20 |
| visitnoise | Lobster | 20 | 19 | +0.20 | -0.08 | 19/20 | 3/20 | +2.12 | -0.18 | 16/20 | 20/20 |
| visitnoise | Sheephead | 20 | 20 | +0.15 | -0.06 | 17/20 | 3/20 | +0.38 | +0.01 | 19/20 | 8/20 |
| visitnoise | Purple urchin | 20 | 18 | -0.10 | +0.02 | 18/20 | 5/20 | -0.75 | +0.12 | 17/20 | 8/20 |
| visitnoise | Red urchin | 20 | 11 | +0.12 | +0.01 | 18/20 | 5/20 | -0.13 | +0.01 | 17/20 | 2/20 |
| visitnoise | Kelp | 20 | 20 | +0.00 | -0.04 | 19/20 | 1/20 | +0.70 | -0.22 | 14/20 | 6/20 |
| coupled | Lobster | 20 | 20 | +0.20 | -0.09 | 18/20 | 3/20 | +2.12 | -0.19 | 17/20 | 20/20 |
| coupled | Sheephead | 20 | 20 | +0.15 | -0.06 | 19/20 | 2/20 | +0.38 | -0.00 | 19/20 | 7/20 |
| coupled | Purple urchin | 20 | 1 | -0.10 | +0.03 | 19/20 | 2/20 | -0.75 | +0.11 | 19/20 | 5/20 |
| coupled | Red urchin | 20 | 8 | +0.12 | -0.01 | 18/20 | 4/20 | -0.13 | -0.02 | 18/20 | 2/20 |
| coupled | Kelp | 20 | 10 | +0.00 | -0.03 | 18/20 | 2/20 | +0.70 | -0.20 | 15/20 | 6/20 |

`bimodal` — mass: coverage 0.95 of 100 (rule ≥ 0.84); biomass: coverage 0.88 of 100 (rule ≥ 0.84). **PASS.**

`visitnoise` — mass: coverage 0.91 of 100 (rule ≥ 0.84); biomass: coverage 0.83 of 100 (rule ≥ 0.84). **FAIL.**

`coupled` — mass: coverage 0.92 of 100 (rule ≥ 0.84); biomass: coverage 0.88 of 100 (rule ≥ 0.84). **PASS.**

## M5. Placebo in time

KFM surveys before the real designation with a fake designation in 1992. A sound model finds nothing. Kelp uses 1999: the KFM kelp counts in the model start in 1996, so 1992 leaves no before years and the planned runs could not be built (changed 2026-11-12, before any kelp placebo result was seen).

| Taxon | Prior | Density | Mean mass | Biomass | Converged |
|---|---|---|---|---|---|
| Purple urchin | rw | -1.09 [-2.67, +0.62] | -0.66 [-1.36, +0.04] | -1.75 [-3.54, +0.08] | NO (1.016 / 350 / 20) |
| Purple urchin | pw | -0.46 [-1.80, +0.58] | -0.23 [-0.75, +0.20] | -0.69 [-2.06, +0.46] | NO (1.019 / 356 / 9) |
| Red urchin | rw | -0.41 [-1.48, +0.71] | +0.05 [-0.68, +0.79] | -0.36 [-1.65, +1.03] | NO (1.012 / 312 / 3) |
| Red urchin | pw | -0.15 [-0.87, +0.46] | +0.01 [-0.31, +0.35] | -0.14 [-0.91, +0.58] | NO (1.013 / 476 / 10) |
| Kelp | rw | -0.14 [-4.94, +4.46] | -0.01 [-1.30, +1.26] | -0.15 [-5.02, +4.76] | NO (1.012 / 405 / 122) |
| Kelp | pw | -0.01 [-1.83, +1.84] | -0.00 [-0.48, +0.48] | -0.01 [-1.86, +1.91] | NO (1.010 / 343 / 45) |

**PASS.**

## M6. Placebo in space

Protected stations dropped; one reference station per MPA relabelled as protected; ten random relabellings.

| Taxon | Runs | Density: mean, excluding zero | Mean mass: mean, excluding zero | Biomass: mean, excluding zero |
|---|---|---|---|---|
| Lobster | 10 | +0.41, 0/10 ✗ | +0.04, 0/10 | +0.46, 0/10 ✗ |
| Sheephead | 10 | +0.31, 0/10 ✗ | +0.06, 0/10 | +0.37, 0/10 |
| Purple urchin | 10 | -0.01, 0/10 | +0.01, 0/10 | -0.00, 0/10 |
| Red urchin | 10 | -0.14, 0/10 | +0.04, 0/10 | -0.10, 0/10 |
| Kelp | 10 | -0.67, 0/10 ✗ | -0.20, 1/10 | -0.87, 0/10 ✗ |

✗ marks a failure of the rule. The rule is applied to mass and biomass; the density column repeats yesterday's density placebo under today's fits.

**FAIL.**

## M7. Predictive check of the mass part

Mean sizes replicated from the fitted model, compared with the data: spread and tails of group means by program, and two protected-versus-reference contrasts.

| Taxon | Statistics outside their 90% replicate interval | Contrast statistics |
|---|---|---|
| Lobster | 5: sd_log|LTER, q05_log|LTER, q95_log|LTER, sd_log|PISCO, q05_log|PISCO | after_minus_before: observed +0.76, replicates [+0.13, +0.84]; age10plus_minus_age0to5: observed +0.34, replicates [+0.28, +0.47] |
| Sheephead | 5: q05_log|KFM, sd_log|LTER, sd_log|PISCO, q05_log|PISCO, q95_log|PISCO | after_minus_before: observed +0.52, replicates [+0.34, +0.68]; age10plus_minus_age0to5: observed +0.54, replicates [+0.37, +0.58] |
| Purple urchin | 5: sd_log|LTER, q05_log|LTER, sd_log|PISCO, q05_log|PISCO, q95_log|PISCO | after_minus_before: observed +0.52, replicates [+0.42, +0.58]; age10plus_minus_age0to5: observed -0.04, replicates [-0.11, +0.07] |
| Red urchin | 6: sd_log|LTER, q05_log|LTER, q95_log|LTER, sd_log|PISCO, q05_log|PISCO, q95_log|PISCO | after_minus_before: observed +0.94, replicates [+0.82, +1.02]; age10plus_minus_age0to5: observed +0.32, replicates [+0.23, +0.45] |
| Kelp | 3: sd_log|PISCO, q05_log|PISCO, q95_log|PISCO | after_minus_before: observed -0.06, replicates [-0.22, +0.02]; age10plus_minus_age0to5: observed -0.20, replicates [-0.24, -0.10] |

24 of 49 statistics outside (rule: at most 10%, and no contrast statistic outside).

**FAIL.**

## M8. Predicting size data the model was not shown

| Taxon | Split | Held-out groups | 90% intervals cover | 50% intervals cover | By program |
|---|---|---|---|---|---|
| Lobster | 1 | 211 | 0.92 | 0.65 | LTER 0.94 (n=79); PISCO 0.92 (n=132) |
| Lobster | 2 | 215 | 0.95 | 0.63 | LTER 0.99 (n=80); PISCO 0.93 (n=135) |
| Sheephead | 1 | 1,142 | 0.92 | 0.58 | KFM 0.90 (n=69); LTER 0.95 (n=59); PISCO 0.92 (n=1014) |
| Sheephead | 2 | 994 | 0.91 | 0.58 | KFM 0.93 (n=85); LTER 0.89 (n=36); PISCO 0.91 (n=873) |
| Purple urchin | 1 | 402 | 0.91 | 0.58 | KFM 0.91 (n=88); LTER 0.89 (n=217); PISCO 0.95 (n=97) |
| Purple urchin | 2 | 267 | 0.96 | 0.61 | KFM 0.95 (n=94); LTER 0.97 (n=74); PISCO 0.96 (n=99) |
| Red urchin | 1 | 264 | 0.91 | 0.57 | KFM 0.91 (n=88); LTER 0.89 (n=78); PISCO 0.92 (n=98) |
| Red urchin | 2 | 287 | 0.90 | 0.57 | KFM 0.89 (n=90); LTER 0.87 (n=84); PISCO 0.92 (n=113) |
| Kelp | 1 | 583 | 0.91 | 0.49 | KFM 0.87 (n=45); PISCO 0.91 (n=538) |
| Kelp | 2 | 570 | 0.91 | 0.48 | KFM 0.92 (n=50); PISCO 0.91 (n=520) |

**FAIL.**

## M9. Leave one MPA out

| Taxon | Runs | Largest shift: density | mean mass | biomass | Sign reversals of clear main effects | Shifts above 0.25 |
|---|---|---|---|---|---|---|
| Lobster | 13 of 13 | -0.21 (Santa Barbara Island SMR) | +0.07 (Matlahuayl SMR) | +0.22 (Matlahuayl SMR) | none | none |
| Sheephead | 17 of 17 | -0.19 (Point Vicente SMCA) | -0.07 (Santa Barbara Island SMR) | -0.20 (Point Vicente SMCA) | none | none |
| Purple urchin | 13 of 13 | -0.22 (Naples SMCA) | -0.07 (Scorpion SMR) | -0.24 (Scorpion SMR) | none | none |
| Red urchin | 13 of 13 | -0.15 (Point Vicente SMCA) | +0.06 (Harris Point SMR) | +0.17 (Anacapa Island SMR 2003) | none | none |
| Kelp | 13 of 13 | +0.37 (Harris Point SMR) | -0.07 (Campus Point SMCA) | +0.40 (Harris Point SMR) | none | density: Harris Point SMR +0.37; biomass: Harris Point SMR +0.40 |

**PASS.**

## M10. Independent arithmetic

`test_multi.py` on the real inputs:

```
2. kelp<-pur: 1574 of 1674 target station-years linked, 103 stations; 2150 sampled links checked
2. kelp<-red: 1574 of 1674 target station-years linked, 103 stations; 2150 sampled links checked
ALL CHECKS PASSED
3. lob  rw        gamma     mass layer log density: model -8656685.592805  independent -8656685.592805  (1621 groups)
3. lob  piecewise gamma     mass layer log density: model -8460272.906103  independent -8460272.906103  (1621 groups)
3. lob  rw        lognormal mass layer log density: model -144616.880614  independent -144616.880614  (1621 groups)
3. pur  rw        gamma     mass layer log density: model -5694027.188889  independent -5694027.188889  (1941 groups)
3. pur  piecewise gamma     mass layer log density: model -4084472.820527  independent -4084472.820527  (1941 groups)
3. pur  rw        lognormal mass layer log density: model -250993.040969  independent -250993.040969  (1941 groups)
3. shp  rw        gamma     mass layer log density: model -28009735.215175  independent -28009735.215175  (6571 groups)
3. shp  piecewise gamma     mass layer log density: model -14404145.229315  independent -14404145.229315  (6571 groups)
3. shp  rw        lognormal mass layer log density: model -742714.090899  independent -742714.090899  (6571 groups)
3. mass layer matches the independent calculation; worst relative difference 7.537057491258459e-15
ALL CHECKS PASSED (1-4)
```

**PASS.**

## M11. Input integrity

4,035 pools, 8,260 same-record units and 305 KFM sheephead visits reproduce the integrated review's values from source rows.

**PASS.**

