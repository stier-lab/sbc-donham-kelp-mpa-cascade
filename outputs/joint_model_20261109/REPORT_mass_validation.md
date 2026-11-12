# Is the density, mass and biomass model sound? Validation report

Rules were fixed in `MASS_VALIDATION_PLAN.md` before the runs. Review only; nothing is adopted. Effects are log changes at ten years; intervals are 90%.


## Scorecard

| Check | Result |
|---|---|
| M1 convergence | INCOMPLETE |
| M2 seeds agree | INCOMPLETE |
| M3 calibration | INCOMPLETE |
| M4 robust to bimodal | INCOMPLETE |
| M4 robust to visitnoise | INCOMPLETE |
| M4 robust to coupled | INCOMPLETE |
| M5 time placebo | INCOMPLETE |
| M6 space placebo (mass, biomass) | INCOMPLETE |
| M7 predictive check of mass | INCOMPLETE |
| M8 hold-out prediction | INCOMPLETE |
| M9 leave one MPA out | INCOMPLETE |
| M10 independent arithmetic | PASS |
| M11 input integrity | PASS |

Known failures of the density part from 2026-11-10 that this battery does not re-test: sheephead time placebo; urchin count layer misfit; calibration narrowly below its rule; effect size depends on the curve prior.

## M1. Convergence

| Taxon | Random-walk prior, seed 1 | seed 2 | Independent-increment prior |
|---|---|---|---|
| Lobster | yes (1.002 / 2696 / 0) | yes (1.008 / 619 / 0) | NO (1.001 / 3124 / 1) |
| Sheephead | not finished | not finished | not finished |
| Purple urchin | yes (1.004 / 1341 / 0) | yes (1.003 / 1601 / 0) | yes (1.002 / 1626 / 0) |
| Red urchin | yes (1.004 / 1448 / 0) | yes (1.002 / 1518 / 0) | yes (1.005 / 1171 / 0) |
| Kelp | yes (1.002 / 1474 / 0) | yes (1.004 / 1360 / 0) | yes (1.004 / 1197 / 0) |

(R-hat / effective draws / divergent transitions.)

**INCOMPLETE.**

## M2. Reproducibility

Largest difference between two independent fits: Lobster 0.003; Purple urchin 0.030; Red urchin 0.008; Kelp 0.018.

**INCOMPLETE.**

## M3. Calibration when the model is right

Simulated on each taxon's real survey and measurement design, then refitted (2 chains, 1,500 + 1,500).

| Scenario | Taxon | Data sets | R-hat < 1.05 | Mass: truth | mean error | covered | excludes zero | Biomass: truth | mean error | covered | excludes zero |
|---|---|---|---|---|---|---|---|---|---|---|---|
| std | Lobster | 20 | 20 | +0.20 | -0.06 | 18/20 | 5/20 | +2.12 | -0.16 | 17/20 | 20/20 |
| std | Purple urchin | 6 | 3 | -0.10 | -0.04 | 6/6 | 2/6 | -0.75 | +0.10 | 5/6 | 2/6 |

- std, mass: coverage 0.92 of 26 (rule ≥ 0.78)
- std, biomass: coverage 0.85 of 26 (rule ≥ 0.78)

**INCOMPLETE.**

## M4. Robustness when an assumption is wrong

`bimodal`: two size cohorts instead of a Gamma. `visitnoise`: extra visit-to-visit noise in mean size. `coupled`: mean size falls with station density while the model assumes independence.

| Scenario | Taxon | Data sets | R-hat < 1.05 | Mass: truth | mean error | covered | excludes zero | Biomass: truth | mean error | covered | excludes zero |
|---|---|---|---|---|---|---|---|---|---|---|---|

`bimodal` — . **INCOMPLETE.**

`visitnoise` — . **INCOMPLETE.**

`coupled` — . **INCOMPLETE.**

## M5. Placebo in time

KFM surveys before the real designation with a fake designation in 1992. A sound model finds nothing.

| Taxon | Prior | Density | Mean mass | Biomass | Converged |
|---|---|---|---|---|---|
| Purple urchin | rw | -1.09 [-2.67, +0.62] | -0.66 [-1.36, +0.04] | -1.75 [-3.54, +0.08] | NO (1.016 / 350 / 20) |
| Purple urchin | pw | -0.46 [-1.80, +0.58] | -0.23 [-0.75, +0.20] | -0.69 [-2.06, +0.46] | NO (1.019 / 356 / 9) |
| Red urchin | rw | -0.41 [-1.48, +0.71] | +0.05 [-0.68, +0.79] | -0.36 [-1.65, +1.03] | NO (1.012 / 312 / 3) |
| Red urchin | pw | -0.15 [-0.87, +0.46] | +0.01 [-0.31, +0.35] | -0.14 [-0.91, +0.58] | NO (1.013 / 476 / 10) |
| Kelp | rw | not finished | | | |
| Kelp | pw | not finished | | | |

**INCOMPLETE.**

## M6. Placebo in space

Protected stations dropped; one reference station per MPA relabelled as protected; ten random relabellings.

| Taxon | Runs | Density: mean, excluding zero | Mean mass: mean, excluding zero | Biomass: mean, excluding zero |
|---|---|---|---|---|
| Lobster | 7 | +0.46, 0/7 ✗ | +0.06, 0/7 | +0.52, 0/7 ✗ |
| Sheephead | 0 | | | |
| Purple urchin | 0 | | | |
| Red urchin | 0 | | | |
| Kelp | 0 | | | |

✗ marks a failure of the rule. The rule is applied to mass and biomass; the density column repeats yesterday's density placebo under today's fits.

**INCOMPLETE.**

## M7. Predictive check of the mass part

Mean sizes replicated from the fitted model, compared with the data: spread and tails of group means by program, and two protected-versus-reference contrasts.

| Taxon | Statistics outside their 90% replicate interval | Contrast statistics |
|---|---|---|
| Lobster | not finished | |
| Sheephead | not finished | |
| Purple urchin | not finished | |
| Red urchin | not finished | |
| Kelp | not finished | |

0 of 0 statistics outside (rule: at most 10%, and no contrast statistic outside).

**INCOMPLETE.**

## M8. Predicting size data the model was not shown

| Taxon | Split | Held-out groups | 90% intervals cover | 50% intervals cover | By program |
|---|---|---|---|---|---|
| Lobster | 1 | not finished | | | |
| Lobster | 2 | not finished | | | |
| Sheephead | 1 | not finished | | | |
| Sheephead | 2 | not finished | | | |
| Purple urchin | 1 | not finished | | | |
| Purple urchin | 2 | not finished | | | |
| Red urchin | 1 | not finished | | | |
| Red urchin | 2 | not finished | | | |
| Kelp | 1 | not finished | | | |
| Kelp | 2 | not finished | | | |

**INCOMPLETE.**

## M9. Leave one MPA out

| Taxon | Runs | Largest shift: density | mean mass | biomass | Sign reversals of clear main effects | Shifts above 0.25 |
|---|---|---|---|---|---|---|
| Lobster | 0 | | | | | |
| Sheephead | 0 | | | | | |
| Purple urchin | 0 | | | | | |
| Red urchin | 0 | | | | | |
| Kelp | 0 | | | | | |

**INCOMPLETE.**

## M10. Independent arithmetic

`test_multi.py` on the real inputs:

```
2. pur<-shp: 1678 of 1768 target station-years linked, 102 stations; 1622 sampled links checked
2. red<-lob: 1567 of 1768 target station-years linked, 101 stations; 1466 sampled links checked
2. red<-shp: 1678 of 1768 target station-years linked, 102 stations; 1622 sampled links checked
2. kelp<-pur: 1574 of 1674 target station-years linked, 103 stations; 2150 sampled links checked
2. kelp<-red: 1574 of 1674 target station-years linked, 103 stations; 2150 sampled links checked
ALL CHECKS PASSED
3. lob  rw        mass layer log density: model -8656685.592805  independent -8656685.592805  (1621 groups)
3. lob  piecewise mass layer log density: model -8460272.906103  independent -8460272.906103  (1621 groups)
3. pur  rw        mass layer log density: model -5694027.188889  independent -5694027.188889  (1941 groups)
3. pur  piecewise mass layer log density: model -4084472.820527  independent -4084472.820527  (1941 groups)
3. shp  rw        mass layer log density: model -28009735.215175  independent -28009735.215175  (6571 groups)
3. shp  piecewise mass layer log density: model -14404145.229315  independent -14404145.229315  (6571 groups)
3. mass layer matches the independent calculation; worst relative difference 3.3062228408178182e-15
ALL CHECKS PASSED (1-3)
```

**PASS.**

## M11. Input integrity

4,035 pools, 8,260 same-record units and 305 KFM sheephead visits reproduce the integrated review's values from source rows.

**PASS.**

