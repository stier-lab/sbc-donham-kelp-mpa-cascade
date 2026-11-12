# Landsat kelp canopy as its own response

Rules and limits on inference were fixed in `CANOPY_PLAN.md` before the full runs. Review only; nothing is adopted. Effects are log changes in spring canopy density ten years after designation (18 frozen comparisons, 48 polygons, 1984-2021); intervals are 90%.

What this response can be used to say: surface canopy biomass in spring over whole reserves and their reference areas. Not plant numbers, size structure or autumn kelp; weaker below about 10 m (Reed, Rassweiler, Nelson and Kui 2026).

Input check: the quarterly table reproduces Tom Bell's delivered annual-maximum table for 58 non-empty polygons with largest difference 0.0.

## Scorecard

| Check | Result |
|---|---|
| C1 independent arithmetic | PASS |
| C2 convergence of the main fits | PASS |
| C6 seeds agree | PASS |
| C3 calibration | INCOMPLETE |
| C4 shuffled labels: bias (R1a) | INCOMPLETE |
| C4 shuffled labels: interval calibration (R1b) | INCOMPLETE |
| C5 time placebo | PASS |
| C7 leave one MPA out | INCOMPLETE |

## Main result: spring canopy (January-June quarters), all depths, flat baseline

| Curve prior | Effect at ten years | Effect at four years | Effect at eight years | Quarters without canopy | Converged |
|---|---|---|---|---|---|
| random walk | +0.49 [-0.22, +1.15] | +0.18 [-0.16, +0.50] | +0.36 [-0.23, +0.91] | 17% | R-hat 1.002, 0 div., ESS 1102 |
| independent increments | +0.20 [-0.25, +0.80] | +0.11 [-0.25, +0.53] | +0.09 [-0.38, +0.56] | 17% | R-hat 1.003, 0 div., ESS 997 |

MPAs designated in 2012 have nine years of canopy data after designation; the ten-year effect leans on the 2003 MPAs and on the curve.

## C8. Sensitivities, reported side by side (none chosen by its result)

| Variant | Random walk | Independent increments |
|---|---|---|
| Spring, designation year after (main) | +0.49 [-0.22, +1.15] | +0.20 [-0.25, +0.80] |
| Annual maximum (the earlier satellite input) | +0.36 [-0.33, +1.08] | +0.10 [-0.34, +0.60] |
| July-December quarters | +0.32 [-0.31, +0.98] | +0.12 [-0.32, +0.68] |
| All four quarters | +0.35 [-0.28, +1.01] | +0.11 [-0.33, +0.64] |
| Spring, designation year counted as before | +0.50 [-0.16, +1.20] | +0.20 [-0.27, +0.78] |
| Spring, pixels to 10 m deep | +0.36 [-0.24, +1.02] | +0.14 [-0.29, +0.67] |
| Spring, pixels deeper than 10 m | +0.51 [-0.20, +1.23] | +0.22 [-0.31, +0.84] |
| Spring, with pre-trend | +0.86 [+0.06, +1.71] | +0.36 [-0.27, +1.09] |

Estimated pre-trend of protected against reference polygons, per decade: random walk -0.17 [-0.42, +0.07]; independent increments -0.09 [-0.33, +0.14].

## C6. Seeds

Seed 1 +0.49 [-0.22, +1.15]; seed 2 +0.50 [-0.18, +1.16]. **PASS.**

## C3. Simulated data on the real design, noise at real-data levels

| Scenario | Data sets (usable) | True effect | Mean error | Rule | Result |
|---|---|---|---|---|---|
| effect present | 20 (20) | +0.55 | -0.10 | at least 16 of 20 intervals cover the truth | 16 of 20 cover |
| no effect | 9 (9) | +0.00 | +0.05 | at most 4 of 20 intervals exclude zero | 8 of 9 cover |

**INCOMPLETE.**

## C4. Shuffled labels

The protected label is given at random to one polygon of each comparison (the MPA or one of its reference areas).

| Prior | Shuffles (usable) | Mean of shuffled effects | SD | Mean ÷ SE | 90% intervals covering zero | Real effect | Shuffles at least as extreme | Randomization p |
|---|---|---|---|---|---|---|---|---|
| random walk | 26 of 40 (26) | -0.05 | 0.21 | -1.3 | 100% ✗ | +0.49 [-0.22, +1.15] | 0 of 26 | 0.037 |
| independent increments | not finished | | | | | +0.20 [-0.25, +0.80] | | |

✗ marks a failed rule (mean within 2 standard errors of zero; 80-97% of intervals cover zero). **Bias: INCOMPLETE. Interval calibration: INCOMPLETE.**

## C5. Time placebo (years before designation only, fake designation year)

| Fake year | Random walk | Independent increments |
|---|---|---|
| 1990 | -0.26 [-0.68, +0.17] | -0.15 [-0.54, +0.18] |
| 1993 | -0.22 [-0.66, +0.25] | -0.13 [-0.53, +0.21] |
| 1996 | -0.17 [-0.68, +0.31] | -0.11 [-0.51, +0.25] |
| 1999 | -0.15 [-0.73, +0.47] | -0.09 [-0.54, +0.34] |

Intervals covering zero: random walk 4 of 4, independent increments 4 of 4 (rule: at least 3 of 4 each). **PASS.**

## C7. Leave one MPA out (random walk)

| MPA left out | Effect at ten years | Shift from the main fit |
|---|---|---|
| Abalone Cove SMCA | +0.47 [-0.24, +1.16] | -0.02 |
| Cabrillo SMR | +0.56 [-0.11, +1.27] | +0.07 |
| Campus Point SMCA | +0.52 [-0.19, +1.23] | +0.02 |
| Carrington Pt SMR | +0.48 [-0.18, +1.26] | -0.02 |
| Farnsworth Onshore SMCA | +0.50 [-0.20, +1.18] | +0.00 |
| Gull Island SMR | +0.49 [-0.23, +1.21] | -0.00 |
| Harris Point SMR | +0.53 [-0.19, +1.22] | +0.04 |
| Judith Rk SMR | +0.48 [-0.25, +1.21] | -0.01 |
| Naples SMCA | +0.48 [-0.24, +1.15] | -0.01 |
| Point Conception SMR | +0.60 [-0.10, +1.30] | +0.10 |
| Point Dume SMCA | +0.53 [-0.22, +1.21] | +0.03 |
| Point Dume SMR | +0.51 [-0.18, +1.22] | +0.02 |

12 of 18 done; largest shift 0.10. **INCOMPLETE.**

