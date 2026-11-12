# Validation round 2: shuffled labels, evidence split, pre-trend

Rules were fixed in `VALIDATION_PLAN_v2.md` before the runs. Review only; nothing is adopted. Effects are log changes at ten years; intervals are 90%.

## Scorecard

| Check | Result |
|---|---|
| R1a bias | FAIL |
| R1b calibration | FAIL |
| R2 evidence split | INCOMPLETE |
| R3a sheephead time placebo with pre-trend | PASS |
| R3b flat vs pre-trend | INCOMPLETE |
| R3 simulation (pre-trend option recovers the effect) | INCOMPLETE |

## R1. Shuffled-label placebo

Within each MPA and its reference stations the protected label is dealt at random to the same number of stations; all stations and surveys are kept. If the model is unbiased the shuffled effects centre on zero. Density only.

| Taxon | Prior | Shuffles (usable) | Mean of shuffled effects | SD | Mean ÷ SE | 90% intervals covering zero | Real effect | Shuffles at least as extreme | Randomization p |
|---|---|---|---|---|---|---|---|---|---|
| Lobster | random walk | 40 (40) | -0.10 ✗ | 0.27 | -2.5 | 98% ✗ | +1.75 [+1.04, +2.43] | 0 of 40 | 0.024 |
| Lobster | independent increments | 40 (40) | -0.07 ✗ | 0.17 | -2.5 | 100% ✗ | +1.82 [+1.10, +2.53] | 0 of 40 | 0.024 |
| Sheephead | random walk | 40 (40) | +0.14 ✗ | 0.14 | +6.2 | 82% | +0.20 [-0.09, +0.50] | 13 of 40 | 0.341 |
| Sheephead | independent increments | 6 of 40 (6) | +0.03 | 0.06 | +1.4 | 100% ✗ | +0.15 [-0.08, +0.41] | 0 of 6 | 0.143 |
| Purple urchin | random walk | not finished | | | | | -0.91 [-1.80, +0.02] | | |
| Purple urchin | independent increments | not finished | | | | | -0.28 [-0.97, +0.37] | | |
| Red urchin | random walk | not finished | | | | | -0.63 [-1.23, -0.05] | | |
| Red urchin | independent increments | not finished | | | | | -0.27 [-0.76, +0.20] | | |
| Kelp | random walk | not finished | | | | | +0.96 [+0.00, +1.91] | | |
| Kelp | independent increments | not finished | | | | | +0.60 [-0.23, +1.53] | | |

✗ marks a failed rule: R1a, the mean of the shuffled effects is within 2 standard errors of zero; R1b, between 80% and 97% of the shuffled 90% intervals cover zero. The randomization p is reported, not gated.

**R1a: FAIL** (fails: lob_rw, lob_pw, shp_rw). **R1b: FAIL** (fails: lob_rw, lob_pw). 3 of 10 taxon-prior sets complete.

## R2. Effects by kind of evidence

Each taxon fitted on all MPAs (main), on MPAs with before-and-after surveys only, and on post-only MPAs only. Counts are animals counted at protected stations before designation, at the before-and-after MPAs.

| Taxon | Prior | MPAs used | Density | Mean mass | Biomass | Converged (R-hat / divergences) |
|---|---|---|---|---|---|---|
| Lobster | random walk | all (13) | +1.75 [+1.04, +2.43] | +0.38 [+0.15, +0.64] | +2.13 [+1.43, +2.88] | yes (1.002 / 0) |
| Lobster | random walk | before-and-after (6) | +1.72 [+0.94, +2.51] | +0.40 [+0.03, +0.75] | +2.12 [+1.26, +2.99] | yes (1.003 / 1) |
| Lobster | random walk | post-only (7) | +1.11 [+0.04, +2.06] | +0.26 [-0.11, +0.63] | +1.37 [+0.26, +2.41] | yes (1.005 / 1) |
| Lobster | independent increments | all (13) | +1.82 [+1.10, +2.53] | +0.23 [-0.03, +0.49] | +2.05 [+1.26, +2.81] | yes (1.001 / 1) |
| Lobster | independent increments | before-and-after (6) | +1.82 [+0.93, +2.68] | +0.21 [-0.08, +0.56] | +2.03 [+1.07, +2.95] | yes (1.001 / 1) |
| Lobster | independent increments | post-only (7) | +1.26 [+0.15, +2.37] | +0.07 [-0.14, +0.35] | +1.33 [+0.14, +2.41] | yes (1.004 / 3) |
| Sheephead | random walk | all (17) | +0.20 [-0.09, +0.50] | +0.33 [+0.10, +0.56] | +0.53 [+0.14, +0.89] | yes (1.002 / 0) |
| Sheephead | random walk | before-and-after (7) | +0.15 [-0.29, +0.55] | +0.24 [-0.07, +0.51] | +0.39 [-0.13, +0.90] | yes (1.002 / 0) |
| Sheephead | random walk | post-only (10) | +0.19 [-0.27, +0.66] | +0.35 [+0.03, +0.66] | +0.54 [-0.03, +1.10] | yes (1.002 / 0) |
| Sheephead | independent increments | all (17) | +0.15 [-0.08, +0.41] | +0.26 [-0.00, +0.47] | +0.41 [+0.05, +0.76] | yes (1.002 / 0) |
| Sheephead | independent increments | before-and-after (7) | +0.09 [-0.16, +0.42] | +0.14 [-0.07, +0.40] | +0.23 [-0.13, +0.61] | yes (1.002 / 0) |
| Sheephead | independent increments | post-only (10) | +0.10 [-0.18, +0.44] | +0.25 [-0.05, +0.54] | +0.35 [-0.10, +0.77] | yes (1.006 / 0) |
| Purple urchin | random walk | all (13) | -0.91 [-1.80, +0.02] | +0.25 [-0.05, +0.53] | -0.66 [-1.66, +0.25] | yes (1.004 / 0) |
| Purple urchin | random walk | before-and-after (7) | -0.41 [-1.50, +0.63] | +0.24 [-0.11, +0.60] | -0.17 [-1.25, +0.99] | yes (1.004 / 0) |
| Purple urchin | random walk | post-only (6) | -1.51 [-2.93, -0.13] | +0.20 [-0.29, +0.64] | -1.31 [-2.75, +0.21] | yes (1.001 / 0) |
| Purple urchin | independent increments | all (13) | -0.28 [-0.97, +0.37] | +0.09 [-0.10, +0.33] | -0.19 [-0.92, +0.50] | yes (1.002 / 0) |
| Purple urchin | independent increments | before-and-after (7) | -0.07 [-0.61, +0.46] | +0.11 [-0.12, +0.39] | +0.03 [-0.58, +0.61] | yes (1.003 / 0) |
| Purple urchin | independent increments | post-only (6) | -0.60 [-1.82, +0.53] | +0.06 [-0.21, +0.35] | -0.54 [-1.79, +0.65] | yes (1.002 / 0) |
| Red urchin | random walk | all (13) | -0.63 [-1.23, -0.05] | +0.41 [+0.08, +0.73] | -0.22 [-0.90, +0.44] | yes (1.004 / 0) |
| Red urchin | random walk | before-and-after (7) | -0.16 [-0.72, +0.43] | +0.34 [-0.05, +0.72] | +0.18 [-0.49, +0.89] | yes (1.004 / 0) |
| Red urchin | random walk | post-only (6) | -1.13 [-2.22, -0.12] | +0.35 [-0.13, +0.83] | -0.79 [-1.91, +0.38] | yes (1.004 / 0) |
| Red urchin | independent increments | all (13) | -0.27 [-0.76, +0.20] | +0.21 [-0.07, +0.52] | -0.06 [-0.66, +0.52] | yes (1.005 / 0) |
| Red urchin | independent increments | before-and-after (7) | -0.03 [-0.48, +0.41] | +0.11 [-0.13, +0.40] | +0.08 [-0.43, +0.62] | yes (1.004 / 0) |
| Red urchin | independent increments | post-only (6) | -0.32 [-1.12, +0.39] | +0.15 [-0.15, +0.52] | -0.17 [-1.05, +0.66] | yes (1.004 / 0) |
| Kelp | random walk | all (13) | +0.96 [+0.00, +1.91] | +0.01 [-0.20, +0.22] | +0.96 [+0.01, +1.95] | yes (1.002 / 0) |
| Kelp | random walk | before-and-after (7) | not finished | | | |
| Kelp | random walk | post-only (6) | +1.10 [+0.11, +2.14] | -0.07 [-0.41, +0.27] | +1.03 [-0.03, +2.13] | yes (1.004 / 0) |
| Kelp | independent increments | all (13) | +0.60 [-0.23, +1.53] | -0.03 [-0.23, +0.14] | +0.57 [-0.30, +1.52] | yes (1.004 / 0) |
| Kelp | independent increments | before-and-after (7) | not finished | | | |
| Kelp | independent increments | post-only (6) | not finished | | | |

| Taxon | Before-and-after MPAs | Counted before (protected stations) | Largest MPA's share | Of those MPAs, with fewer than 20 counted before | Label under the proposed rule |
|---|---|---|---|---|---|
| Lobster | 6 of 13 | 48 | 83% | 5 | rw: rests on the after-designation contrast; pw: rests on the after-designation contrast |
| Sheephead | 7 of 17 | 884 | 33% | 1 | rw: rests on the after-designation contrast; pw: rests on the after-designation contrast |
| Purple urchin | 7 of 13 | 118,944 | 32% | 0 | rw: rests on the after-designation contrast; pw: rests on the after-designation contrast |
| Red urchin | 7 of 13 | 21,819 | 35% | 0 | rw: rests on the after-designation contrast; pw: rests on the after-designation contrast |
| Kelp | 7 of 13 | 28,594 | 64% | 3 | rw: not finished; pw: not finished |

Proposed label rule (Adrian to confirm): "before-and-after supported" needs the before-and-after-only density interval to exclude zero and at least three of those MPAs with 20 or more animals counted before designation at protected stations.

## R3. Pre-trend sensitivity

`--pretrend` lets protected stations at before-and-after MPAs trend against their references before designation and extends that line through the after years. The flat baseline stays the main specification.

### Time placebos (KFM years before the real designation, fake designation year)

| Taxon | Fake year | Flat baseline | With pre-trend | Pre-trend per decade | Converged flat / pre-trend |
|---|---|---|---|---|---|
| Sheephead | 1990 | -1.51 [-2.94, -0.21] | -0.75 [-2.84, +1.23] | -0.66 [-2.00, +0.72] | yes (1.002 / 0) / yes (1.002 / 0) |
| Sheephead | 1993 | -2.43 [-4.44, -0.53] | -1.85 [-4.35, +0.37] | -0.48 [-1.62, +0.67] | yes (1.003 / 0) / yes (1.002 / 0) |
| Sheephead | 1996 | -3.21 [-5.82, -0.68] | -2.57 [-5.28, +0.36] | -0.48 [-1.40, +0.52] | yes (1.002 / 0) / yes (1.002 / 0) |
| Sheephead | 1999 | -3.82 [-7.85, +0.07] | -3.09 [-7.00, +0.96] | -0.74 [-1.61, +0.05] | yes (1.002 / 0) / yes (1.002 / 0) |
| Purple urchin | 1990 | -0.97 [-2.39, +0.44] | -0.85 [-2.88, +1.11] | -0.11 [-1.27, +1.15] | yes (1.003 / 0) / yes (1.004 / 0) |
| Purple urchin | 1993 | -1.27 [-3.06, +0.38] | -0.97 [-3.27, +1.15] | -0.21 [-1.26, +0.90] | yes (1.002 / 0) / yes (1.004 / 0) |
| Purple urchin | 1996 | -1.21 [-3.65, +1.20] | -0.59 [-3.18, +2.22] | -0.41 [-1.40, +0.52] | yes (1.001 / 0) / yes (1.005 / 0) |
| Purple urchin | 1999 | -1.92 [-5.47, +1.62] | -1.34 [-4.96, +2.32] | -0.42 [-1.31, +0.42] | yes (1.002 / 0) / yes (1.007 / 0) |
| Red urchin | 1990 | -0.35 [-1.32, +0.60] | +0.62 [-0.93, +2.00] | -0.73 [-1.72, +0.20] | yes (1.009 / 0) / yes (1.006 / 0) |
| Red urchin | 1993 | -0.42 [-1.65, +0.78] | +0.42 [-1.20, +1.98] | -0.57 [-1.37, +0.24] | yes (1.005 / 0) / yes (1.005 / 0) |
| Red urchin | 1996 | -0.76 [-2.29, +0.92] | -0.03 [-1.96, +1.94] | -0.41 [-1.09, +0.33] | yes (1.004 / 0) / yes (1.007 / 0) |
| Red urchin | 1999 | -1.20 [-3.82, +1.42] | -0.45 [-3.43, +2.25] | -0.37 [-1.01, +0.22] | yes (1.004 / 0) / yes (1.007 / 0) |

Sheephead: flat-baseline interval covers zero at 1 of 4 fake years; with pre-trend at 4 of 4 (rule R3a: at least 3 of 4). **PASS.**

### Real data: flat baseline against pre-trend

| Taxon | Prior | Density, flat | Density, pre-trend | Biomass, flat | Biomass, pre-trend | Pre-trend per decade | Depends on the baseline? | Converged |
|---|---|---|---|---|---|---|---|---|
| Lobster | random walk | +1.75 [+1.04, +2.43] | +1.26 [+0.49, +1.99] | +2.13 [+1.43, +2.88] | +1.64 [+0.87, +2.43] | +0.70 [+0.10, +1.31] | no | yes (1.002 / 0) |
| Lobster | independent increments | +1.82 [+1.10, +2.53] | +1.29 [+0.53, +2.08] | +2.05 [+1.26, +2.81] | +1.52 [+0.70, +2.36] | +0.71 [+0.09, +1.31] | no | yes (1.002 / 0) |
| Sheephead | random walk | +0.20 [-0.09, +0.50] | not finished | +0.53 [+0.14, +0.89] | | | | |
| Sheephead | independent increments | +0.15 [-0.08, +0.41] | not finished | +0.41 [+0.05, +0.76] | | | | |
| Purple urchin | random walk | -0.91 [-1.80, +0.02] | -0.96 [-2.02, +0.03] | -0.66 [-1.66, +0.25] | -0.71 [-1.79, +0.35] | +0.07 [-0.68, +0.83] | no | yes (1.006 / 0) |
| Purple urchin | independent increments | -0.28 [-0.97, +0.37] | -0.23 [-1.01, +0.39] | -0.19 [-0.92, +0.50] | -0.14 [-0.93, +0.54] | -0.18 [-0.89, +0.55] | no | yes (1.006 / 0) |
| Red urchin | random walk | -0.63 [-1.23, -0.05] | -0.80 [-1.51, -0.15] | -0.22 [-0.90, +0.44] | -0.39 [-1.13, +0.37] | +0.25 [-0.26, +0.79] | no | yes (1.003 / 0) |
| Red urchin | independent increments | -0.27 [-0.76, +0.20] | -0.32 [-0.94, +0.19] | -0.06 [-0.66, +0.52] | -0.11 [-0.79, +0.54] | +0.08 [-0.43, +0.58] | no | yes (1.016 / 0) |
| Kelp | random walk | +0.96 [+0.00, +1.91] | not finished | +0.96 [+0.01, +1.95] | | | | |
| Kelp | independent increments | +0.60 [-0.23, +1.53] | not finished | +0.57 [-0.30, +1.52] | | | | |

"Depends on the baseline" (rule R3b): the two density estimates differ in sign, or the pre-trend interval is more than twice as wide.
Flagged: none so far.

### Simulated data with a true pre-trend (+0.5 per decade, SD 0.3 across MPAs)

| Taxon | Fitted | Data sets (usable) | True effect | Mean error | 90% intervals covering the truth |
|---|---|---|---|---|---|
| Sheephead | flat baseline | 20 (20) | +0.22 | +0.21 | 85% |
| Sheephead | with pre-trend | 20 (20) | +0.22 | -0.08 | 90% |
| Purple urchin | flat baseline | 12 (12) | -0.65 | +0.57 | 75% |
| Purple urchin | with pre-trend | 11 (11) | -0.65 | +0.26 | 100% |

Rule: coverage with the pre-trend option between 80% and 99%. **INCOMPLETE.**

