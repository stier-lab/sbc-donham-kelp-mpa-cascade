# Density, mean mass and biomass in one model

2026-11-11. Review only. Nothing here is adopted, and the paper and the frozen pipeline are untouched.


## 1. What the model is

Each taxon has two linked parts fitted together. Density: the single-taxon count model under the recommended structure. Mean individual mass: station level, regional year effects (independent, shared trend, regional trend), a station AR(1) process, a siting difference, and its own response curve to protection. A group of n measured animals enters through its mean mass, Gamma(n × shape, mean M). Biomass is density × mean mass, so on the log scale the biomass effect is the density effect plus the mass effect. All effects are at ten years after designation; intervals are 90%.

Only animals measured at the station are used. No sizes are borrowed from other sites or programs and none are imputed; station-years with counts but no sizes inform density only.

## 2. Inputs: rebuilt from source rows and checked

`code/joint_model/prep_mass.py` recomputes every mass group from the registered source rows and compares it with the integrated review's own values:

- 4,035 size pools (urchins, PISCO lobster): number measured and mean mass equal the review's pool inventory (largest relative difference 1e-15; none unmatched).
- 8,260 survey units with sizes on the count record (sheephead, LTER lobster, LTER urchin groups): mean mass equals stored biomass ÷ density (largest relative difference 7e-16).
- KFM sheephead: 522 roving-diver size surveys at 33 sites (9,533 fish, 2007–2023). The 305 visits used in the November 3 KFM sheephead review are reproduced exactly.
- No group comes from an excluded site.

| Taxon | Program | Kind | Groups | Animals | Station-years |
|---|---|---|---|---|---|
| M. franciscanus | KFM | pool | 550 | 99,312 | 540 |
| M. franciscanus | LTER | unit | 409 | 1,809 | 56 |
| M. franciscanus | PISCO | pool | 875 | 14,712 | 416 |
| M. pyrifera | KFM | pool | 403 | 31,247 | 397 |
| M. pyrifera | PISCO | unit | 3,153 | 87,158 | 777 |
| P. interruptus | LTER | unit | 802 | 8,436 | 48 |
| P. interruptus | PISCO | pool | 1,335 | 6,697 | 712 |
| S. pulcher | KFM | roving | 522 | 9,533 | 487 |
| S. pulcher | LTER | unit | 214 | 761 | 60 |
| S. pulcher | PISCO | unit | 5,944 | 17,199 | 1,111 |
| S. purpuratus | KFM | pool | 549 | 106,856 | 538 |
| S. purpuratus | LTER | unit | 891 | 8,332 | 80 |
| S. purpuratus | PISCO | pool | 726 | 15,727 | 368 |

(All sites the sources hold; the model uses the groups at its own stations.)

## 3. Simulated data: does the model recover known mass and biomass effects?

Counts and group mean masses were simulated on each taxon's real survey and measurement design at known values, then refitted (4 chains, 1,500 + 1,500).

| Taxon | Data sets | Converged at R-hat < 1.05 | Mass: truth | Mass: mean error, coverage | Biomass: truth | Biomass: mean error, coverage |
|---|---|---|---|---|---|---|
| Lobster | 4 | 4 | +0.20 | -0.00, 4/4 | +2.12 | +0.02, 4/4 |
| Sheephead | 8 | 8 | +0.15 | -0.01, 7/8 | +0.38 | -0.04, 8/8 |
| Purple urchin | 4 | 3 | -0.10 | +0.10, 3/4 | -0.75 | +0.32, 4/4 |
| Red urchin | 4 | 4 | +0.12 | -0.01, 3/4 | -0.13 | +0.21, 4/4 |
| Kelp | 4 | 3 | +0.00 | +0.01, 4/4 | +0.70 | -0.47, 3/4 |

**PASS.** Coverage 0.92 over 48 intervals (rule: at least 0.81).

## 4. Real data: effect of protection at ten years

### Random-walk curve prior (main)

| Taxon | Density | Mean mass | Biomass | Biomass as % change | Converged (R-hat / eff. draws / divergences) |
|---|---|---|---|---|---|
| Lobster | +1.75 [+1.04, +2.43] | +0.38 [+0.15, +0.64] | +2.13 [+1.43, +2.88] | +740% [+318%, +1687%] | yes (1.002 / 2696 / 0) |
| Sheephead | +0.20 [-0.09, +0.50] | +0.33 [+0.10, +0.56] | +0.53 [+0.14, +0.89] | +70% [+15%, +143%] | yes (1.002 / 1512 / 0) |
| Purple urchin | -0.91 [-1.80, +0.02] | +0.25 [-0.05, +0.53] | -0.66 [-1.66, +0.25] | -48% [-81%, +28%] | yes (1.004 / 1341 / 0) |
| Red urchin | -0.63 [-1.23, -0.05] | +0.41 [+0.08, +0.73] | -0.22 [-0.90, +0.44] | -20% [-59%, +56%] | yes (1.004 / 1448 / 0) |
| Kelp | +0.96 [+0.00, +1.91] | +0.01 [-0.20, +0.22] | +0.96 [+0.01, +1.95] | +162% [+1%, +602%] | yes (1.002 / 1474 / 0) |

### Independent-increment curve prior

| Taxon | Density | Mean mass | Biomass | Biomass as % change | Converged (R-hat / eff. draws / divergences) |
|---|---|---|---|---|---|
| Lobster | +1.82 [+1.10, +2.53] | +0.23 [-0.03, +0.49] | +2.05 [+1.26, +2.81] | +675% [+253%, +1562%] | NO (1.001 / 3124 / 1) |
| Sheephead | +0.15 [-0.08, +0.41] | +0.26 [-0.00, +0.47] | +0.41 [+0.05, +0.76] | +51% [+5%, +113%] | yes (1.002 / 2210 / 0) |
| Purple urchin | -0.28 [-0.97, +0.37] | +0.09 [-0.10, +0.33] | -0.19 [-0.92, +0.50] | -17% [-60%, +65%] | yes (1.002 / 1626 / 0) |
| Red urchin | -0.27 [-0.76, +0.20] | +0.21 [-0.07, +0.52] | -0.06 [-0.66, +0.52] | -6% [-48%, +69%] | yes (1.005 / 1171 / 0) |
| Kelp | +0.60 [-0.23, +1.53] | -0.03 [-0.23, +0.14] | +0.57 [-0.30, +1.52] | +77% [-26%, +356%] | yes (1.004 / 1197 / 0) |

Seeds: largest difference between two independent fits 0.030 (rule: below 0.05) — **PASS**.

Mass data behind each estimate (model stations only):

| Taxon | Station-years with measured animals | Station-years with density | Groups by program | MPAs with before-and-after mass data |
|---|---|---|---|---|
| Lobster | 489 | 1,419 |  | 0 |
| Sheephead | 1,554 | 1,892 | KFM 413, LTER 214, PISCO 5,944 | 3 |
| Purple urchin | 812 | 1,515 |  | 4 |
| Red urchin | 804 | 1,515 |  | 3 |
| Kelp | 1,032 | 1,421 |  | 3 |

## 5. The open biomass issues, one by one

Each row re-fits the model with one thing changed and gives the shift in the effect at ten years against the main fit (random-walk prior). Shifts above 0.10 are marked ⚠.

### 5.1 Allometry

The length-to-mass exponent is fixed at the published value. These fits move it by −0.1 and +0.1 (for an exponent near 2.9 to 3.0 that is a large error). The intercept cancels in every effect.

| Taxon | Change | Density | Mean mass | Biomass | Shift in biomass effect | Converged |
|---|---|---|---|---|---|---|
| Lobster | exponent −0.1 | +1.74 [+1.03, +2.42] | +0.37 [+0.12, +0.60] | +2.11 [+1.35, +2.81] | -0.02 | NO (1.002 / 2829 / 1) |
| Lobster | exponent +0.1 | +1.75 [+1.07, +2.45] | +0.40 [+0.16, +0.66] | +2.15 [+1.44, +2.90] | +0.02 | yes (1.002 / 2490 / 0) |
| Sheephead | exponent −0.1 | +0.20 [-0.09, +0.49] | +0.31 [+0.11, +0.55] | +0.51 [+0.15, +0.88] | -0.02 | yes (1.002 / 1606 / 0) |
| Sheephead | exponent +0.1 | +0.20 [-0.10, +0.48] | +0.34 [+0.12, +0.58] | +0.54 [+0.16, +0.90] | +0.01 | yes (1.002 / 1535 / 0) |
| Purple urchin | exponent −0.1 | -0.92 [-1.82, -0.00] | +0.25 [-0.03, +0.53] | -0.67 [-1.66, +0.26] | -0.01 | yes (1.004 / 1389 / 0) |
| Purple urchin | exponent +0.1 | -0.91 [-1.82, -0.01] | +0.26 [-0.03, +0.56] | -0.65 [-1.63, +0.28] | +0.01 | yes (1.004 / 1350 / 0) |
| Red urchin | exponent −0.1 | -0.62 [-1.20, -0.02] | +0.40 [+0.10, +0.72] | -0.22 [-0.90, +0.43] | +0.00 | yes (1.008 / 1276 / 0) |
| Red urchin | exponent +0.1 | -0.63 [-1.19, -0.01] | +0.42 [+0.08, +0.74] | -0.21 [-0.88, +0.47] | +0.01 | yes (1.003 / 1424 / 0) |

Average shift in log mean mass from the exponent change (it moves levels a lot and contrasts little): M. franciscanus -0.42 / +0.42; P. interruptus -0.44 / +0.44; S. pulcher -0.35 / +0.35; S. purpuratus -0.36 / +0.36.

### 5.2 One size definition for urchins (25 mm and larger)

PISCO counts and measures urchins of 25 mm and up; KFM and LTER count all sizes. These fits put KFM and LTER on the PISCO definition: mean mass from measured animals of 25 mm and up, and density scaled by the share of measured animals of 25 mm and up at that station and year. Shares of measured animals at 25 mm and up: M. franciscanus KFM 83%; M. franciscanus LTER 97%; S. purpuratus KFM 56%; S. purpuratus LTER 88%.

| Taxon | Change | Density | Mean mass | Biomass | Shift in biomass effect | Converged |
|---|---|---|---|---|---|---|
| Purple urchin | 25 mm and up, rw | -0.95 [-1.88, -0.06] | +0.19 [-0.08, +0.44] | -0.76 [-1.69, +0.20] | -0.10 ⚠ | yes (1.002 / 939 / 0) |
| Purple urchin | 25 mm and up, piecewise | -0.31 [-1.04, +0.35] | +0.06 [-0.10, +0.25] | -0.26 [-1.01, +0.45] | -0.07 | yes (1.003 / 1050 / 0) |
| Red urchin | 25 mm and up, rw | -0.65 [-1.24, -0.08] | +0.38 [+0.06, +0.72] | -0.26 [-0.92, +0.41] | -0.04 | yes (1.004 / 1390 / 0) |
| Red urchin | 25 mm and up, piecewise | -0.27 [-0.79, +0.18] | +0.22 [-0.06, +0.52] | -0.06 [-0.63, +0.55] | +0.00 | yes (1.006 / 1066 / 0) |

### 5.3 LTER urchin sizes are group means

LTER records a mean size for each group of urchins, so the spread inside a group is unknown. These fits drop LTER urchin sizes from the mass part entirely (LTER still informs density).

| Taxon | Change | Density | Mean mass | Biomass | Shift in biomass effect | Converged |
|---|---|---|---|---|---|---|
| Purple urchin | no LTER sizes | -0.93 [-1.89, -0.06] | +0.23 [-0.09, +0.54] | -0.70 [-1.70, +0.24] | -0.04 | yes (1.004 / 1286 / 0) |
| Red urchin | no LTER sizes | -0.62 [-1.21, -0.04] | +0.42 [+0.06, +0.75] | -0.20 [-0.87, +0.49] | +0.03 | yes (1.003 / 1337 / 0) |

### 5.4 KFM sheephead sizes

The main sheephead fit now includes KFM roving-diver sizes (2007–2023). This row leaves them out, as before today. Roving-diver sizes and transect counts are different samples of the fish at a site; that mismatch is not removed. KFM lobster has no sizes at all and informs density only.

| Taxon | Change | Density | Mean mass | Biomass | Shift in biomass effect | Converged |
|---|---|---|---|---|---|---|
| Sheephead | without KFM sizes | +0.20 [-0.08, +0.50] | +0.29 [+0.06, +0.51] | +0.49 [+0.13, +0.85] | -0.04 | yes (1.002 / 2144 / 0) |

### 5.5 Density and mass treated as independent

These fits let mean mass at a station move with that year's density fluctuation there (for example a recruitment pulse raising numbers and lowering mean size).

| Taxon | Change | Density | Mean mass | Biomass | Shift in biomass effect | Converged |
|---|---|---|---|---|---|---|
| Lobster | mass linked to density fluctuation | +1.76 [+1.04, +2.42] | +0.39 [+0.14, +0.63] | +2.15 [+1.43, +2.92] | +0.02 | yes (1.002 / 1820 / 0) |
| Sheephead | mass linked to density fluctuation | +0.20 [-0.08, +0.50] | +0.34 [+0.10, +0.57] | +0.54 [+0.18, +0.88] | +0.01 | yes (1.002 / 1923 / 0) |
| Purple urchin | mass linked to density fluctuation | -0.92 [-1.81, -0.04] | +0.25 [-0.02, +0.51] | -0.67 [-1.62, +0.17] | -0.01 | yes (1.006 / 1221 / 0) |
| Red urchin | mass linked to density fluctuation | -0.64 [-1.20, -0.05] | +0.40 [+0.09, +0.72] | -0.23 [-0.89, +0.39] | -0.01 | yes (1.006 / 1241 / 0) |
| Kelp | mass linked to density fluctuation | +0.96 [+0.01, +1.94] | -0.02 [-0.23, +0.19] | +0.94 [-0.02, +1.91] | -0.02 | yes (1.003 / 1273 / 0) |

Link coefficient (change in log mean mass per log unit of density fluctuation): Lobster +0.08 [-0.16, +0.32]; Sheephead -0.26 [-0.37, -0.14]; Purple urchin -0.15 [-0.19, -0.11]; Red urchin -0.21 [-0.28, -0.13]; Kelp -0.05 [-0.09, -0.01].

### 5.6 KFM reference sites flagged for fish

Scott Gabara flagged Webster's Arch, Miracle Mile and Wyckoff Ledge as questionable references for fish. These fits leave each out, then all three. They are sensitivities; no site is dropped from the main fit.

| Taxon | Change | Density | Mean mass | Biomass | Shift in biomass effect | Converged |
|---|---|---|---|---|---|---|
| Sheephead | without Webster's Arch | +0.20 [-0.09, +0.49] | +0.32 [+0.10, +0.55] | +0.52 [+0.13, +0.87] | -0.01 | yes (1.003 / 1929 / 0) |
| Sheephead | without Miracle Mile | +0.20 [-0.10, +0.48] | +0.33 [+0.09, +0.55] | +0.53 [+0.16, +0.90] | -0.00 | yes (1.002 / 2111 / 0) |
| Sheephead | without Wyckoff Ledge | +0.24 [-0.05, +0.53] | +0.33 [+0.11, +0.57] | +0.57 [+0.20, +0.94] | +0.04 | yes (1.002 / 1703 / 0) |
| Sheephead | without all three | +0.23 [-0.06, +0.52] | +0.34 [+0.12, +0.58] | +0.57 [+0.20, +0.95] | +0.04 | yes (1.002 / 1612 / 0) |

### 5.7 Prior on the mass response

Prior scale of the mass response doubled.

| Taxon | Change | Density | Mean mass | Biomass | Shift in biomass effect | Converged |
|---|---|---|---|---|---|---|
| Lobster | mass prior ×2 | +1.75 [+1.03, +2.42] | +0.40 [+0.15, +0.66] | +2.15 [+1.42, +2.90] | +0.02 | NO (1.002 / 2752 / 2) |
| Sheephead | mass prior ×2 | +0.20 [-0.08, +0.49] | +0.31 [+0.08, +0.56] | +0.51 [+0.13, +0.87] | -0.02 | yes (1.002 / 1998 / 0) |
| Purple urchin | mass prior ×2 | -0.91 [-1.82, -0.00] | +0.27 [-0.03, +0.58] | -0.65 [-1.63, +0.29] | +0.01 | yes (1.004 / 1375 / 0) |
| Red urchin | mass prior ×2 | -0.63 [-1.23, -0.06] | +0.44 [+0.10, +0.78] | -0.18 [-0.86, +0.50] | +0.04 | yes (1.004 / 1398 / 0) |
| Kelp | mass prior ×2 | +0.96 [+0.04, +1.94] | -0.01 [-0.22, +0.21] | +0.95 [+0.01, +1.95] | -0.02 | yes (1.004 / 1175 / 0) |

**Sensitivities that move the biomass effect by more than 0.10:** Purple urchin: 25 mm and up, rw (-0.10).

## 6. What is still not fixed

- Allometry coefficients have no published uncertainty in the files we hold; the exponent check above is a bound, not a propagated uncertainty.
- KFM lobster has no size data. Lobster mass and biomass rest on LTER and PISCO.
- KFM sheephead sizes come from roving divers and counts from transects.
- Urchin mass has few MPAs with before-and-after size data, so the urchin mass effect leans on post-designation contrasts and the siting-difference prior.
- Kelp: LTER counts stipes with a fixed converted mass per stipe, so LTER kelp informs density only; KFM and PISCO count plants.
- The density model's known weak points carry over: the urchin count layer misfit and the sheephead time placebo.
- NPS permission for redistributing KFM records is a release question and is unchanged.

