# Results Summary: CA MPA Kelp Forest Analysis

**Generated:** 2026-11-01 09:25:22
**Pipeline version:** Modular pBACIPS v2.0

---

## Overview

This document summarizes the key results from the California MPA kelp forest
pBACIPS/CI analysis. The primary synthesis combines before/after step contrasts,
fitted changes between surveyed post-designation years, and mean log response ratios over retained CI input years.
All trend models use elapsed calendar years; contrasts span at most eleven years from the first post-designation survey.
Observed comparison durations can differ, and CI means are not uniformly post-only.
The legacy RR column exponentiates these contrasts; it is not a common MPA/reference level ratio.
The kelp Biomass category combines diver biomass and satellite canopy-index contrasts.

**Interpretation:**
- **Estimate > 0 / exp(Estimate) > 1** = positive pooled design-specific contrast
- **Estimate < 0 / exp(Estimate) < 1** = negative pooled design-specific contrast
- **Estimate = 0 / exp(Estimate) = 1** = zero pooled contrast
- **% Change** = (exp(Estimate) - 1) x 100, not a common percentage increase inside MPAs
- Significance threshold: p < 0.05

---

## Meta-Analysis Results (Table 2)

**Estimate** = pooled log-scale contrast. **RR** = exp(Estimate), retaining the legacy column name.
Values above or below 1 denote positive or negative pooled contrasts, not a common inside/reference difference.

### Biomass

| Taxa | k | Estimate | RR | % Change | SE | t | p-value | p (FDR) | 95% CI (RR) | Effect |
|------|---|----------|----|----------|----|----|---------|---------|-------------|--------|
| S. purpuratus | 10 | -0.5842 | 0.558 | -44.2% | 0.2716 | -2.1512 | 0.0348* | 0.0522 | [0.324, 0.958] | negative |
| M. franciscanus | 14 | 0.571 | 1.77 | +77% | 0.2556 | 2.2342 | 0.0286* | 0.0514 | [1.063, 2.946] | positive |
| M. pyrifera | 31 | 1.1379 | 3.12 | +212% | 0.2631 | 4.3246 | 0.0000*** | 0.0004 | [1.847, 5.272] | positive |
| P. interruptus | 5 | 0.4163 | 1.516 | +51.6% | 0.3136 | 1.3274 | 0.1886 | 0.2425 | [0.811, 2.833] | positive |
| S. pulcher | 17 | 0.9298 | 2.534 | +153.4% | 0.261 | 3.5631 | 0.0007*** | 0.0029 | [1.506, 4.263] | positive |

### Density

| Taxa | k | Estimate | RR | % Change | SE | t | p-value | p (FDR) | 95% CI (RR) | Effect |
|------|---|----------|----|----------|----|----|---------|---------|-------------|--------|
| S. purpuratus | 13 | -1.5595 | 0.21 | -79% | 0.4697 | -3.3203 | 0.0015** | 0.0045 | [0.082, 0.538] | negative |
| M. franciscanus | 16 | -0.4648 | 0.628 | -37.2% | 0.4637 | -1.0023 | 0.3201 | 0.3601 | [0.249, 1.588] | negative |
| P. interruptus | 17 | 1.0485 | 2.853 | +185.3% | 0.4649 | 2.2555 | 0.0276* | 0.0514 | [1.127, 7.226] | positive |
| S. pulcher | 20 | -0.2481 | 0.78 | -22% | 0.4632 | -0.5355 | 0.5942 | 0.5942 | [0.309, 1.97] | negative |

*Significance: \*p<0.05, \*\*p<0.01, \*\*\*p<0.001 (uncorrected). p (FDR) = Benjamini-Hochberg adjusted p-values.*

### Heterogeneity (Joint Model)

tau2 and pseudo-I2 are shared across all taxa within each response type.

| Response | tau2 | I2 (%) |
|----------|------|--------|
| Biomass | 0.3445 | 84.1 |
| Density | 0.888 | 95.6 |

---

## Replicate-Level Summary

- **Significant positive contrasts:** 60
- **Significant negative contrasts:** 33
- **Not significant:** 50
- **Total MPA-level effect sizes:** 143

### Per-Taxa Breakdown

| Taxa | n | Significant | Raw Mean lnRR |
|------|---|-------------|---------------|
| P. interruptus | 22 | 13 (59%) | 1.036 |
| S. pulcher | 37 | 23 (62%) | 0.549 |
| S. purpuratus | 23 | 18 (78%) | -1.3 |
| M. franciscanus | 30 | 21 (70%) | -0.159 |
| M. pyrifera | 31 | 18 (58%) | 1.22 |

*Raw Mean lnRR is the arithmetic mean of individual MPA effect sizes (bio + den combined), not the meta-analytic estimate. See Table 2 above for meta-analytic means.*

Full replicate-level effect sizes are in `outputs/replicate_effects.csv`.

---

## Cross-Taxa Meta-Regressions (Table 3)

Knapp-Hartung meta-regressions testing whether predator MPA effects predict prey effects across sites.

| Relationship | k | Slope | SE | p | Significant? |
|-------------|---|-------|----|---|-------------|
| S. purpuratus density -> M. pyrifera biomass | 11 | -0.328 | 0.29 | 0.288 | no |
| M. franciscanus density -> M. pyrifera biomass | 10 | -0.21 | 0.238 | 0.403 | no |
| P. interruptus density -> S. purpuratus density | 11 | -0.434 | 0.2 | 0.058 | no |
| S. pulcher density -> S. purpuratus density | 11 | -0.091 | 0.785 | 0.911 | no |
| P. interruptus biomass -> S. purpuratus biomass | 4 | -0.555 | 0.111 | 0.038 | **yes** |
| S. pulcher biomass -> S. purpuratus biomass | 9 | 0.093 | 0.174 | 0.611 | no |

---

## Variance Components

| Response | Component | Estimate | 95% CI |
|----------|-----------|----------|--------|
| Biomass | MPA | 0.0204 | [0, 0.3346] |
| Biomass | Source | 0.3181 | [0.0448, 3.2481] |
| Density | MPA | 0.3408 | [0.0354, 1.642] |
| Density | Source | 0.8032 | [0, 18.7515] |

---

## Temporal Recovery Slopes (Pooled lmer Model)

Species-specific slopes from `lmer(lnRR ~ time*species + (1+time|MPA) + (1|source))`.

| Species | Slope (lnRR/yr) | SE | p |
|---------|-----------------|-----|---|
| P. interruptus | 0.1106 | 0.0173 | < 0.001 |
| S. pulcher | -0.0229 | 0.0133 | 0.085 |
| S. purpuratus | -0.1187 | 0.016 | < 0.001 |
| M. franciscanus | -0.0717 | 0.0139 | < 0.001 |
| M. pyrifera | 0.1595 | 0.0181 | < 0.001 |

---

## Outlier Sensitivity

The primary analysis uses the **full dataset** (no outlier removal, k = 143).
A global Cook's D threshold of 4/n flagged **93 of 143** observations (65%) as influential.
These flagged observations span all taxa and arise primarily from between-taxon divergence
(the expected trophic cascade signal) rather than within-taxon anomalies.

A per-taxon Cook's D threshold (4/k) flagged only **12 of 143** observations,
confirming the global threshold's over-flagging.

See Table S9 for 4-method sensitivity comparison and Figures S13-S15 for visualizations.

---

For known limitations, see `docs/ANALYSIS_REVISIONS.md`.

*This summary was auto-generated by the analysis pipeline.*
*For questions, contact Emily Donham or Adrian Stier.*
