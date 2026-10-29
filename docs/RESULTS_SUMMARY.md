# Results Summary: CA MPA Kelp Forest Analysis

**Generated:** 2026-10-29 10:35:23
**Pipeline version:** Modular pBACIPS v2.0

---

## Overview

This document summarizes the key results from the California MPA kelp forest
pBACIPS/CI analysis. The primary synthesis combines before/after step contrasts,
fitted changes from model time 0 to 11, and mean log response ratios over retained CI input years.
Most diver models use survey indices; Landsat and Campus Point lobster biomass use calendar time.
The numerical horizon does not standardize calendar duration, and CI means are not uniformly post-only.
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
| S. purpuratus | 10 | -0.6978 | 0.498 | -50.2% | 0.2419 | -2.8854 | 0.0051** | 0.0115 | [0.307, 0.806] | negative |
| M. franciscanus | 14 | 0.5248 | 1.69 | +69% | 0.2247 | 2.3358 | 0.0222* | 0.0400 | [1.08, 2.644] | positive |
| M. pyrifera | 33 | 1.1303 | 3.096 | +209.6% | 0.2345 | 4.8189 | 0.0000*** | 0.0001 | [1.94, 4.941] | positive |
| P. interruptus | 5 | 0.3663 | 1.442 | +44.2% | 0.2899 | 1.2637 | 0.2103 | 0.2704 | [0.81, 2.57] | positive |
| S. pulcher | 17 | 0.96 | 2.612 | +161.2% | 0.2305 | 4.1649 | 0.0001*** | 0.0004 | [1.65, 4.134] | positive |

### Density

| Taxa | k | Estimate | RR | % Change | SE | t | p-value | p (FDR) | 95% CI (RR) | Effect |
|------|---|----------|----|----------|----|----|---------|---------|-------------|--------|
| S. purpuratus | 13 | -1.9902 | 0.137 | -86.3% | 0.5692 | -3.4963 | 0.0009*** | 0.0026 | [0.044, 0.426] | negative |
| M. franciscanus | 16 | -0.6105 | 0.543 | -45.7% | 0.5636 | -1.0832 | 0.2829 | 0.3183 | [0.176, 1.676] | negative |
| P. interruptus | 17 | 0.9515 | 2.59 | +159% | 0.5643 | 1.6861 | 0.0968 | 0.1452 | [0.838, 8.001] | positive |
| S. pulcher | 20 | -0.3195 | 0.727 | -27.3% | 0.5631 | -0.5673 | 0.5725 | 0.5725 | [0.236, 2.239] | negative |

*Significance: \*p<0.05, \*\*p<0.01, \*\*\*p<0.001 (uncorrected). p (FDR) = Benjamini-Hochberg adjusted p-values.*

### Heterogeneity (Joint Model)

tau2 and pseudo-I2 are shared across all taxa within each response type.

| Response | tau2 | I2 (%) |
|----------|------|--------|
| Biomass | 0.2742 | 79.7 |
| Density | 1.358 | 97.2 |

---

## Replicate-Level Summary

- **Significant positive contrasts:** 64
- **Significant negative contrasts:** 33
- **Not significant:** 48
- **Total MPA-level effect sizes:** 145

### Per-Taxa Breakdown

| Taxa | n | Significant | Raw Mean lnRR |
|------|---|-------------|---------------|
| P. interruptus | 22 | 13 (59%) | 1.08 |
| S. pulcher | 37 | 26 (70%) | 0.529 |
| S. purpuratus | 23 | 18 (78%) | -1.393 |
| M. franciscanus | 30 | 21 (70%) | -0.348 |
| M. pyrifera | 33 | 19 (58%) | 1.266 |

*Raw Mean lnRR is the arithmetic mean of individual MPA effect sizes (bio + den combined), not the meta-analytic estimate. See Table 2 above for meta-analytic means.*

Full replicate-level effect sizes are in `outputs/replicate_effects.csv`.

---

## Cross-Taxa Meta-Regressions (Table 3)

Knapp-Hartung meta-regressions testing whether predator MPA effects predict prey effects across sites.

| Relationship | k | Slope | SE | p | Significant? |
|-------------|---|-------|----|---|-------------|
| S. purpuratus density -> M. pyrifera biomass | 11 | -0.41 | 0.247 | 0.131 | no |
| M. franciscanus density -> M. pyrifera biomass | 10 | -0.146 | 0.133 | 0.304 | no |
| P. interruptus density -> S. purpuratus density | 11 | -0.221 | 0.256 | 0.410 | no |
| S. pulcher density -> S. purpuratus density | 11 | 0.541 | 0.458 | 0.268 | no |
| P. interruptus biomass -> S. purpuratus biomass | 4 | -0.562 | 0.115 | 0.040 | **yes** |
| S. pulcher biomass -> S. purpuratus biomass | 9 | 0.14 | 0.176 | 0.454 | no |

---

## Variance Components

| Response | Component | Estimate | 95% CI |
|----------|-----------|----------|--------|
| Biomass | MPA | 0.2087 | [0.0176, 0.9097] |
| Biomass | Source | 0.5125 | [0, 20.0609] |
| Density | MPA | 0.2706 | [0.0663, 1.0657] |
| Density | Source | 0.6521 | [0.0279, 14.5858] |

---

## Temporal Recovery Slopes (Pooled lmer Model)

Species-specific slopes from `lmer(lnRR ~ time*species + (1+time|MPA) + (1|source))`.

| Species | Slope (lnRR/yr) | SE | p |
|---------|-----------------|-----|---|
| P. interruptus | 0.1155 | 0.0231 | < 0.001 |
| S. pulcher | -0.0201 | 0.0205 | 0.343 |
| S. purpuratus | -0.1106 | 0.0223 | < 0.001 |
| M. franciscanus | -0.0672 | 0.0211 | 0.006 |
| M. pyrifera | 0.1645 | 0.0239 | < 0.001 |

---

## Outlier Sensitivity

The primary analysis uses the **full dataset** (no outlier removal, k = 145).
A global Cook's D threshold of 4/n flagged **92 of 145** observations (63%) as influential.
These flagged observations span all taxa and arise primarily from between-taxon divergence
(the expected trophic cascade signal) rather than within-taxon anomalies.

A per-taxon Cook's D threshold (4/k) flagged only **13 of 145** observations,
confirming the global threshold's over-flagging.

See Table S9 for 4-method sensitivity comparison and Figures S13-S15 for visualizations.

---

For known limitations, see `docs/ANALYSIS_REVISIONS.md`.

*This summary was auto-generated by the analysis pipeline.*
*For questions, contact Emily Donham or Adrian Stier.*
