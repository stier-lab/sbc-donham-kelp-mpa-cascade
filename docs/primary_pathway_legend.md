# Primary Consumer-Resource Pathways

## Caption

The six primary cross-taxa meta-regressions. Each point is one MPA, with primary source-specific effect sizes averaged equally within taxon and response; source units are never averaged. Panels relate (a) lobster biomass to purple urchin biomass, (b) lobster density to purple urchin density, (c) sheephead biomass to purple urchin biomass, (d) sheephead density to purple urchin density, (e) purple urchin density to kelp biomass and (f) red urchin density to kelp biomass. Lobster is Panulirus interruptus, sheephead is Semicossyphus pulcher, purple urchin is Strongylocentrotus purpuratus, red urchin is Mesocentrotus franciscanus and giant kelp is Macrocystis pyrifera.

Axes show primary MPA effect sizes as unitless log contrasts, not physical biomass or density. Their estimand follows the selected source-specific model, as recorded in the membership export. For Mean CI histories, it is the mean log contrast over retained CI input years, which can include pre-designation and designation-year observations; these histories are not uniformly post-only. Blue circles denote islands; orange triangles denote mainland MPAs. Whiskers are +/-1 pooled sampling SE in both axes. Solid lines are metafor rma fits with REML and Knapp-Hartung tests; grey bands are nominal 95% confidence intervals for the mean relationship, not prediction intervals. Fits are displayed only within the observed predictor range; panel scales differ. k counts MPAs, not source rows or annual observations. The four points in panel (a) are labeled because only four MPAs contribute.

For n contributing sources, the effect is their arithmetic mean and its sampling variance is sum(SE^2)/n^2, matching the primary Table 3 implementation. This arithmetic assumes independent source estimates. Predictor uncertainty is displayed but is not modeled by rma. Predictor and response can draw on different source sets; kelp includes satellite canopy as well as diver-derived biomass indices. These are estimated cross-sectional associations, not direct observations of predation or evidence that an MPA caused a trophic cascade.

MPAs are not wholly independent: Naples SMCA and Campus Point SMCA share mapped PISCO Arroyo Quemado reference reefs, including in purple-urchin biomass. Programme footprints can overlap, and satellite reference polygon overlap is unverified. Equal-source averaging does not remove these dependencies. The primary regression and its nominal intervals are retained without spatial/shared-reference adjustment; no new primary estimates or significance claims are introduced.

## Reproduction And Verification

Run Rscript code/R/31_primary_pathway_display.R, or source the script; all working objects and RNG/options changes remain local or are restored. Both run_review_concerns.R and run_all.R include this module. Standalone execution loads no full pipeline, raw data or external workspace.

Every fit must match Table 3 k, the intercept, slope, slope SE and slope CI endpoints to their four-decimal export rounding, plus the unrounded slope p-value within 1e-8, before any figure or table is written. A mismatch fails loudly rather than overwriting the primary table.

Export dimensions: 170 x 210 mm; PNG 600 dpi; six facets in two columns.

- Primary input: data/sumstats_final.csv; MD5 153842a6bab8d85d1152a4be5a05b29c
- Verification input: tables/table_03_cross_taxa_meta_regression.csv; MD5 5ffc08496576b6aba43ed0f6965af2cd
- Outputs: plots/fig_s_primary_pathway.pdf, plots/fig_s_primary_pathway.png, tables/table_s_primary_pathway_membership.csv, docs/primary_pathway_legend.md

The membership CSV has one row per plotted MPA/panel, with both taxa/responses, individual primary source keys, source lists, source models/designs, source counts, equal-source means, variances and SEs. verified_slope and table3_slope document concordance. All primary values are read-only.

| Panel | k | Verified slope | Table 3 slope |
|---|---:|---:|---:|
| a | 4 | -0.56160199 | -0.5616 |
| b | 11 | -0.22103263 | -0.2210 |
| c | 9 | 0.13946706 | 0.1395 |
| d | 11 | 0.54127678 | 0.5413 |
| e | 11 | -0.41042310 | -0.4104 |
| f | 10 | -0.14575539 | -0.1458 |

Available runtime: R version 4.6.2 (2026-10-29); dplyr 1.2.1; ggplot2 4.0.3; ggrepel 0.9.8; metafor 5.2.1; ragg 1.5.2; scales 1.4.0.
