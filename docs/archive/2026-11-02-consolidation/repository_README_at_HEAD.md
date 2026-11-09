# Kelp Resilience and Trophic Responses in California MPAs

[![R](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)](https://www.r-project.org/)
[![Status](https://img.shields.io/badge/Status-In%20Preparation-yellow.svg)]()

> **Manuscript Title:** Kelp resilience and trophic responses in a southern California marine protected-area network
>
> **Target journal:** Journal of Applied Ecology (Research Article)
>
> *Scientific overview updated: 2026-10-20. The primary synthesis retains 145 histories from 21 MPAs. Population patterns, trophic links and designation-related disturbance responses are separate questions; current framing is in `docs/SCIENTIFIC_APPROACH.md`.*

---

> **For co-authors — the short version.** You don't need to run anything to write the
> paper. Figures are in `plots/`, tables in `tables/`, and plain-language results
> summaries in `docs/` (start with `docs/RESULTS_SUMMARY.md`). The main paper analyses
> are scripts `08`-`13`, the in-pipeline resilience modules `14/15/16/19/21/22/23`,
> and the PISCO-reference sensitivity in script `26`; Figure 5 comes from script `22`.
> Scripts `27`-`30` add matched design/geography comparisons, a same-pair
> historical heatwave comparison, the complete pair-history atlas and arithmetic checks.
> Script `34` adds Gull Island KFM omissions and a common-cohort combined-predator check. See `docs/COAUTHOR_SENSITIVITY_INTEGRATION_2026-10-12.md` for methods, results and limits. Sheephead display names use Bodianus; historical analysis identifiers are unchanged.
>
> Modules `36` and `37` check exact primary-contrast and Figure 4 uncertainty. Following the singleton correction, calendar-HAC checks preserve all nine pooled directions, but adjusted support for purple-urchin biomass depends on bandwidth. Figure 4 uses 831 biomass observations in 63 retained diver histories; its full random-slope/source model is nonsingular. Same-row calendar AR1 retains a positive kelp slope, but the purple-urchin slope's interval includes zero. See `docs/BIOMASS_SINGLETON_REFIT_2026-10-18.md`; the retired pooled DW reconstruction is not primary-model validation.
> Module `42` compares programme-separated and averaged pathways on identical MPAs and post-designation years. The broader PISCO/KFM comparison has five island MPAs for four density-based links; all full-cohort intervals include zero, with site-sensitive omissions. Strict primary-history and broader annual-record scopes are distinct; sparse LTER and biomass comparisons are explicitly unestimable. See `docs/PROGRAMME_PATHWAY_REFITS_2026-10-17.md`. Primary results are unchanged.
> Module `44` tests the supplied PISCO sheephead lookup's 10 cm young-of-year size proxy, retaining the exact primary sampling and model routes. Removing fish <10 cm or <=10 cm preserves all nine pooled directions and BH classifications; local annual contributions can be larger than the pooled raw shares. These are size proxies, not verified ages or an adult-only cohort. See `docs/SHEEPHEAD_YOY_2026-10-19.md`.
> Module `45` implements the non-target sessile-invertebrate comparability check on fixed legacy KFM reefs. Four core and three secondary taxa are matched to three kelp or four urchin MPA calendars. Sparse histories and heterogeneous starting differences do not establish reef equivalence; the check is not a siting correction or a network-wide test. See `docs/NONFISHED_COMPARISON_PLAN.md` and `docs/NONFISHED_COMPARISON_RESULTS.md`. Run `Rscript code/R/45_nonfished_comparators.R` and `Rscript code/tests/test_nonfished_comparators.R`; optional raw paths use `NONFISHED_KFM_COUNTS` and `NONFISHED_KFM_SITES`. Primary outputs are unchanged.
> Scripts `31`-`32` check the six primary pathway displays and all 37 retained
> pBACIPS fitted patterns. The October 10 refit leaves pooled lobster density
> inconclusive; see `docs/STEP_DIRECTION_REFIT_2026-10-10.md`. The kelp disturbance
> and resilience results are unchanged. Progressive models have source-specific
> survey-index or calendar-time clocks, not one common eleven-year duration.
> Script `24` checks that the resilience tables, figures, and headline values still
> match the manuscript and supplement. A plain tour is in `docs/RESILIENCE_SYNTHESIS.md`,
> and the audit trail is in
> `docs/RESILIENCE_PIPELINE_INTEGRATION.md`. The `CLAUDE.md` file in the root is a technical config
> for an AI coding assistant — you can ignore it. Questions on any of it — just ask Adrian.

---

## Overview

This repository contains the analysis pipeline for comparing predator, urchin and kelp populations in southern California marine protected areas (MPAs) and reference reefs. The primary pBACIPS/CI synthesis combines different population contrasts; it is not a common percentage effect of designation. Separate analyses examine whether trophic responses covary, how diver-survey kelp performed during and after the 2014–2016 heatwave, and whether satellite-observed disturbance responses differed before and after designation.

Some population directions match cascade predictions, and baseline-relative diver kelp biomass was higher inside MPAs during and after the heatwave. Cross-taxa links and a designation-related improvement in disturbance responses remain unresolved. Non-random siting, heterogeneous measurements and incomplete satellite-reference reconstruction limit attribution. These results do not establish either a complete predator-mediated cascade or its absence. See `docs/SCIENTIFIC_APPROACH.md` for the evidence hierarchy and remaining decisions.

This is the public analysis/code repository. Current Word manuscripts, cover
letters, Supporting Information construction, citation-audit working files, and
publisher PDF source copies live outside this public repo in the private
manuscript workspace.

### Key Findings

The analysis examines:
- **Predator recovery** (lobster, sheephead) inside vs. outside MPAs
- **Urchin population patterns**, including direct protection of harvested red urchins
- **Kelp forest biomass** trajectories following MPA establishment
- **Covariation among trophic responses**, distinct from a causal test of a complete cascade

---

## Authors

Adrian Stier<sup>1</sup>\*, Emily Donham<sup>1</sup>, Tom Bell<sup>2</sup>, Jennifer E. Caselle<sup>3</sup>, Scott Gabara<sup>4</sup>,
Scott Hamilton<sup>5</sup>, Li Kui<sup>3</sup>, Hunter Lenihan<sup>8</sup>, Robert Miller<sup>3</sup>, Dan Pondella<sup>6</sup>,
and Peter Raimondi<sup>7</sup>

<sup>1</sup> Department of Ecology, Evolution, and Marine Biology, University of California Santa Barbara, Santa Barbara, CA 93106, USA
<sup>2</sup> Department of Estuarine and Ocean Sciences, University of Massachusetts Dartmouth, New Bedford, MA 02744, USA
<sup>3</sup> Marine Science Institute, University of California Santa Barbara, Santa Barbara, CA 93106, USA
<sup>4</sup> National Park Service
<sup>5</sup> Moss Landing Marine Laboratories, San Jose State University, Moss Landing, CA 95039, USA
<sup>6</sup> Occidental College, Los Angeles, CA 90041, USA
<sup>7</sup> Department of Ecology and Evolutionary Biology, University of California Santa Cruz, Santa Cruz, CA 95060, USA
<sup>8</sup> Bren School of Environmental Science & Management, University of California Santa Barbara, Santa Barbara, CA 93106, USA

\* Corresponding author. Author order follows the current coauthor-review package; confirm
final affiliations and CRediT roles with the author team before submission.

Dan Reed is acknowledged for detailed comments and ecological advice, in accordance with his contributor preference; he is not listed as an author.

---

## Repository Structure

This is the **analysis repo**. It consumes harmonized CSVs produced by the sibling data-processing repo ([Donham-Stier-CA-MPA-Data-2026](https://github.com/stier-lab/Donham-Stier-CA-MPA-Data-2026)). The harmonized CSVs (~1 MB) are tracked in git, so no raw monitoring data is needed to run the main analysis.

```
Donham-Stier-CA-MPA-2026/
│
├── code/
│   └── R/                            # Analysis scripts
│       ├── 00_libraries.R            # Package dependencies
│       ├── 00b_color_palette.R       # Color scheme & ggplot theme
│       ├── 00c_analysis_constants.R  # Named constants & site exclusions
│       ├── 01_utils.R                # Utility functions
│       ├── 02_pBACIPS_function.R     # Core pBACIPS methodology
│       ├── 03_load_harmonized_data.R # Load harmonized CSVs
│       ├── 08_effect_sizes.R         # Calculate effect sizes
│       ├── 09_meta_analysis.R        # Multilevel meta-analysis
│       ├── 10_temporal_analysis.R    # Temporal dynamics (SI Figs S3-S5)
│       ├── 11_figures.R              # Publication figures (Figs 1-4, S1-S2, S7a-e)
│       ├── 12_results_summary.R      # Generate results summaries
│       ├── 13_additional_analyses.R  # Moderator analyses (SI Fig S6)
│       ├── 14-23_*.R                 # Resilience / heatwave robustness suite
│       ├── resilience_modules.R      # Shared 14-23 module membership
│       ├── 24_resilience_pipeline_check.R # Paper/SI resilience concordance gate
│       ├── 25_prereview_sensitivities.R # Isolated temporal and permutation sensitivity checks
│       ├── 26_pisco_control_sensitivity.R # PISCO reference-control dependence check
│       ├── 27_design_geography_sensitivity.R # Same-history change versus level; source/region checks
│       ├── 28_historical_heatwave_bacips.R # 1997-1998 versus 2014-2016 at the same Landsat pairs
│       ├── 29_pair_history_atlas.R    # All 145 primary histories, kept source-specific
│       ├── 30_review_concern_checks.R # Coverage, arithmetic and integration assertions
│       ├── 31_primary_pathway_display.R # Primary consumer-prey fits and membership
│       ├── 32_pbacips_shape_display.R # Verified selected fits and nonlinear comparison
│       ├── run_review_concerns.R     # Rebuild these checks from tracked inputs
│       ├── run_all.R                 # Pipeline orchestration
│       ├── run_resilience.R          # Full resilience suite
│       └── run_figures_only.R        # Fast cached regeneration for core figures
│
├── data/
│   ├── harmonized/                   # Analysis-ready CSVs (tracked in git)
│   └── cache/                        # Bootstrap results & intermediates
│
├── tables/                           # Generated manuscript tables (CSV)
├── plots/                            # Generated figures (PDF & PNG at 600 DPI)
├── outputs/                          # Filter audits, data flow summaries
├── docs/                             # Documentation & manuscript drafts
└── README.md                         # This file
```

---

## Data Sources

This analysis integrates long-term monitoring data from three major programs:

| Source | Full Name | Coverage | Data Type |
|--------|-----------|----------|-----------|
| **PISCO** | Partnership for Interdisciplinary Studies of Coastal Oceans | Statewide | Fish, invertebrate, algae surveys |
| **KFM/MBON** | Kelp Forest Monitoring / Marine Biodiversity Observation Network | Channel Islands | Long-term ecological monitoring |
| **LTER** | Santa Barbara Coastal Long-Term Ecological Research | Santa Barbara Channel | Comprehensive ecosystem data |
| **Landsat** | NASA/USGS Satellite Imagery | Statewide | Kelp canopy cover (remote sensing) |

### Data

**No raw data is needed to run this analysis.** The harmonized CSVs in `data/harmonized/` are tracked in git (~1 MB total). Raw monitoring data (~1.3 GB) lives in the sibling [data-processing repo](https://github.com/stier-lab/Donham-Stier-CA-MPA-Data-2026).

PISCO and LTER sheephead biomass use the same verified individual-mass equation:
`W(g) = 0.0144 * TL(cm)^3.04` (Hamilton & Caselle 2015,
https://doi.org/10.1098/rspb.2014.1817). Size-class counts weight the masses before
aggregation. The original PISCO spreadsheet is retained; a documented override
is applied during harmonization. `data/harmonized/sheephead_allometry.json`
records the correction and input fingerprints, checked by module 35 before
analysis and by the figures-only runner. See
`docs/SHEEPHEAD_ALLOMETRY_REFIT_2026-10-12.md` for the before/after audit.

---

## Setup Instructions

### Prerequisites

- **R** (>= 4.0)
- **RStudio** (recommended)

### 1. Clone the Repository

```bash
git clone https://github.com/stier-lab/sbc-2026-donham-kelp-mpa-cascade.git Donham-Stier-CA-MPA-2026
cd Donham-Stier-CA-MPA-2026
```

### 2. Install R Dependencies

```r
# Recommended: restore exact package versions from lockfile
renv::restore()

# Alternative: install packages manually
source("code/R/00_libraries.R")
```

Package versions are captured in `renv.lock` for reproducibility. Key packages include:
`tidyverse`, `here`, `metafor`, `emmeans`, `lme4`, `cowplot`, `patchwork`, `sf`

### 3. Run the Analysis

```r
# Full pipeline (runtime depends on sensitivity and diagnostic fits)
source(here::here("code", "R", "run_all.R"))

# Or regenerate the cached core figure suite (selected figures are faster)
source(here::here("code", "R", "run_figures_only.R"))
```

No raw monitoring data setup is required for the main pipeline. The small public
Kumagai comparator inputs needed by scripts `15`, `16`, `18`, and `24` are
tracked under `data/external/kumagai2024/` with provenance. The full resilience
suite can also use optional raw PISCO inputs for the standalone Eisaguirre and
SSWD modules; set `DONHAM_DATA_REPO_PISCO` if those files are not in the default
local development location.

---

## Running the Analysis

### Full Pipeline

```r
# Run complete analysis (runtime depends on sensitivity and diagnostic fits)
source(here::here("code", "R", "run_all.R"))
```

The pipeline loads harmonized CSVs, calculates effect sizes, runs meta-analysis, generates figures, runs the in-pipeline resilience subset, checks resilience concordance against the manuscript/SI claims, and produces results summaries. Prereview sensitivities run after the primary meta-analysis; the fitted-history displays verify reconstructed contrasts against the primary results before exporting. Progress is logged to the console and a timestamped file in `logs/`.

### Individual Scripts

```r
# Fast cached figure regeneration for the core figure suite
source(here::here("code", "R", "run_figures_only.R"))

# Render only specific figures
RENDER_FIGURES <- c("fig03", "fig04")
source(here::here("code", "R", "run_figures_only.R"))
```

---

## Analysis Pipeline

```
  ┌─────────────────────────────────┐
  │  data/harmonized/ (4 CSVs)      │
  │  Produced by data-processing    │
  │  repo; tracked in git           │
  └───────────────┬─────────────────┘
                  ▼
  ┌─────────────────────────────────┐
  │ 03_load_harmonized_data.R       │
  │ Load & validate CSVs            │
  └───────────────┬─────────────────┘
                  ▼
  ┌─────────────────────────────────┐
  │ 08_effect_sizes.R               │
  │ pBACIPS effect sizes            │
  │ 4 candidate models per MPA      │
  └───────────────┬─────────────────┘
                  ▼
  ┌─────────────────────────────────┐
  │ 09_meta_analysis.R              │
  │ Multilevel meta-analysis        │
  │ (metafor::rma.mv)               │
  └───────────────┬─────────────────┘
                  ▼
  ┌─────────────────────────────────┐
  │ 10_temporal_analysis.R          │
  │ Temporal dynamics, SI Figs S3-S5│
  └───────────────┬─────────────────┘
                  ▼
  ┌─────────────────────────────────┐
  │ 11_figures.R                    │
  │ Publication figures             │
  │ Figs 1-4, S1-S2, S7a-e         │
  └───────────────┬─────────────────┘
                  ▼
  ┌─────────────────────────────────┐
  │ 13_additional_analyses.R        │
  │ Moderator analyses, SI Fig S6   │
  └───────────────┬─────────────────┘
                  ▼
  ┌─────────────────────────────────┐
  │ 14/15/16/19/21/22/23            │
  │ Resilience + heatwave robustness│
  └───────────────┬─────────────────┘
                  ▼
  ┌─────────────────────────────────┐
  │ 24_resilience_pipeline_check.R  │
  │ Figure 5 + SI concordance gate  │
  └───────────────┬─────────────────┘
                  ▼
  ┌─────────────────────────────────┐
  │ 26_pisco_control_sensitivity.R  │
  │ PISCO reference-dependence check│
  └───────────────┬─────────────────┘
                  ▼
  ┌─────────────────────────────────┐
  │ 12_results_summary.R            │
  │ Auto-generate results CSVs      │
  └─────────────────────────────────┘
```

---

## Outputs

### Main Text Figures

| Figure | Description | Output File |
|--------|-------------|-------------|
| **Figure 1** | Map of MPAs with Channel Islands | `plots/fig_01_mpa_map.pdf` |
| **Figure 2** | Trophic cascade case studies: 3×3 grid (predators/urchins/kelp × 3 sites), before/after MPA | `plots/fig_02_cascade_case_studies.pdf` |
| **Figure 3** | Meta-analytic mean effect sizes by taxa, RR-scaled axis | `plots/fig_03_mean_effects.pdf` |
| **Figure 4** | Recovery trajectories: 3×2 trophic grid (biomass), lmer predictions with 95% CI | `plots/fig_04_recovery_curves.pdf` |
| **Figure 5** | Giant-kelp resistance and recovery through the 2014-16 marine heatwave | `plots/fig_kelp_resilience.pdf` |

### Supplemental Figures

| Figure | Description | Output File |
|--------|-------------|-------------|
| **Figure S1** | Data processing pipeline (raw → standardized → lnRR) | `plots/fig_s01_data_processing.pdf` |
| **Figure S2** | Forest plot: effect sizes by MPA, RR-scaled axis | `plots/fig_s02_forest_plot.pdf` |
| **Figure S3** | Species-level GAM recovery curves with lmer overlay (t≤15) | `plots/fig_s04_recovery_curves.pdf` |
| **Figure S4** | Species-pair phase portraits of trophic cascade | `plots/fig_s05_cascade_phase.pdf` |
| **Figure S5** | Per-MPA slope comparison and cascade consistency | `plots/fig_s06_slope_comparison.pdf` |
| **Figure S6** | Combined moderator comparisons: SMR vs SMCA, Islands vs mainland | `plots/fig_s09_moderator_comparisons.pdf` |
| **Figure S7a-e** | Site-level appendix (5 taxa-specific panels) | `plots/fig_s08_appendix_*.pdf` |
| **Figure S8** | DHARMa residual diagnostics for NLS models | `plots/fig_s11_dharma_diagnostics.pdf` |
| **Figure S9** | Funnel plots + Egger's test for publication bias | `plots/fig_s12_funnel_plots.pdf` |
| **Figure S10** | lmer residual diagnostics (4-panel) | `plots/fig_s13_lmer_residuals.pdf` |
| **Figure S11** | NLS model selection frequency + DW statistics | `plots/fig_s14_model_selection.pdf` |
| **Figure S12** | Cook's D values + 4-method sensitivity comparison | `plots/fig_s15_sensitivity_summary.pdf` |
| **Figure S13** | Outlier removal rationale: 4-panel visualization | `plots/fig_s16_outlier_rationale.pdf` |
| **Figure S14** | Temporal trajectories by outlier status (lnRR) | `plots/fig_s17_temporal_outlier_trajectories.pdf` |
| **Figure S15** | Raw data trajectories by outlier status | `plots/fig_s18_raw_trajectories_outlier_status.pdf` |

Note: SI renumbered figures but disk filenames retain original numbering.

### Resilience / Robustness Figures

The in-pipeline resilience suite also produces SI and internal robustness figures: `fig_heatwave_cascade`, `fig_heatwave_cascade_regression`, `fig_heatwave_replication`, `fig_resistance_recovery`, `fig_kelp_resilience_paired`, `fig_ecological_memory`, `fig_resilience_stability`, `fig_s_methods_multiverse`, and `fig_s_env_moderators`. The full standalone suite can additionally produce `fig_compound_disturbance`, `fig_s_effectiveness_predictors`, and the Eisaguirre reproduction panels when external raw-data/comparison inputs are available. See `docs/RESILIENCE_SYNTHESIS.md` for manuscript status.

### Tables

| Table | Description | Output File |
|-------|-------------|-------------|
| **Table 1** | Average density/biomass by taxa | `tables/average_responses.csv` (data repo) |
| **Table 2** | Meta-analysis summary statistics | `tables/table_02_meta_analysis.csv` |
| **Table 3** | Cross-taxa trophic cascade meta-regressions | `tables/table_03_cross_taxa_meta_regression.csv` |
| **Table S1** | MPA data availability + taxa coverage (merged) | `tables/table_s_mpa_data_availability.csv` |
| **Table S1b** | Sample sizes by taxa × response × source | `tables/table_s_sample_sizes.csv` |
| **Table S2** | Variance components with CIs | `tables/table_s_variance_components.csv` |
| **Table S3a/b** | Source random effect sensitivity | `tables/table_s_source_sensitivity_*.csv` |
| **Table S4** | Temporal meta-regression (pooled + bio + den) | `tables/table_s_temporal_meta_regression.csv` |
| **Table S5** | Cascade consistency scores by MPA | `tables/table_s_cascade_consistency.csv` |
| **Table S6** | Moderator meta-regression + diagnostics | `tables/table_s_moderator_meta_regression.csv` |
| **Table S7** | Temporal diagnostics: GAM linearity + AR1 | `tables/table_s_gam_linearity_diagnostics.csv` |
| **Table S8** | Kelp leave-one-out sensitivity | `tables/table_s_kelp_leave1out.csv` |
| **Table S9** | Outlier sensitivity (4 methods) | `tables/table_s_outlier_sensitivity.csv` |

Resilience outputs checked by `code/R/24_resilience_pipeline_check.R`:

| Output | Description | Script |
|--------|-------------|--------|
| `tables/table_resistance_recovery.csv` | Main Figure 5 resistance/recovery values for kelp and cascade taxa | 22 |
| `tables/table_resistance_recovery_pairs.csv` + `tables/table_resistance_recovery_pair_counts.csv` | Pair-level support for the 9/10 and 8/10 reserve-count statements | 22 |
| `tables/table_resistance_recovery_sensitivity.csv` | Baseline-window, paired-test, CI, and leave-one-reserve-out checks for the kelp headline | 22 |
| `tables/table_heatwave_replication.csv` | Two-heatwave repeatability support | 19 |
| `tables/table_s_methods_multiverse.csv` + `tables/table_s_methods_decomposition.csv` | Kumagai-style-to-ours method bridge and data-effect support | 15 |
| `tables/table_s_env_moderators.csv` | Tests of environmental covariates mapped by Kumagai et al. | 16 |
| `tables/table_resilience_stability.csv` + `tables/table_ecological_memory.csv` | Stability and reserve-level memory support | 21, 23 |
| `tables/table_s_pisco_control_source_removal.csv` + `tables/table_s_pisco_control_leave_pair_out.csv` | PISCO reference-control dependence sensitivity | 26 |

---

## Statistical Methods

### pBACIPS Methodology

We implement the **progressive-change Before-After-Control-Impact-Pairs** approach following [Thiault et al. (2017)](https://doi.org/10.1111/2041-210X.12655). This method:

1. Pairs each MPA site with a nearby reference site
2. Calculates the log response ratio (ln RR) between paired sites
3. Fits 4 candidate models and selects best via AICc:
   - **Step**: Abrupt shift after MPA implementation
   - **Linear**: Gradual constant-rate recovery
   - **Asymptotic**: Fast initial change that saturates
   - **Sigmoid**: Delayed onset followed by rapid change

### Meta-Analysis

Effect sizes are synthesized using multilevel meta-analysis (`metafor::rma.mv`) with:
- Taxa as fixed effect moderator
- Random effects for MPA and data source
- Heterogeneity assessment (I², τ²)
- Separate models for biomass and density responses

### Quality Control

- Before-period trend screening where applicable; a non-detected pretrend does not verify parallel untreated trajectories
- Bootstrap uncertainty quantification for biomass estimates
- Cross-validation of effect sizes across data sources

---

## Key Dependencies

| Package | Purpose |
|---------|---------|
| `tidyverse` | Data wrangling and visualization |
| `here` | Reproducible file paths |
| `metafor` | Meta-analysis models |
| `emmeans` | Estimated marginal means |
| `lme4` | Mixed-effects models |
| `investr` | Nonlinear model inference |
| `cowplot` | Figure composition |
| `patchwork` | Multi-panel figure assembly |
| `sf`, `terra`, `tidyterra`, `marmap`, `elevatr`, `rnaturalearthdata` | Geographic mapping |

---

## Documentation

- **`docs/ANALYSIS_REVISIONS.md`** - Comprehensive comparison of V5 archived code vs current pipeline
- **`docs/RESULTS_SUMMARY.md`** - Auto-generated results summary with all meta-analysis values
- **`docs/RESILIENCE_SYNTHESIS.md`** - Plain-language synthesis of scripts 14-23
- **`docs/RESILIENCE_PIPELINE_INTEGRATION.md`** - Paper/SI resilience claim-to-output audit trail
- **`docs/README.md`** - Public documentation map and notes on removed manuscript/PDF artifacts

Generated browser HTML, manuscript drafts, and source-PDF libraries are
intentionally not tracked in this public analysis repo. Regenerate paper-facing
manuscript products from the private manuscript repo.

### Resilience analysis suite

Scripts 14–23 form an integrated **resilience module**. The main pipeline runs
`14/15/16/19/21/22/23`, then `24_resilience_pipeline_check.R` verifies that the
generated tables and figures still match the manuscript and Supporting Information.
Run `code/R/run_resilience.R` for the broader standalone suite; it adds the external
input modules `17/18/20` when the raw PISCO or comparison inputs are available.

The suite tests resilience across facets: core heatwave response (14), repeatability
across two heatwaves (19), resistance/recovery on the state variable (22), ecological
memory/reserve-consistency (23), temporal stability (21), method-vs-data robustness
(15), environmental moderators (16), generality across stressor types / sea-star
wasting (20), predictability of effectiveness (18), and cross-study context (17). The modern diver-kelp advantage is an observed baseline-relative pattern, not a demonstrated effect of designation. Baseline, test, cohort and model choices matter. The historical satellite comparison and cross-taxa tests leave attribution unresolved. Disease and environmental-moderator checks do not remove confounding merely because they are non-significant. Full primary and sensitivity results are separated in `docs/SCIENTIFIC_APPROACH.md` and `docs/RESILIENCE_SYNTHESIS.md`.

### Related work (compare & contrast)

Kumagai et al. (2024) and Ortiz-Villa et al. (2025) provide context for California kelp responses to protection and marine heatwaves. Our paired population contrasts, survey-based disturbance responses and historical-event comparison address related but distinct questions. Differences in observations, footprints, periods, response definitions and analytical aggregation prevent attributing agreement or disagreement to one methodological choice.

The method/data bridge (module 15) and environmental-moderator analysis (module 16) are supporting comparisons, not evidence that the result is method-invariant or free of environmental confounding. Non-detection after multiple-test correction does not establish that an environmental gradient has no effect or that omitted covariates cannot bias a contrast.

The Eisaguirre comparison (module 17) examines design-level patterns and size-structured interpretations using available PISCO records; differences in data, collinearity and specification limit mechanistic attribution. It is not a new demonstration of a complete causal cascade.

The sunflower-star comparison (module 20) is small and observational. No detected moderation by pre-disease star abundance does not rule out disease effects. The second-heatwave comparison (module 19) has a consistent kelp direction, but does not clearly distinguish the second event from the intervening period; it is not an independent causal replication of a protection benefit.

Finally, we test whether the among-MPA variation in **effectiveness** is *predictable* (`code/R/18_mpa_effectiveness_predictors.R`, doc `docs/mpa_effectiveness_predictors.md`): assembling environmental, reserve-design, and trophic predictors per MPA and validating with out-of-sample checks. The current result is weak and method-dependent: no univariate predictor survives FDR, random-forest out-of-bag R² is negative for kelp and purple urchin, purple-urchin LOO-CV R² is strongly negative, but a small-k top-three linear LOO model gives kelp LOO-CV R² = 0.389. The defensible statement is therefore the average MPA response plus its meta-analytic heterogeneity, with cautious language about any attempt to predict which reserves are most effective.

---

## Known Limitations

Key methodological caveats documented throughout the pipeline (see `docs/ANALYSIS_REVISIONS.md` for full details):

- **Kelp biomass calibration**: Diver biomass uses the twelve monthly dry-foliar
  slopes in SBC 127.5 Macrocystis biomass relationships, model 1, applied before annual averaging. KFM
  retains site-year matching of size and abundance surveys, which are not
  necessarily from the same month. This is dry foliar mass, not whole-plant wet
  mass or the satellite canopy index. See `docs/KELP_FIXED_COEFFICIENT_REFIT_2026-10-29.md`.
- **Dependence**: Density and biomass derive from overlapping observations but are modelled separately. Shared references and reused cross-programme size-frequency pools can still induce dependence that MPA/source random effects do not fully resolve.
- **Cross-taxa attenuation bias**: Cross-taxa meta-regressions (Table 3) use estimated effect sizes as moderators, introducing errors-in-variables bias
- **Scaling and offsets**: Each paired history is divided by a common maximum, not a site total. Common scaling cancels in an uncorrected log ratio; the added zero offset affects comparisons near zero and does not remove habitat or survey differences.
- **Cross-programme biomass bootstrap**: Missing sizes are drawn from matched MPA/status/year/species pools assembled across programmes. Separate count collection does not ensure independent biomass estimates.
- **Climate/stressor attribution**: The core pBACIPS analysis is paired but observational; heatwave, environmental-moderator, and sea-star-wasting checks are now included in scripts 14–23, but ENSO and other overlapping regional drivers cannot be fully isolated
- **Near-zero references**: Some large contrasts reflect very low reference kelp. Weighting and omission checks describe their influence but do not establish suitable reference habitat or eliminate offset sensitivity.
- **Reproducibility**: Package versions captured in `renv.lock`; tracked harmonized CSVs are the public analysis inputs; raw-data archiving is separate from this analysis repo

---
## Citation

If you use this code or data, please cite:

> Stier, A., Donham, E., et al. (in preparation). Kelp resilience and trophic responses in a southern California marine protected-area network.

Working title for the current Journal of Applied Ecology draft: *Kelp resilience and trophic responses in a southern California marine protected-area network*.

---

## License

This project is shared for academic collaboration and code review. No formal
open-source reuse license has been assigned yet; please contact the authors
before reusing code, data products, or figures for other purposes.

The provider files in `data/external/sbc11215/` and `data/external/sbc1275/` retain their separate CC BY 4.0
license; see that directory's README for attribution and reuse terms.

---

## Acknowledgments

We thank the Partnership for Interdisciplinary Studies of Coastal Oceans (PISCO), the National Park Service Kelp Forest Monitoring Program, and the Santa Barbara Coastal LTER for data access. Funding acknowledgements will accompany the final publication.

---

<p align="center">
  <i>University of California, Santa Barbara</i><br>
  <i>Department of Ecology, Evolution, and Marine Biology</i>
</p>

Reference-selection history and remaining evidence gaps: [assignment audit](docs/REFERENCE_ASSIGNMENT_AUDIT.md), with one ledger row per primary history.
