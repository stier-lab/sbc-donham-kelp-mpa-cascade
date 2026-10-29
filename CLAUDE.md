# CA MPA Kelp Forest Analysis - Assistant Notes

This file is a public-safe guide for coding assistants working in this analysis
repository. Human collaborators should start with `README.md`.

## Repository Boundary

This repository is the analysis source of truth for the Donham-Stier California
MPA kelp-forest manuscript. Keep it limited to:

- R analysis code under `code/R/`
- tracked analysis-ready inputs under `data/harmonized/`
- small repo-controlled comparator inputs under `data/external/`
- generated tables, figures, and analysis documentation

Do not add Word manuscripts, cover letters, generated manuscript/SI HTML, or
publisher PDF libraries. Those belong in the private manuscript workspace.

## Entry Points

Use these commands from the repository root:

```r
source(here::here("code", "R", "run_all.R"))
source(here::here("code", "R", "run_figures_only.R"))
source(here::here("code", "R", "run_resilience.R"))
```

`run_all.R` is the manuscript pipeline. It runs the core pBACIPS/meta-analysis
workflow, the in-pipeline resilience modules, and
`24_resilience_pipeline_check.R`.

`run_figures_only.R` uses `data/cache/figures_snapshot.rds` to regenerate figures
without rerunning the full analysis.

`run_resilience.R` runs the broader resilience suite. Scripts `17` and `20`
require optional raw PISCO inputs; they should skip cleanly when those inputs are
absent.

## Data Inputs

The 18 October 2026 singleton-resampling correction follows the monthly-kelp
and sheephead releases. `lobster_bootstrap_correction.json` and
`urchin_bootstrap_correction.json` pin the scoped input changes. Sampling must
use vector positions: `a[sample.int(length(a), n, replace = TRUE)]`, never
`sample(a, n, ...)` for a possible numeric singleton. The corrected primary
cohort still has 145 histories; five biomass contrasts change. The subsequent SBC127.5 update retains all nine pooled directions but adds a
negative purple-urchin biomass result, sensitive to selection of kelp contrast
forms rather than the tested HAC bandwidths. The four-MPA lobster/urchin biomass association has
a nominal negative interval but exact permutation P = 0.125. See
`docs/BIOMASS_SINGLETON_REFIT_2026-10-18.md`. Do not restore stale bootstrap
caches or claim every annual-uncertainty conclusion is unchanged.

SBC 127.5 Macrocystis biomass relationships, model 1, is the adopted monthly dry foliar calibration.
Harmonization applies recorded-month slopes before annual averaging; KFM
resamples whole dated size records and retains site-year abundance matching.
Do not replace the lookup with an unweighted mean or claim that monthly
conversion preserves annual lnRRs. `33_kelp_calibration_check.R` verifies the
adopted input fingerprints before loading. Source bytes are pinned in
`data/external/sbc1275/`; see `docs/KELP_FIXED_COEFFICIENT_REFIT_2026-10-29.md`.
The older SBC112.15 bytes remain only as a verified ancestor in the release chain.

PISCO and LTER sheephead biomass use Hamilton & Caselle (2015),
`W(g) = 0.0144 * TL(cm)^3.04`, multiplied by size-class abundance.
The provider's original PISCO lookup is preserved; its SPUL coefficients are
overridden during source-data processing. `data/harmonized/sheephead_allometry.json`
records the sourced correction and its lineage from the monthly kelp release.
`35_sheephead_allometry_check.R` blocks stale inputs, and the figures-only runner
also requires the cache to match current input fingerprints.

Adrian approved retaining Naples LTER sheephead biomass on 12 October 2026.
The primary set has 145 histories (108 CI, 37 pBACIPS). The added history has
23 paired years and retains the selected linear model; raw inputs and the
other 144 effects are unchanged. Module 34 exports its former exclusion as
a joint-model sensitivity. Do not restore the old omission or change Naples
density as part of this decision. See `docs/NAPLES_BIOMASS_INCLUSION_2026-10-12.md`.

The main pipeline is self-contained for code review:

- `data/harmonized/*.csv` are the analysis-ready inputs.
- `data/cache/*.rds` are tracked so figure-only regeneration works.
- `data/external/kumagai2024/` contains the small Kumagai comparator files needed
  by scripts `15`, `16`, `18`, and `24`.

Optional external inputs can be overridden with environment variables:

- `KUMAGAI_SUBTIDAL_CSV`
- `KUMAGAI_CS_GRID`
- `KUMAGAI_HUMAN_GRAVITY`
- `DONHAM_DATA_REPO_PISCO`
- `SBC_KELP_ENV_NETCDF`

Do not hard-code user-specific absolute paths in new code.

## Active Pipeline Map

- `00_libraries.R`: package loading
- `00b_color_palette.R`: figure theme and colors
- `00c_analysis_constants.R`: shared constants
- `01_utils.R`: helper functions and data loaders
- `02_pBACIPS_function.R`: pBACIPS model fitting
- `03_load_harmonized_data.R`: tracked harmonized-data import
- `35_sheephead_allometry_check.R`: sourced sheephead formula/units and release checks before import
- `08_effect_sizes.R`: per-MPA effect sizes
- `09_meta_analysis.R`: multilevel meta-analysis
- `10_temporal_analysis.R`: retained-diver trajectories; nonsingular candidate hierarchy and exact-contrast Satterthwaite predictions
- `36_primary_uncertainty.R`: exact 145-history calendar-HAC covariance check, isolated after 09
- `37_temporal_uncertainty.R`: same-row biomass calendar-AR1 sensitivity, isolated after 10
- `11_figures.R`: main figures and most supplement figures
- `13_additional_analyses.R`: reserve-design moderator checks
- `14/15/16/19/21/22/23`: in-pipeline resilience modules
- `24_resilience_pipeline_check.R`: manuscript/SI concordance gate
- `12_results_summary.R`: generated result summaries

Module membership for the resilience suite is defined in
`code/R/resilience_modules.R`. Update that file first if the suite changes.

## Resilience Guardrails

The main manuscript uses one numbered resilience display: Figure 5 from
`22_resistance_recovery.R`. The broader suite supports interpretation and
robustness:

- `14`: AR1 heatwave-window model
- `15`: method-vs-data comparison against Kumagai-style endpoints
- `16`: environmental moderators from the Kumagai PCA covariate set
- `19`: repeatability across two heatwaves
- `20`: SSWD/sunflower-star alternative mechanism
- `21`: temporal stability
- `22`: resistance/recovery state-variable analysis
- `23`: reserve-level ecological memory
- `18`: among-MPA effectiveness prediction
- `17`: Eisaguirre paired-design reproduction

Use `docs/RESILIENCE_SYNTHESIS.md`,
`docs/RESILIENCE_PIPELINE_INTEGRATION.md`, and
`docs/CONCORDANCE_AUDIT.md` as the current prose-to-output map.

## Hygiene Rules

- Prefer `here::here()` or repo-relative paths.
- Use environment variables for optional external sources.
- Keep historical code and scratch material out of the tracked public head.
- If manuscript claims change, update the relevant table/figure, the prose-facing
  docs, and `24_resilience_pipeline_check.R` in the same commit.
