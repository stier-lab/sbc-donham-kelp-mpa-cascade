# Craig pre-review - fixes for Donham-Stier CA MPA analysis repo

(Simulated Craig pre-review - for internal cleanup before external review. Not
his actual words.)

## Fix first

1. Do not let the repo tell two different stories about authorship. The README
   and orientation file now match the current coauthor-review package: Adrian is
   lead/corresponding author and Emily is second author.
2. Do not call `run_figures_only.R` the all-figure pipeline. That was misleading
   once Figure 5 and the resilience/support checks moved into `run_all.R`. The
   docs now describe it as cached core figure regeneration.
3. Do not leave Figure 1 dependent on a hidden runtime install of
   `rnaturalearthhires`. The map now uses `rnaturalearthdata::states50`, and the
   full pipeline renders Figure 1 from the current code.
4. Do not bury the PISCO reference-control check. Script `26` now runs in
   `run_all.R`, is listed in the README, and has explicit outputs.
5. Do not put chart captions inside the figures. The source no longer assigns
   chart-level titles or subtitles; panel labels and data labels remain.

## Jargon & Plain Language

- "Resilience" is now grounded in the pipeline: Figure 5 comes from script `22`,
  and script `24` checks the Figure 5 values and supporting SI products.
- Kumagai is now framed as a comparator/context analysis, not as a claim that we
  merely reran their work.
- The PISCO-control sensitivity language distinguishes "dependence on PISCO
  reference controls" from the stronger claim of proving a reference reef is or
  is not well matched.

## Figures

- Main figures regenerated in the final full run:
  `fig_01_mpa_map`, `fig_02_cascade_case_studies`, `fig_03_mean_effects`,
  `fig_04_recovery_curves`, and `fig_kelp_resilience`.
- Visual check: no chart-level title/subtitle blocks remain.
- Direct labels retained where they carry information: site names, trophic group
  labels, map labels, the MHW event label, axes, legends, and panel tags.

## Reproducibility

- Local full run passed: `logs/pipeline_20260929_191058.log`.
- Resilience integration gate passed: `24_resilience_pipeline_check.R`.
- PISCO reference-control sensitivity passed: `26_pisco_control_sensitivity.R`.
- Remaining concern: `renv.lock` exists, but this checkout has no activated
  project renv library. Before external code review, run a clean clone plus
  `renv::restore()` and document the result.

## Statistics & Claims

- Keep the central claim narrow: MPAs are associated with stronger giant-kelp
  resistance and recovery during/after the 2014-2016 heatwave.
- Keep the mechanism language cautious: predator recovery and urchin patterns
  support partial trophic-cascade recovery, but the cross-taxa links remain weak
  or specification-sensitive.
- Keep the predictor/effectiveness language cautious: among-reserve prediction
  remains weak and method-dependent.

## Still Worth Fixing Before Broad Release

1. Capture and classify the "50 or more warnings" from a full run. If they are
   expected model/plot warnings, suppress or document them at source.
2. Run a fresh public clone after these uncommitted changes are committed and
   pushed.
3. Decide on license/DOI/archive policy before presenting the repo as a public
   reusable artifact.
4. If old manuscript/PDF artifacts must never have been public, rewrite history
   or create a fresh public repository; current-head cleanup alone is not enough.
