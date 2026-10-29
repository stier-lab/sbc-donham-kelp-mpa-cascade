# Figure Visual Audit

Date: 2026-09-29

Scope: final visual render/check requested for the information-design pass:
figures should not carry chart-level titles or subtitles inside the artwork.

## What Changed

- Removed `ggtitle()` calls from manuscript and resilience figure source.
- Removed `labs(title = ...)` and `labs(subtitle = ...)` calls from figure source.
- Removed titled `plot_annotation()` calls from patchwork figures.
- Retained necessary interpretive labels: panel tags, axis labels, legends, map
  labels, direct labels, category headers, and in-panel event labels.
- Replaced the Figure 1 coastline path with `rnaturalearthdata::states50` so the
  map no longer tries to install `rnaturalearthhires` at render time.

## Verification

Commands/checks run:

```sh
Rscript -e 'source(here::here("code", "R", "run_all.R"))'
rg -n "ggtitle\\(|labs\\([^\\n]*(title|subtitle)|plot_annotation\\([^\\n]*(title|subtitle)" code/R
```

The final full run completed successfully:

- Log: `logs/pipeline_20260929_191058.log`
- Elapsed time: 2.89 minutes
- Core manuscript figures regenerated: `plots/fig_01_mpa_map.pdf`,
  `plots/fig_02_cascade_case_studies.pdf`, `plots/fig_03_mean_effects.pdf`,
  `plots/fig_04_recovery_curves.pdf`
- Main-text resilience figure regenerated: `plots/fig_kelp_resilience.pdf`
- Resilience gate passed: `24_resilience_pipeline_check.R`
- PISCO sensitivity passed: `26_pisco_control_sensitivity.R`

Main-figure PDF text was also checked for removed title/subtitle phrases such as
"MPAs increase", "Estimated MPA effect", "Raw kelp biomass", "Standardized",
"Resistance & recovery", and "Giant-kelp resilience in each reserve"; no matches
were found.

## Visual Read

The five main figure PNGs were inspected after rendering. None has a
chart-level title or subtitle block. The remaining large labels are functional
figure content: site names, trophic group labels, map/site labels, the marine
heatwave event label, axes, legends, and panel tags.

## Remaining Polish

The pipeline still emits nonfatal warnings. Those warnings did not stop figure
generation or the resilience concordance gate, but they should be classified as
expected modeling/plotting warnings or cleaned before formal public code review.
