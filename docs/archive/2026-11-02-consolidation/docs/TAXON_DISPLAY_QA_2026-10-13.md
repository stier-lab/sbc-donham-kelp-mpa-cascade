# Sheephead display names

The current accepted display name is Bodianus pulcher (B. pulcher). Historical source names, factor levels, statistical keys and joins remain unchanged. The authority record is `data/taxon_display_names.csv`.

The supplementary plotting routes in scripts 10, 11 and 13 now use `taxon_display()` for discrete axes, facet labels and phase-portrait axis labels. The helper preserves R plotmath expressions while changing species strings, so italics and mathematical notation remain intact.

Eight affected supplementary assets were regenerated with the existing cached figure runner: `fig_s04`, `fig_s05`, `fig_s06`, `fig_s09`, `fig_s14`, `fig_s15`, `fig_s16` and `fig_s17`. All files in `tables/` remained byte-identical to the pre-render copies under `local_archive/pondella_taxon_display_2026-10-13/before_tables`. No main figure was selected for regeneration. No data inclusion, model specification or inference rule changed.

`Rscript code/tests/test_taxon_display.R` checks full and abbreviated labels, compound labels, plotmath expressions, unchanged purple-urchin labels, unchanged source factor levels and unchanged plotted numerical values. The cached runner completed successfully; it retained existing warnings from diagnostic plots. The updated sheephead labels were checked visually in all eight PNGs and against text extracted from the PDF figures referenced by the legacy supplement. This is a label-correction check, not a new statistical or comprehensive legacy-figure readiness audit.

The manuscript repository rebuilds the legacy supplement's generated tables and narrative display strings through the same helper. The current paired readability manuscript and core supplement remain the coauthor reading package.
