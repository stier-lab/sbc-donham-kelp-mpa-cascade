# Landsat annual index input

This tracked analysis-ready derivative retains observed zeros in the delivered annual Landsat canopy product. Six verified empty second-reference slots are removed before aggregation; these are not zero-biomass observations. It contains one value per MPA, status and complete paired calendar year. Remaining replicate annual values are averaged before scaling by the common maximum across both statuses within each MPA, using the data-processing repository's existing calculate_proportions() implementation.
The export checks that index + log_offset reproduces every mpa/reference value in harmonized_landsat_rr.csv to 1e-10. It does not reclassify sites, change polygons, infer physical units, or extend the record beyond 2021. The historical-event analysis uses window means of the index before adding the zero offset. Its dimensionless ratios are invariant to the common scaling factor.

Source filename: `MPA_Runs_new.csv`
Source MD5: `503ee8bc832a8abef46cb0bbe8d9d895`
Data-processing utility MD5: `896cc42a328364ec7206e98d965f1760`
Export rows: 1436; MPA labels: 20; observed years: 1984-2021. The source has 60 rows, of which 54 are retained. The pinned release marker is `landsat_reference_correction.json`; it records the exact six slots, source checksum and both corrected-file fingerprints. Diver inputs were unchanged.

Rebuild explicitly from the original source workspace:

```sh
Rscript code/R/export_landsat_annual_index.R /path/to/Donham-Stier-CA-MPA-Data-2026
```

The explicit rebuild regenerates both the annual index and the ratio handoff through the canonical corrected data-processing exporter. The ordinary analysis pipeline reads the tracked derivatives and validates the correction marker; it needs neither the original source file nor the archived data-processing workspace.

Any changed primary inputs or fitted histories also invalidate the pinned reference
registry. After the primary refit, explicitly run
`Rscript code/R/export_reference_registry.R /path/to/Donham-Stier-CA-MPA-Data-2026`
and rerun the full pipeline. Script 27 must reject a stale registry rather than
silently use footprints from the previous release.
