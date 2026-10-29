# Source reconciliation diagnostics

See `docs/SOURCE_RECONCILIATION_2026-10-17.md` for conclusions, source links,
rebuild commands and limitations. These files do not replace primary inputs.

- `report_naples_verified.csv`: two manually transcribed, visually verified
  published coefficients, with PDF and printed page numbers.
- `naples_report_style_diagnostics.csv`: calendar-year regressions on current
  corrected PISCO inputs through 2020; not an exact report reproduction.
- `satellite_missingness_diagnostics.csv`: every same-MPA/status candidate
  comparison between workbook pixel groups and delivered replicates.
- `satellite_missing_years.csv`: missing delivered years that are reconstructable
  for groups whose paired finite annual values match exactly.
- `satellite_reference_centroid_geometry.csv`: reference-centroid and pixel
  distances from the mean coordinate of the matching MPA's workbook pixels.
- `survey_source_counts.csv`: row counts by file kind, programme and method.
- `survey_nearby_date_candidates.csv`: site-coordinate pairs within 1 km and
  overlapping valid dates. Not verified duplicate observations. Duplicate
  site-to-MPA mappings are collapsed for coordinate lookup only.
- `survey_audit_manifest.json`: source CSV and recovered report hashes, scope,
  distance threshold, source URL and caveats. CSVs are decoded as Latin-1
  because the supplied PISCO files contain non-UTF-8 notes; audited keys are
  ASCII. No malformed lines are skipped.

Workbook/NetCDF identifiers and their earlier integrity checks are retained
under `../source_metadata_2026-10-16/`. Original restricted input tables and
the downloaded report are not redistributed in this directory.
