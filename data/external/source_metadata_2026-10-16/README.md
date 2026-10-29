# Source evidence recovered 16 October 2026

These are audit outputs, not revised primary inputs. The current primary cohort
and results remain unchanged. Start with
`docs/NO_TAKE_SATELLITE_SOURCE_AUDIT_2026-10-16.md`.

`no_take_review.csv` separates current CDFW ordinary-fishing rules from inherited
analysis inclusion. It is not Jenn's missing programme classification list and
does not verify every regulation in every study year. Scientific collecting,
tribal exemptions and other special permits are not encoded as ordinary fishing.
Site location and exposure to fishing still matter even where take is limited.

Recovered attachments and NetCDF remain in ignored
`local_archive/metadata_audit_2026-10-16/`, not this public evidence directory.
Gmail message IDs and Drive file ID identify the read-only original retrievals:

- 2020-07-02, Tom Bell, MPA Kelp Data: `1731158abac79364`,
  `MPA_Sites_Controls.xlsx`.
- 2020-12-19, Tom Bell, Re: prop 84 project: `1767de09ce2383a5`, `BACI.zip`.
- 2022-08-12, Tom Bell, Re: Kelp biomass inside/outside mpas:
  `182945f103c3bafb`, linked `CAkelpCanopyEnv_2021.nc`, Drive
  `1353POLzjcAZ5lrw9oHd8D9hLpSLe_AlW` (207,796,876 bytes).
- 2022-10-21, Tom Bell, Re: Follow up about LANDSAT/PISCO data:
  `183fa8c993c82b45`, `MPA_Ref_Locations.xlsx`.
- 2026-10-15, Jenn Caselle's manuscript attachment comment 24:
  the monitoring groups' de facto rocky-reef SMR list reportedly supplied to Emily.

Reproduce workbook retrieval and comparison with the bundled Python runtime:
`python3 code/tools/audit_source_metadata.py`. This requires authorized `gws`
Gmail read access. It reuses cached original workbooks if present.

Reproduce pixel/NetCDF checks with:
`Rscript code/tools/check_satellite_pixels.R /path/to/CAkelpCanopyEnv_2021.nc`.
The NetCDF can be recovered through authorized `gws drive files get` with the
observed file ID above, `alt=media` and a local output path. R requires ncdf4,
readxl and jsonlite. No EDI credential is required for this Drive retrieval.

The pixel audit sums each quarter over valid biomass pixels, then takes the
annual maximum. It also tests spatial means. It does not reproduce a minimum
valid-pixel fraction filter, fill missing years, infer polygon equivalence, or
replace delivered data. Exact agreement on paired finite years is reported
separately from agreement in missingness. Candidate reference labels may be
swapped between workbook and CSV; compare values and pixel sets, not rep alone.

The NetCDF's `Units` attribute (capital U) identifies wet kilograms per 900 m2
pixel for biomass, versus square metres of canopy per pixel for area. The actual
analysis derivative remains dimensionless, and wet canopy biomass must not be
called diver-derived dry foliar biomass.
