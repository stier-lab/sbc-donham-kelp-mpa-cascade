# Reference identity registry

These small, analysis-ready derivatives are the only reference-registry inputs to routine script 27. The ordinary pipeline never reads the original PISCO fish/swath files, site mappings or data-processing code. Rebuilds are an explicit one-time audit, not part of run_all.R.

- reference_sites.csv: one declared reference identity per primary series, including excluded sites. reference_years records surveyed PISCO site/year footprints intersected with the finite harmonized primary history. For MBON and Landsat, years describe harmonized coverage under declared membership, not reconstructed individual reference contributions. included, excluded, exclusion_reason, identity_verified, contribution_audit, coverage_kind and registry_version document these distinctions. Mapping, footprint and exclusion-definition MD5 columns record the source snapshot.
- reference_primary_keys.csv: every current primary series, including MPA, Taxa, Source, Resp, Region, AnalysisType, Model and exact observed/before/after year lists. Missing declarations remain explicit unmapped identities, not dropped series.
- mapping_provenance.csv: role, kind, path, md5. analysis_input rows bind the registry to both harmonized lnRR files and sumstats_final.csv. derived_registry rows checksum the vendored site/footprint and key files. source_audit rows retain original mapping, survey and exclusion/processing-code checksums for provenance only; routine analyses do not open these paths.

Script 27 stops before fitting if a required input or derived-registry MD5 differs, a current primary key/model/year classification is absent or changed, or reference-year flags/coverage are inconsistent. It cannot detect new changes in an unavailable upstream workspace: an upstream data refresh requires an explicit re-export, followed by updating the repository input-checksum manifest. No raw survey rows are vendored.

PISCO uses exact declared site IDs, swath campus/zone and May-October filters, and fish campus/non-CAN filters (no fish seasonal filter). Source exclusion constants are parsed individually; the upstream pipeline is never sourced. KFM exclusions apply to non-sheephead series; LTER uses its declared exclusion list. Exact shared IDs connect groups; coordinates or outcome similarity never create identities. Landsat MPA labels are not global polygon IDs, so sharing/overlap remains unverified. PISCO biomass-specific missingness and MBON individual site contributions are not reconstructed.

Registry rows: 861; primary series: 145.

```sh
Rscript code/R/export_reference_registry.R /path/to/Donham-Stier-CA-MPA-Data-2026
```
