# Source-first pipeline implementation

The data repository now starts from 19 checksum-verified delivery snapshots, rather than silently treating the combined size file as an original. It rebuilds KFM measurements from the historical individual files and Scott's workbook and recovers UCSB lobster sizes from the separately delivered swath surveys. Source rows retain file hashes, original row numbers and dates. Researcher exclusions are preserved.

[Full workflow, calculations, donor rules and limitations](../../../sbc-donham-kelp-mpa-data/docs/RAW_TO_BIOMASS.md)

## Executed and verified

- Source intake copied and SHA256-verified 19 original deliveries or explicitly labeled provider/curator products. The snapshots are read-only; original Drive files were not changed.
- Source reconstruction read 253,228 KFM red/purple measurement rows. The current approved group roster retains 158,559 rows representing 206,028 measured animals. This count includes no source-only additions.
- All selected KFM group animal totals agree with the independent source audit. Hare Rock uses 102 corroborated measurements. Pelican 1995 remains held.
- Reconstructed 2,109 focal UCSB lobster size rows from swath records, preserving multiplicity and converting cm to mm. They match the old combined subset on physical measurements/counts/calendar dates. Source survey-year labels are preserved separately.
- Quarantined 25,426 focal standalone PISCO/LTER size rows whose independent originals or merge recipe were not found. They are visible in the source report, not automatically added to ready inputs.
- Executed source-to-mass ledgers for PISCO and LTER sheephead, UCSB and LTER lobster, KFM urchin sizes, KFM kelp plants, PISCO kelp and LTER kelp. LTER lobster imputation includes 968 recipient–donor links across the original source file. These are source-level calculations before final analytical eligibility, not new manuscript estimates.
- Instrumented production animal resampling to retain exact donor source IDs, recipient context, conditional mass/SD, formula and iteration count. A fixed-seed test verifies that logging does not change results or random draws. Invalid donor/count states stop rather than disappear.
- Added hash-aware cache dependencies and forced regeneration for audited production runs, so legacy caches cannot bypass donor-use records.
- Seven Python source tests and R tests for lineage, source identity/permutation and missing-lobster-length behavior pass. All R modules parse. The main production entry point was run and stopped at the intended source-completeness guard before harmonization.
- All 487 protected primary files remain unchanged. No manuscript/supplement changes, email sends or push.

## Not yet complete

This is not an end-to-end newly harmonized primary dataset. Completion is blocked by standalone PISCO/LTER size-source documentation, Pelican 1995, and source coverage for the four retained but currently unused 2023 KFM groups. Explicit required-year gating could avoid unused source holds, but has not been silently substituted for the full-input guard.

The current pooled donor rule is documented; the recent visit-matching/joint-uncertainty alternatives are not silently promoted. Existing downstream zero-fill, sampling-unit and potential LTER urchin area-normalization behavior still needs numerical reconciliation before primary promotion. Those issues were identified by code review, not demonstrated new primary-result changes.

## Reproduce

In the data repository:

```sh
make sources
make trace
make check
make harmonize
```

The first three succeed with the available snapshots. The last deliberately reports the source holds. `make trace` creates a fresh run directory under `output/lineage/`; this run completed at `20261102T141017`. Complete readme and source/mass tables are local/private. The ledger columns expose original versus estimated measurements and retain unsupported positive observations as missing mass, distinct from observed zeros.

Independent code review identified and repaired stale-trace validation, absent approved-group readiness, remaining mutable input paths, fractional donor indices, misleading zero/missing-count donor links and unrestricted receipt metadata. Targeted gws searches found no independent PISCO/LTER size attachment or construction recipe in the checked correspondence; this is not proof that those originals never existed. No emails were sent.
