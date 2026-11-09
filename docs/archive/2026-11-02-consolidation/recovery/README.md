# README recovery and cleanup audit

November 2, 2026. Both missing README working versions have been reconstructed and preserved as historical documents. The current consolidated guides remain in place. No further file loss was found within the checks below.

## Recovered documents

- [Analysis README before cleanup](../repository_README_before_cleanup.md): 527 lines; 13 recorded edits recovered.
- [Data README before cleanup](../../../../../sbc-donham-kelp-mpa-data/docs/archive/2026-11-02-consolidation/repository_README_before_cleanup.md): 100 lines; five recorded edits recovered.

The original failure was a destination collision: each root README was copied to the same path later used for its archive index. Inspection of the cleanup script's destination scheme identifies exactly those two collisions. The other 47 snapshots were not overwritten.

## How recovery was verified

The long-running project chat is stored under its original August 31 session date, not the date of the latest turn. Its records contain the edit commands and the pre-cleanup README readouts. We reconstructed the documents from pinned Git bases and the 18 successful README edits. A failed earlier command was excluded because it stopped before reaching its README write. No historical shell commands or unrelated edits were rerun.

The recovered analysis README matches 5,378 characters of independently displayed pre-cleanup text. The data README matches its 3,393-character pre-cleanup readout. The reconstruction preserves the remainder from the Git bases and recorded edits, including appended analysis notes. Replaying the saved recipe produces identical recovered files.

These are reconstructed working versions, not originals recovered from a disk snapshot. There is no independent whole-file checksum from immediately before cleanup, so whole-file byte identity cannot be independently certified. Time Machine reported no local snapshots. The earlier statement that the local edits were unavailable is superseded by this recovery.

## Other cleanup checks

- All 47 unaffected archive snapshots match their pre-recovery manifest checksums.
- All 49 archive entries now have verified files and unique destinations distinct from archive indexes. The original documentation paths still exist; all archive-pointer links resolve.
- All 487 protected analysis/manuscript files match their pre-cleanup protection manifest.
- All 75 current decision-evidence hashes match.
- All 25 registered raw-source files match their source-registry checksums.
- All 28 source-build products and five recorded configurations match their build manifest; the builder, source registry and two KFM selection/mapping files also match their recorded hashes.
- Neither the analysis nor data repository has deleted tracked files.
- Current guide navigation and the decision-table checker pass. No analyses, raw data, manuscript text or numerical outputs were changed during recovery.

There is no complete pre-cleanup checksum inventory of every ignored or untracked file. These checks support **no additional loss found**, rather than an absolute guarantee about every file on disk.

## Evidence and reproduction

- [Verification results](verification.json)
- [Recovered analysis edits](analysis_recovered_edits.diff) and [recovered data edits](data_recovered_edits.diff)
- [Recorded README transformations and intermediate hashes](reconstruction_steps.json)
- [Analysis pre-cleanup readout](analysis_pre_cleanup_excerpt.txt)
- [Data pre-cleanup readout](data_pre_cleanup_excerpt.txt)
- [Analysis manifest before recovery](analysis_manifest_before_recovery.csv)
- [Data manifest before recovery](data_manifest_before_recovery.csv)

From the analysis repository, run `python3 code/tools/check_readme_recovery.py`, `python3 code/tools/check_documentation_consolidation.py` and `python3 code/tools/check_current_decisions.py`. The README checker replays only the reviewed text transformations into temporary files and validates the archived results. It does not replace the current guides or run models.

The archive checker now rejects duplicate destinations or an archived document occupying the archive index/manifest path. Prior Git-only copies remain available as reconstruction bases; neither they nor the current guides were overwritten by this recovery.
