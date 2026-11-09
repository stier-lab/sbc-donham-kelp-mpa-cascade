# Original files, older records and pipeline safeguards

## What Emily's folder contains

Inventoried 86 local files in the specified Kelp MPA folder, excluding R-project internals. Inspected both combined-size CSVs, original V10 code and R history, code notes, conceptual notes, site-year workbook, site/geolocation tables and Drive file/revision metadata. The saved R workspace was also inventoried (576 objects); its size objects are already processed and do not supply the missing assembly script. No original natural-habitat export or combined-file assembly script was found there. The two size CSVs have identical bytes and match our current input.

Drive lists an original-project copy created January22,2026 with a retained revision carrying a July26,2024 modification date; its checksum matches today's file. The working copy was uploaded February2,2026 and has the same checksum. These records strengthen the evidence that the bad combined input predates the June code change. A preserved modification date is not proof of the exact assembly date or its author. Only one revision is exposed for the original-project copy, so no earlier clean version was recovered.

Emily's code notes explicitly say KFM conducts separate size-frequency surveys and that these sizes also support PISCO biomass estimates. That supports the intended cross-program pooling, but does not validate the erroneous site labels. `Old/Sites_Years.xlsx` records site/source/status and annual coverage, not the physical-site conversion that assembled sizes. Source metadata and original files remain unchanged.

## Older-record audit

Compared 904 delivered KFM red/purple site-species-year distributions through2011 against the legacy original measurements, using exact millimeter frequencies:

- 873 match their intended physical site.
- 15 match the already-known Wyckoff plus Miracle Mile mixture.
- 14 match the already-known Southeast Reef substitution.
- Two remain inconsistent: Hare Rock purple1987 has103animals versus102in the legacy file, with differences at22,23,28,29and30mm; Pelican Bay red1995 has135versus140, missing one each at51,63,105,134and184mm relative to the legacy source. Neither matches another site's full distribution. These may reflect source revisions, transcription or filtering; the evidence does not choose among them.

The recent audit supplies558further verified groups. Twenty2023site/species groups are outside Scott's2012–2022 export and are unverified, not demonstrated errors. The full registry contains1482groups,1460verified and22unverified. Existing researcher exclusions remove16of those unverified2023groups from required donors. Six site-retained raw groups still block a complete guarded production build: the two old discrepancies plus red/purple2023at Scorpion and Arch Point. No researcher exclusion was relaxed. The two old groups are used in the existing expanded review pools (Harris Point purple and Scorpion red); the four2023groups have no uses in that review pool. Thus six raw-data holds do not imply six affected current estimates.

## Pipeline changes

In the data-processing repository, `size_identity_utils.R` now validates the raw-source checksum, source-verified site/year/species registry and replacement-file checksums. Missing evidence, altered inputs, duplicate identity keys or unsupported identities stop processing. Module03 applies validated allowed-site replacements before pooling. The local original-measurement bundle is git-ignored; raw files are untouched.

Module04 replaces the first-match lookup with explicit physical-site/comparison memberships for KFM/LTER. Each source measurement appears at most once per comparison, with its source ID retained across comparisons. Conflicting roles, missing memberships and duplicate within-comparison use stop processing. PISCO's existing metadata path is retained with duplicate/conflict guards. Researcher site exclusions remain unchanged.

The three shared Anacapa reference sites now carry all three documented memberships, independent of metadata row order. The diagnostic verified subset has9959size records at those three sites (5649Admiral's Reef,1604East Fish Camp,2706Lighthouse); each is represented once in each mapped comparison, not three times in one pool. This does not reinclude excluded protected sites or promote Anacapa1978 into the paper. Table: [Anacapa membership coverage](anacapa_memberships.csv).

Urchin bootstrap cache names now incorporate the source/registry/replacement, geolocation, PISCO site-table and code/constants fingerprints, preventing old identity-dependent caches from bypassing the correction.

## What changes in results

No new manuscript or primary numerical results were produced. The previous six-history [corrected biomass results](../kfm_scott_reconciliation_20261102/index.html) remain the latest complete review of the known alias corrections. Here, the real-source guard correctly stops the full build at six unresolved retained groups. A clearly labeled diagnostic subset, excluding those unresolved groups only for testing, verifies the new membership implementation; it is not a new inclusion decision or result set. Resolving these six source groups or explicitly limiting the build to the intended surveyed analysis years is required before a full rebuild. Missing size evidence is not turned into zero, and density source files remain unchanged.

## Verification and reproduction

Run from the analysis repository:

- `python3 code/tools/audit_older_kfm_identity.py`
- `python3 code/tools/build_kfm_identity_bundle.py` (builds and installs the private checked bundle)
- `Rscript code/tools/check_kfm_pipeline_guard.R`

Run `Rscript code/tests/test_size_identity.R` from the data repository. Tests cover altered sources, invalid/duplicate registry keys, missing site membership, conflicting comparison roles, duplicated source measurements, and invariance to metadata row ordering. All five edited R modules parse. The real complete input raises the expected unverified-source error; the diagnostic covered subset passes with no excluded donor IDs. Existing primary and manuscript outputs remain unchanged. No email or push.

## Correction to the earlier Git-origin interpretation

The saved workspace and master site table show that some KFM sizes already had MPA metadata before a backfill. In particular, Wyckoff’s source key matches the existing Harris Point reference assignment. The June commit is not evidence that all KFM sizes were previously unused. Treat its broad commit-message explanation as an unverified historical claim; trace each source key separately. The known physical-size mixing remains proven by exact source frequencies, independent of this timing correction.

The saved size pool specifically contains13964frequency-counted animals under Wyckoff/Harris Point and4844under the Southeast Reef source key with unmatched MPA metadata. Thus the two errors had different downstream paths: Wyckoff already matched the original site table; the mislabeled Southeast Reef key required the later site-ID backfill. Saved-state totals describe that processed workspace, not the full original export or the final analysis sample.
