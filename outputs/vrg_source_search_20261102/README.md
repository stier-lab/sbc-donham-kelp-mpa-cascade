# VRG size records: archived files and remaining source gap

November 2, 2026. Updated with Adrian’s acceptance of the contributor-combined VRG delivery. No measured values, ecological exclusions, primary results or paper files changed.

## What we could reconcile

The project’s combined file contains 18,750 focal VRG size-frequency records representing 28,487 animals. Both copies in Emily’s Drive project—under `data` and `Github/MPA_Kelp/data`—are byte-identical to our combined file. This verifies that we preserved her delivered combined table; it does not independently verify the original VRG measurements or the steps used to assemble it.

Emily’s hidden saved R workspace contains 18,284 VRG rows in `SizeFreq.south`. Every one matches the combined file exactly, including repeated rows. The remaining 466 records are entirely at three reference sites she explicitly excluded in her code following researcher advice:

| Physical site | Records removed |
|---|---:|
| SBI_GRAVEYARD_CANYON | 97 |
| SBI_SUTIL | 283 |
| SCAI_SHIP_ROCK | 86 |

These exclusions remain in place. All 18,750 focal rows have grouping labels consistent with the curated physical-site table; this check found no VRG counterpart to the Painted Cave alias problem.

## What the archive establishes about assembly

The visible analysis script starts by reading `ALL_sizefreq_2024.csv`. It then selects taxa, excludes reference sites, standardizes species labels, assigns comparison memberships and draws sizes for biomass. Its saved RStudio source copy follows the same downstream workflow. Neither provides the earlier import, unit conversion or row-binding steps that assembled the combined file.

The workspace contains derived size tables and biomass objects. Those are useful checks of filtering and use, but are not independent source deliveries. The temporary workspace `.RDataTmp` cannot be read as a complete R save file. It was left untouched. The available command history and RStudio source buffers did not reveal an earlier VRG import or assembly recipe.

The archived code describes VRG lobster lengths as collected from 2011 onward on transects and roving dives, with sizes separate from transect counts. Emily’s notes likewise distinguish these measurements from UCSB swath lengths. The delivered VRG swath file has no recorded sizes for the focal lobsters or adult red/purple urchins, so it cannot reconstruct the separate size-frequency table.

## Search scope

We inspected the local Emily Donham project tree, including hidden RStudio files, both data directories, the GitHub scripts and saved workspace. Live Drive filename searches included sizefreq, VRG, Vantuna and broader size/frequency CSV/Excel/archive candidates; the broader search completed two pages (148 files). The separate historical urchin-size files found in Bart DiFiore’s archived project contain UCSB data, not VRG. They do not fill this gap. This is a search of the accessible material, not proof that the original VRG delivery no longer exists.

## What remains unresolved

We still need the original VRG size-frequency delivery that supplied these records, its units and missing-value definitions, and any script or spreadsheet steps used to combine it with the other programs. The record ranges are red urchins 2004–2023, purple urchins 2007–2023 and lobsters 2011–2023.

The practical follow-up is to ask Emily for that original delivery or its source email, plus any assembly script. If she no longer has it, the original VRG contributor may be able to provide a matching export and metadata. No email was sent during this audit.

## Pipeline decision — accepted by Adrian

Adrian explicitly approved using Emily’s combined file on November 2, 2026. The pipeline now accepts its 18,750 focal VRG records as a contributor-combined delivery, retaining their file hash and row identities. It removes the same 466 researcher-excluded records and passes the remaining 18,284 to candidate size inputs. Existing comparison exclusions remain in force. This clears the VRG source hold; no further original-file request is required to proceed under this decision.

The data are labeled `user_accepted_contributor_combined`, not independently verified original measurements. Millimeters retain the interpretation used by Emily’s biomass workflow. Acceptance is limited to the VRG subset and exact file checksum; it does not restore a fallback for other programs. The original-delivery/assembly gap remains documented as a limitation rather than a blocker. Separate KFM source holds and biomass-matching/uncertainty work remain. No new model fit or paper change was made.

[Recorded acceptance and scope](acceptance_decision.json)

## Tables and reproduction

- [Coverage by species](vrg_coverage.csv)
- [Reconciled workspace exclusions](workspace_exclusions.csv)
- [File hashes and evidence roles](source_inventory.csv)
- [Workspace object inventory](workspace_inventory.csv)
- [Verification summary](verification.json)

Reproduce with `Rscript --vanilla code/tools/inspect_vrg_workspace.R`, followed by `python3 code/tools/reconcile_vrg_archive.py`. Workspace deserialization loads saved model namespaces but does not run the archived analysis scripts; unavailable package warnings are recorded during inspection. Source files are not rewritten.
