# Resolution of the 21 unmatched site assignments

November 1, 2026. Review only; paper, source data, models, reference choices and previous inventory remain unchanged.

## Main finding

Nineteen assignment rows (17 distinct group labels) are explicit grouped-site identifiers in the source master table. They link through `site_new` to the east/west/center/north/south station codes actually recorded in the fish and swath files. These are one-to-many group links, not interchangeable names for one physical station. They should be navigation/grouping entries, not additional sampled stations or missing datasets.

The previous audit checked direct station IDs and legacy-name fields but did not follow the master table's explicit `site_new` hierarchy. That made these 19 assignments appear unresolved. No observations were recovered or added by resolving the hierarchy. All linked PISCO constituent stations were already present in the inventory.

The other two labels, SCI_TRAILER and SCI_WILLOWS, are unassigned map entries: no MPA, status or coordinates, no constituent station mapping, and no exact or name-fragment records in the supplied fish, swath or size sources. They are unresolved historical/map labels, not confirmed ecological exclusions and not evidence that surveys never occurred.

## Identity and exclusion are separate

Of the 19 mapped assignments, 14 have no whole-group PISCO exclusion identified; 5 are associated with real exclusions: Anacapa West Isle in the SMCA; East Fish Camp and Lighthouse Reef in the 1978 reserve assignment; Yellowbanks; and Tyler Bight. East Fish Camp and Lighthouse Reef each also appear in the 2003 assignment, where the same reference data have a different comparison role. Counts therefore refer to 21 assignments, not 21 independent sites.

Yellowbanks CEN/W and Tyler Bight E/W are explicitly excluded in source code. Yellowbanks E has observed records in only two years and is not in that explicit exclusion vector; the ecological advice still does not justify restoring it as a local Gull reference. Anacapa SMCA and the 1978 reserve are whole-comparison exclusions, not missing records. Finding source data does not reverse any exclusion.

KFM codes also map to many of these geographic group names. These remain a separate network. KFM exclusions at Cavern Point, Coche Point, Little Scorpion, Potato Pasture and San Pedro Point must not be transferred automatically to PISCO constituent stations. Jenn's October 19 reply explicitly distinguishes separate KFM and PISCO observations even at nearby sites/dates. Geographic links do not authorize pooling or deduplicating across programs. Shared size samples remain a separate, previously verified issue.

## Every assignment

| MPA | assigned_group | raw_children | analysis_status |
| --- | --- | --- | --- |
| Anacapa Island SMCA | ANACAPA_WEST_ISLE | ANACAPA_WEST_ISLE_CEN;ANACAPA_WEST_ISLE_E;ANACAPA_WEST_ISLE_W | Real whole-MPA exclusion |
| Anacapa Island SMR 1978 | ANACAPA_EAST_FISH_CAMP | ANACAPA_EAST_FISH_CAMP_CEN;ANACAPA_EAST_FISH_CAMP_E;ANACAPA_EAST_FISH_CAMP_W | Real whole-MPA exclusion |
| Anacapa Island SMR 1978 | ANACAPA_LIGHTHOUSE_REEF | ANACAPA_LIGHTHOUSE_REEF_CEN;ANACAPA_LIGHTHOUSE_REEF_E;ANACAPA_LIGHTHOUSE_REEF_W | Real whole-MPA exclusion |
| Anacapa Island SMR 2003 | ANACAPA_EAST_FISH_CAMP | ANACAPA_EAST_FISH_CAMP_CEN;ANACAPA_EAST_FISH_CAMP_E;ANACAPA_EAST_FISH_CAMP_W | No whole-group exclusion found for this PISCO assignment |
| Anacapa Island SMR 2003 | ANACAPA_LIGHTHOUSE_REEF | ANACAPA_LIGHTHOUSE_REEF_CEN;ANACAPA_LIGHTHOUSE_REEF_E;ANACAPA_LIGHTHOUSE_REEF_W | No whole-group exclusion found for this PISCO assignment |
| Gull Island SMR | SCI_GULL_ISLE | SCI_GULL_ISLE_E;SCI_GULL_ISLE_W | No whole-group exclusion found for this PISCO assignment |
| Gull Island SMR | SCI_YELLOWBANKS | SCI_YELLOWBANKS_CEN;SCI_YELLOWBANKS_E;SCI_YELLOWBANKS_W | Real reference exclusion / short remaining branch |
| Harris Point SMR | SMI_HARRIS_PT_RESERVE | SMI_HARRIS_PT_RESERVE_E;SMI_HARRIS_PT_RESERVE_W | No whole-group exclusion found for this PISCO assignment |
| Harris Point SMR | SMI_TYLER_BIGHT | SMI_TYLER_BIGHT_E;SMI_TYLER_BIGHT_W | Real reference exclusion |
| Unassigned | SCI_TRAILER |  | Not an eligible study assignment; reason for map entry unresolved |
| Unassigned | SCI_WILLOWS |  | Not an eligible study assignment; reason for map entry unresolved |
| Scorpion SMR | SCI_CAVERN_POINT | SCI_CAVERN_POINT_E;SCI_CAVERN_POINT_W | No whole-group exclusion found for this PISCO assignment |
| Scorpion SMR | SCI_COCHE_POINT | SCI_COCHE_POINT_E;SCI_COCHE_POINT_W | No whole-group exclusion found for this PISCO assignment |
| Scorpion SMR | SCI_LITTLE_SCORPION | SCI_LITTLE_SCORPION_E;SCI_LITTLE_SCORPION_W | No whole-group exclusion found for this PISCO assignment |
| Scorpion SMR | SCI_PELICAN | SCI_PELICAN_CEN;SCI_PELICAN_E;SCI_PELICAN_FAR_WEST;SCI_PELICAN_W | No whole-group exclusion found for this PISCO assignment |
| Scorpion SMR | SCI_POTATO_PASTURE | SCI_POTATO_PASTURE_E;SCI_POTATO_PASTURE_W | No whole-group exclusion found for this PISCO assignment |
| Scorpion SMR | SCI_SAN_PEDRO_POINT | SCI_SAN_PEDRO_POINT_E;SCI_SAN_PEDRO_POINT_W | No whole-group exclusion found for this PISCO assignment |
| South Point SMR | SRI_CHICKASAW | SRI_CHICKASAW_E;SRI_CHICKASAW_W | No whole-group exclusion found for this PISCO assignment |
| South Point SMR | SRI_CLUSTER_POINT | SRI_CLUSTER_POINT_N;SRI_CLUSTER_POINT_S | No whole-group exclusion found for this PISCO assignment |
| South Point SMR | SRI_SOUTH_POINT | SRI_SOUTH_POINT_E;SRI_SOUTH_POINT_W | No whole-group exclusion found for this PISCO assignment |
| South Point SMR | SRI_TRANCION_CANYON | SRI_TRANCION_CANYON_E;SRI_TRANCION_CANYON_W | No whole-group exclusion found for this PISCO assignment |

## Evidence and remaining questions

- Original master site table: exact `site_new` group links, retaining program, MPA and role. Leading/trailing whitespace is normalized only for map matching; original source fields remain unchanged.
- Raw fish, swath and size files: constituent record counts and observed years independently checked; all source hashes retained. See constituent_station_evidence.csv. These are raw years before seasonal/response filtering, not final analysis coverage.
- Source exclusion vectors: PISCO references, whole MPAs and KFM stations are recorded separately. No sign/significance rule is introduced.
- Scott Hamilton's October 29 email supports retaining the Tyler Bight and Yellowbanks habitat-based exclusions despite the existence of useful data.
- Emily's October 29 email describes intended before/after station continuity, with fish and benthic coverage allowed to differ. It does not explain every KFM station-specific omission.
- Jenn's October 19 email says KFM and PISCO survey observations are separate. Her reply does not establish that their physical footprints are identical.
- Targeted coauthor email searches found no explanation of Trailer or Willows and no contrary instruction to collapse these group labels. Search absence is not evidence that no private discussion occurred. Original email bodies remain private; durable IDs are recorded below without copying correspondence.

Only Trailer and Willows still need an identity/history answer: what did those labels refer to, and were observations ever collected under another name? If worth pursuing, ask the PISCO mapping/data custodian (Jenn or the current data manager); do not ask coauthors to reconfirm the 19 source-resolved group links. The missing count/length/frond measurements identified in other audits remain separate; this name audit did not recover them.

## Files and reproduction

- assignment_resolution.csv: all 21 assignments, identity classification, separate exclusion status and recommended treatment.
- constituent_station_evidence.csv: every exact group-to-station link by source, raw years, exclusion flags and existing inventory membership.
- remaining_name_searches.csv: bounded raw-name searches for Trailer and Willows.
- verification.json / source_hashes.json: assertions, source fingerprints and preservation checks.

Run `python3 code/tools/resolve_unmatched_site_assignments.py`. The next inventory revision should represent parents as group records linked to station records, preserving response-specific coverage. No analytical joins should substitute a parent for a constituent station or count both as samples.

Private-source message IDs checked: Scott Hamilton 1a188f9a2394d9e3; Emily Donham 1a189059ab077678; Jenn Caselle 1a15475000fd8c4a. No emails sent.
