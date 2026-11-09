# Pelican source choice and PISCO missing-record repairs

November 2, 2026. Pipeline and review tables updated; paper and frozen results unchanged.

## Pelican: original 140 accepted

The source-selection policy now includes all 140 original Pelican Bay red-urchin measurements from August 7, 1995, as Adrian requested. The combined-file version with 135 animals remains available for sensitivity; no original file is overwritten. The source build contains 140 more selected rows than the preceding held build, but only five more animals than the older combined-file analysis.

The isolated, previously verified before–after sensitivity is retained: Scorpion KFM red-urchin estimate 1.253 with the 135-animal input versus 1.270 with the original 140. Both conditional intervals remain positive. This is not a prediction of the result after every pending pipeline change or a complete uncertainty refit.

## Cavern Point West: missing count recovered

The cached PISCO 1.11 swath release passes its provider SHA256 checksum and records **219 red urchins** for UCSB, Cavern Point West, September 16, 2013, INNER zone, transect 2. Emily’s two delivered swath copies have a missing count on that same keyed record.

The builder applies this one keyed recovery to a derived swath table. It preserves original source row IDs and records the supporting newer-release row ID. The provider file and metadata are registered immutable snapshots. The pipeline consumes the reconciled table; it does not substitute the entire newer release or add years/sites. A changed key, missing recipient, duplicate candidate or unexpected supporting count stops the build.

| Same station/year comparison | Previously zero-filled | Corrected |
|---|---:|---:|
| Cavern Point West red urchins, 2013 | 1.319/m² | 1.928/m² |
| Mean of six surveyed protected Scorpion stations, 2013 | 3.468/m² | 3.569/m² |

These use the inherited transect-to-zone-to-station aggregation on the same delivered observations. They are diagnostic density comparisons, not new pooled temporal-model fits. Biomass must be regenerated with the corrected count and its documented donor pool before reporting its final change.

## Kelp records: unknown stipes remain unknown

Little Dume West, September 15, 2011, OUTER transect 1 records two plants without stipe sizes. PISCO 1.11 retains that missing size. Underwater Arch, July 13, 2023, INNER transect 2 records two plants without sizes; the checked UCSB/UCSC provider release does not contain that VRG record. No replacement size was recovered. Plant counts alone cannot determine stipes or kelp mass.

The prior production branch already marked the affected role-year kelp biomass unavailable, but its stipe-density output could still retain the zero-filled calculation. The revised branch marks stipe totals/density unavailable too, without averaging around missing stations. Whole role-year handling preserves the existing conservative station-weight policy.

| Site-year stipe-density diagnostic | Previous zero-filled value | Revised status |
|---|---:|---|
| Little Dume West, 2011 | 2.622 stipes/m² | Unknown |
| Underwater Arch, 2023 | 0 stipes/m² | Unknown, not absent kelp |

The original plant counts remain intact in the immutable source. These are missing measurements, not ecological site exclusions or grounds to remove unrelated years/species.

## General safeguard

Survey-grid rows now retain a source-presence marker. Only a roster-generated species absence gets zero. An explicit source row with unknown count or unknown stipes keeps its missing value. Before role-year export, affected density and biomass are protected from legacy intermediate zero filling and available-station averaging. This addresses the audited cases; it does not claim that every missing-value convention across all programs has been resolved.

## Current pipeline status

The required-source gate passes: zero required holds. Four unused 2023 KFM groups remain deferred and cannot supply supported biomass for an extension. Researcher exclusions, VRG acceptance and current analytical scope remain intact. Next run a complete review build, verify donor/aggregation diagnostics and compare temporal-model results before updating the paper. No full harmonization/model rerun is claimed here.

## Evidence and reproduction

- [Recovered source count and both row identities](pisco_count_recovery.csv)
- [Site-year changes](site_year_changes.csv)
- [Scorpion station-by-station comparison](scorpion_2013_station_density.csv)
- [Retained Pelican sensitivity](scorpion_primary_step_sensitivity.csv)
- [Verification](verification.json)

Run `make check` in the data repository and `python3 code/tools/review_pelican_pisco_recovery.py` in the analysis repository. The audit tests original Pelican row count, checksum-verified recovery, unchanged exclusions, source-grid absence versus explicit unknowns, and role-year missingness. Frozen manuscript outputs remain unchanged.
