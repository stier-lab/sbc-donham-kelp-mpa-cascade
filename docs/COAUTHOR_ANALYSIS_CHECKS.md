# Design, geography and historical disturbance checks

Added during coauthor revision, 2026-10-09. These analyses complement the primary
progressive-change meta-analysis; they do not silently replace its estimands.

## Rebuild

From this repository root:

```sh
Rscript code/R/run_review_concerns.R
Rscript code/tests/test_review_concerns.R
Rscript code/tests/test_review_concern_integration.R
```

The full `code/R/run_all.R` pipeline also runs scripts 27-30. Script 28 is listed
in both resilience orchestrators. All routine inputs are tracked: current primary
keys, harmonized annual responses, site metadata and the uncorrected Landsat index.
`data/harmonized/LANDSAT_ANNUAL_INDEX.md` records the one-time upstream reconstruction.
It reproduces the existing offset-corrected satellite histories exactly, while
preserving annual zeros for baseline-relative event calculations.

The integration test uses a temporary tracked-input fixture, compares standalone
and shared-environment CSVs, rejects masked package lookups, guards against
external raw-workspace reads, preserves inherited pipeline objects and exercises
deliberate corruption failures. It does not overwrite repository outputs.

## Questions and outputs

| Question | Module | Products |
|---|---|---|
| Does before-period information change the result on the same pairs? | 27 | `table_s_design_matched_*`; `fig_s_design_matched.*` |
| Do source and island/mainland differences persist with comparable coverage? | 27 | `table_s_design_geography.csv`, source support, post-only and common-calendar checks |
| Do shared references or individual MPAs dominate the new bridge? | 27 | Shared-reference mapping, group uncertainty and MPA/group influence outputs |
| Did relative disturbance response improve after designation? | 28 | `table_s_historical_heatwave_*`; `fig_s_historical_heatwave_bacips.*` |
| What are the observed histories of every analysed pair? | 29 | 145 histories on 37 four-facet pages, combined PDF, `table_s_pair_history_atlas_index.csv` |
| Are calculations and cohorts coherent? | 30 / tests | Fail-fast assertions and synthetic estimand/calendar/zero tests |

## Interpretation

The matched bridge compares a post-only level with an after-minus-before change.
Different estimates can arise from pre-existing inside/reference differences;
these are not two interchangeable causal estimates. Equal-source averages within
MPA prevent programme replication from being counted as extra independent MPAs.
Shared-reference sensitivities approximate uncertainty for known control overlap;
unequal group sizes, incomplete individual-site footprints and unverified satellite
polygon overlap limit calibration. They do not account for all spatial covariance.
Several endpoints are sparse. PISCO has no matched island
before/after histories; a common-calendar post-only comparison is narrower than
a causal geographic test. Secondary GLS checks disclose optimizer and covariance
instability rather than changing primary inclusion.

The historical event comparison uses the same satellite pairs and uncorrected
window means, with at least two paired years per window. Primary resistance has
14 pairs and a modern-minus-historical log advantage of 0.422 (95% CI -0.575 to
1.419); recovery has 15 pairs and 0.505 (-0.387 to 1.398). Both corrected p-values
are 0.377. These broad intervals do not establish a designation-associated
improvement, and also do not establish no benefit. The 2013-only baseline excludes
Cabrillo because its mean is zero. Common-cohort tables separate membership from
baseline changes; influence checks retain, rather than selectively remove, Cabrillo.

The three focal modern-event diver comparisons, two primary historical endpoints,
eight historical-window alternatives and geographic checks are distinct families.
Do not select an uncorrected significant sensitivity as the headline result.
Sign-flip tests assume symmetry, which is doubtful for historical resistance.

MPA placement is nonrandom, fishing effort is unavailable, and the historical
events differ in duration, intensity, oceanography and disease. The modern
baseline overlaps mainland designation. Satellite canopy and diver biomass
have different footprints. Predation, grazing rates and predator-size mechanisms
were not tested. Fishing restrictions are not protection from every stressor.

## Manuscript connection

The private manuscript source reports the new Methods, Results and Discussion
qualifications; SI S1.7 and S2.14 show methods and results. Figure S15 is the
complete history atlas; S16 is the matched-design bridge; S17 is the historical
event comparison. Five main figures are retained, with physical units now shown
for diver-derived kelp in Figure 5. Figure legends are also emitted by each module.
