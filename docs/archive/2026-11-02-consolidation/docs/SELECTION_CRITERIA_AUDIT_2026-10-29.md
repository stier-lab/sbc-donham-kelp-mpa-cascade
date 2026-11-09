# Selection criteria and observed coverage

The primary cohort and estimates are unchanged. This audit separates reference
selection, observed response coverage, and model assignment. Those are different
decisions and cannot be represented by one before-year threshold.

## Retained histories

All 145 primary histories span at least ten calendar years and contain at least
four paired years after the recorded designation year. Twenty-four diver
histories at Anacapa and South Point lack observations within one year of that
date. They are CI comparisons. A universal within-one-year inclusion requirement
therefore misdescribes the analyzed cohort.

Anacapa's retained KFM record starts in 2005; the analysis uses the 2003 designation
year. These dates must not be conflated. The [CDFW reserve overview, page 1](https://nrm.dfg.ca.gov/FileHandler.ashx?DocumentID=100887)
confirms establishment in 2003 and re-establishment in the statewide network in
2012. This audit does not change either the data or the analysis clock.

## PISCO station continuity

`audit_selection_coverage.R` reconstructs exact station-year survey footprints from
the supplied fish/swath tables, current source filters, explicit excluded site
codes, and MPA assignments. It does not fuzzy-match or combine differently named
stations. Fish and benthic surveys are considered separately. It intersects
survey years with years observed for every retained response in the respective
branch, then counts years shared by the same protected and reference stations.

These are coverage upper bounds: a surveyed station-year does not prove that every
species has a usable biomass estimate there. Source processing includes a
five-year station filter across the combined fish/benthic data; that filter is
distinct from any minimum number of baseline years for a paired response.

For all retained benthic responses together, no station pair reaches three matched
baseline years, whether the designation year is included or excluded. Fish
coverage differs:

| Fish comparison | Maximum matched years strictly before designation | Including designation year | Pairs with at least four baseline and four after years, including designation |
|---|---:|---:|---:|
| Abalone Cove | 2 | 3 | 0 |
| Campus Point | 3 | 4 | 4 |
| Naples | 3 | 4 | 6 |
| Point Vicente | 5 | 6 | 7 |

These pairs share stations and must not be counted as independent replications.
An additional fish-only check removes station-years with invalid counts or missing
or nonpositive lengths for positive sheephead counts. It preserves every threshold
count in the table. Eight source rows flag six station-years; none of the matched
baseline years underlying the four-year fish comparisons is affected. The flagged
raw site-years are recorded separately, without asserting that all contribute to
primary estimates or that a survey roster proves complete ecological sampling.
Point Vicente has four pairs with at least four strictly prior years and four
after years. Thus neither a three- nor four-year coverage threshold explains all
PISCO CI assignments in the current files. The primary analysis uses CI for all
PISCO histories as a program-wide rule. The separate fixed-cohort routing checks
in `52_routing_clock_sensitivity.R` already test response-level alternatives;
this audit adds station continuity evidence, not new refits.

The implemented non-PISCO candidate screen uses at least three paired baseline
observations and tests baseline trends. Baseline includes the designation year.
Historical recollections of a four-year screening rule are not an implemented
four-year threshold. They do not establish ecological equivalence of stations or
justify a delayed protection date.

The source mapping also contains a blank/reference conflict for POINT_LOMA_CEN
under South La Jolla. It is outside the retained PISCO cohort; no retained mapping
has conflicting statuses. The conflict is recorded, not silently resolved.

## Judith Rock satellite exclusion

The original analysis explicitly attributes Judith Rock's exclusion to overlap
with the San Miguel Island SC extraction. The delivered coordinate workbook
confirms 364 shared MPA pixels (965 Judith Rock pixels and 1,082 San Miguel pixels).
Their Ref_1 sets share 1,889 pixels; all 726 Judith Ref_2 pixels occur in the San
Miguel Ref_2 set. These are extraction-mask comparisons, not verification of legal
MPA boundaries.

Both records are currently excluded. Consequently, overlap explains the original
comment but does not by itself justify continuing to exclude Judith Rock after
San Miguel is removed. Confirmation of the intended Judith/San Miguel footprints
and the final inclusion decision remains necessary. The audit does not relabel
San Miguel Island SC as San Miguel Island SMR or add either to the primary cohort.

## Naples reference basis

Program investigators confirm that the PISCO Arroyo Quemado references were
selected using experienced divers' assessments and community-composition nMDS
comparisons during the decadal evaluation. They replaced Campus Point in that
program's pairing. No close match to Naples was identified. The supplied mapping
and the technical report agree on Arroyo Quemado East/West; the current analysis
already uses them. The ordination itself has not been recovered or independently
verified here. This evidence does not establish why the separate LTER references
were selected. The manuscript states the qualitative basis and limitation without
recounting the history of the reassignment.

## What remains unresolved

- Individual LTER reference choices: assignments and consultation are documented,
  but per-site habitat/selection reasons have not been recovered.
- Judith Rock's current eligibility and San Miguel's intended reserve footprint.
- The historical three-versus-four-year screening recollection. The current
  three-observation code rule and program-wide PISCO exception are explicit.
- Historical boundaries or individual reference choices lacking direct evidence.

KFM fixed-station checks, Gull-reference omission and matched routing/clock
sensitivities remain in the pipeline. Their stable pooled directions do not prove
that reference habitats are equivalent or eliminate sensitivity of weaker claims.

## Reproduce

```sh
python3 code/tools/audit_selection_criteria.py
Rscript code/tools/audit_selection_coverage.R /path/to/sbc-donham-kelp-mpa-data
```

Python requires openpyxl. The supplied satellite workbook and source PISCO tables
must be available locally. The scripts write `table_s_selection_*` audit tables,
record input hashes, and assert that source inputs remain unchanged. Private email
evidence and manuscript editing records remain in the private manuscript repo.
