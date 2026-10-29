# Reference assignment audit

The assignment ledger covers all 145 primary histories and 39 MPA–program combinations. Rebuild it from the checksum-bound reference registry with:

```sh
Rscript code/R/export_reference_registry.R /path/to/data-processing-repo
Rscript code/R/27_design_geography_sensitivity.R
python3 code/audit_reference_assignments.py
```

Output: `tables/table_s_reference_assignment_audit.csv`. This is a documentation audit, not a new matching procedure. It separates identifiable sites, the documented rule used to select them, and ecological comparability. The last column of evidence cannot be inferred from the first two. The ledger's independent-validation flag describes this audit, not whether a monitoring program originally designed matched comparisons.

## Documented rules

The original `pBACIPS_PISCO_V10.R` records the PISCO suitability exclusions as advice from Scott Hamilton (lines 24–28 and 570 onward). Its mapping-table comment (line 49) records decadal-review changes and deliberate reuse of controls across MPAs. These repeated assignments must not be counted as independent references.

For KFM, the original script distinguishes shorter records where longer histories exist and sites on the wrong side of an island where a better reference exists (lines 18–22, 1124–1128 and 1459 onward). It does not assign one of those reasons to every excluded station. These are not all old Anacapa reserve stations. Current processing comments now state the supported rule classes.

For LTER, original line 2160 explicitly excludes Arroyo Hondo because kelp was essentially absent, attributing that decision to Bob Miller. The rule is stronger evidence for excluding that station than for ecological equivalence of the retained references. `a-l-02` and `AHND` identify the excluded station; both identifier systems must be checked. Other LTER aliases in the mapping likewise do not represent additional independent reefs.

Satellite footprints use a separate geometric selection procedure; see `docs/SOURCE_RECONCILIATION_2026-10-17.md` and the October20 source metadata. Matching at the MPA label does not mean the diver and satellite surveys sampled the same footprint. The documented Naples–Campus Point reference overlap includes 320 satellite pixels; the registry's pair labels alone do not encode that overlap.

## Correction and verification

The registry exporter had omitted the standalone LTER exclusion-code list, marking six AHND rows included even though the numerical processing excludes them. It now reads both lists. After regenerating the registry and script27, all numerical design-analysis columns and primary/harmonized numerical inputs were unchanged. The changed entries are exclusion flags, reference identifiers and file hashes. No reference was reassigned, no observations added or removed, and no primary model changed.

## Influence versus suitability

- `table_s_coauthor_gull_omission.csv`: omitting Gull Island KFM preserves all nine pooled directions and adjusted classifications, retaining the separate PISCO series. Keep the inherited primary/omission decision explicit; stability does not validate the contested Yellow Banks pairing.
- `table_s_pisco_control_source_removal.csv` and `table_s_pisco_control_leave_pair_out.csv`: assess dependence on a source and on individual PISCO comparisons, not ecological similarity. Removing a source changes coverage.
- `table_s_design_reference_groups.csv` and `table_s_design_leave_reference_group_out.csv`: represent known shared references. Campus Point and Naples share PISCO Arroyo Quemado East/West. Different site codes across programs must not be assumed to prove physical independence.
- `table_s_design_matched_estimates.csv`: before/after checks on matched histories address baseline differences without establishing parallel untreated trajectories.

Remaining ecological selection rationales are recorded per series. New alternatives should be chosen from defensible habitat, depth, exposure, location and pre-designation criteria, before evaluating their post-designation effects. Preserve the current cohort and distinguish any future reference change from analysis-routing or time-clock changes. Private correspondence and detailed collaborator questions remain in the manuscript repository.

Malakhoff and Miller (2021) included Gull Island and Yellow Banks in their four-reserve versus ten-control design. Their methods and archived author code establish membership, not the suitability of a one-to-one Gull–Yellow Banks pairing.
