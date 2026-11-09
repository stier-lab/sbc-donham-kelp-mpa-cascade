# Consolidated numerical verification — November 2, 2026

Passed:206 candidate keys,175 retained histories,31 exclusions; all3225 approved annual rows;175 station-support histories;1839 paired-date summaries;256 field comparisons covering46492 cells. Consolidated CSVs, supplement preview and embedded viewer values agree with their dated evidence, including missingness. Coverage year sets, gaps, before/after counts, reference unions/intersections, omission shifts and interval arithmetic independently recalculated.

All950 diagnostic coefficients, standard errors, interval endpoints, sample sizes and degrees of freedom independently refitted in R; maximum numerical difference8.08e-12. All447 reference-omission estimates and intervals independently checked; maximum difference3.12e-14. Original inventory checks pass for1856 site rows,32927 annual metadata rows,41 unchanged input hashes and487 protected manuscript/result files. Viewer tabs/search checks pass.

No numerical result discrepancies found. Standardized empty strings and missing cells in the viewer and increased embedded numeric precision to15 decimal digits. No biological data or manuscript results changed.

Limits: agreement with source records does not establish that historical protection metadata, unmatched habitats, seasonal choices or sparse/shared size samples are scientifically resolved. Those remain explicitly flagged; this is not a claim of independent validation of every original field observation. No new email/permission check was performed.

Reproduce with verify_consolidated_models.R, verify_consolidated_inventory.py, verify_dataset_inventory.py, verify_all_program_omissions.R and check_dataset_inventory_viewer.js. model_accuracy_checks.csv records every fit; consolidated_accuracy_checks.csv and consolidated_accuracy.json record the table checks.

Biomass extension:94 histories now linked;104125 size-donor occurrences checked against source records. Two isolated missing-stipe sensitivity fits independently verified in R. See ../biomass_support_20261102/README.md.
