# California MPA kelp-forest analysis

Manuscript: *Kelp resilience and trophic responses in a southern California marine protected-area network* (target: Journal of Applied Ecology).

This repository analyzes predator, urchin and kelp histories in protected and reference areas. It contains analysis code, review diagnostics, figures and reproducibility records. Contributor-data processing lives in the companion `sbc-donham-kelp-mpa-data` repository; manuscript production is separate.

Start with the [current analysis guide](docs/ANALYSIS_GUIDE.md) and [searchable decision tables](outputs/dataset_inventory_20261101/index.html).

The source repairs and approved survey/model rules are combined in a reproducible [integrated review pipeline](outputs/integrated_review_20261103/index.html). Run `python3 code/tools/run_whole_dataset_biomass_review.py` to rebuild registered sources, comparisons, supported biomass, joint uncertainty, matched-year pathway checks and conditional heatwave refits. The [whole-dataset review](outputs/whole_dataset_biomass_update_20261105/index.html) shows the corrected results and change steps. The paper and frozen numerical outputs remain unchanged; review results are not released manuscript results.

- [Workflow status and reproducible entry points](outputs/dataset_inventory_20261101/workflow_status.csv)
- [Inclusion gates](outputs/dataset_inventory_20261101/inclusion_gates.csv) and [biomass calculations](outputs/dataset_inventory_20261101/biomass_calculations.csv)
- [Inclusion, measurement and borrowing rules](outputs/dataset_inventory_20261101/current_rules.csv)
- [Comparison decisions and exclusion reasons](outputs/dataset_inventory_20261101/comparison_decisions.csv)
- [Site coverage and survey characteristics](outputs/dataset_inventory_20261101/site_response_coverage.csv)
- [Remaining questions](outputs/dataset_inventory_20261101/unresolved_issues.csv)
- [Technical source-to-biomass workflow](../sbc-donham-kelp-mpa-data/docs/RAW_TO_BIOMASS.md)
- [Output folder status](outputs/README.md) and [script index](code/README.md)
- [Documentation catalog and archives](docs/README.md)

To regenerate the current inventory from existing audit evidence, run `python3 code/tools/report_dataset_inventory.py`, then `python3 code/tools/check_current_decisions.py`. These commands do not rerun ecological models. Python dependencies include pandas and numpy. Use the existing R environment for analysis scripts; historical model reproduction instructions remain with each dated report and the [preserved repository README](docs/archive/2026-11-02-consolidation/repository_README_before_cleanup.md).

Do not promote new data, cohorts or outputs simply because a source is now available. Follow the guide's review and release steps. Raw data access and redistribution remain subject to contributor conditions.
