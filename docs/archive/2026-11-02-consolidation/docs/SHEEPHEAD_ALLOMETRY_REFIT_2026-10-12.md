# Sheephead Allometry Correction and Refit

Implemented 12 October 2026 following approval to correct the PISCO sheephead conversion.

This is the preserved 144-history allometry phase. Adrian subsequently approved
adding Naples LTER sheephead biomass; the current 145-history analysis and its
separate omission sensitivity are documented in `NAPLES_BIOMASS_INCLUSION_2026-10-12.md`.
The phase-specific allometry audit should be replayed on its preserved checkpoint,
not presented as a current membership audit after that addition.

## Adopted Equation

Both PISCO and LTER use **W(g) = 0.0144 x TL(cm)^3.04**, from Hamilton & Caselle (2015), Methods page 2, DOI https://doi.org/10.1098/rspb.2014.1817. Individual mass is multiplied by recorded size-class abundance before aggregation. KFM supplies sheephead density only.

The provider's original species table was not rewritten. Its MD5 remains `6ccfacd18fff07208f510384e914a3b2`. The harmonization helper overrides only the in-memory SPUL row; other fish retain their supplied parameters. The unsupported 0.0238/2.8847 pair remains documented as historical input, not attributed to Hamilton or Williams. No MPA exclusion, survey filter or site pairing changed, including the separate Naples LTER exclusion.

## Rebuild and Checks

The staged data pipeline completed in 2.24 minutes (`pipeline_20261012_124447.log`). The scoped promotion retains the existing CSV order and all non-target record bytes. Its source, old/new parameters, units, count weighting and file hashes are recorded in `data/harmonized/sheephead_allometry.json`; individual changed records are in `SHEEPHEAD_ALLOMETRY_INPUT_CHANGES_2026-10-12.json`.

| Export | Records retained | Changed PISCO sheephead/Bio records | Non-target records byte-preserved |
| --- | ---: | ---: | ---: |
| Raw annual/status responses | 6,280 | 398 of 408 | 5,872 |
| Annual response ratios | 3,341 | 204 of 204 | 3,137 |

The ten unchanged raw target records have zero biomass. Other taxa, all densities, LTER sheephead, site metadata and satellite files are unchanged. The monthly kelp calibration marker is retained; its hashes are the parent of the explicitly recorded sheephead correction, not rewritten as if this were the original kelp-only release.

The initial model run stopped at the stale reference-registry guard. The registry was explicitly re-exported from the same mappings and survey footprints: all 857 reference records, 144 primary identities, exclusions and reference years are retained. The restarted **full pipeline completed with exit 0 in 7.03 minutes**, log `logs/pipeline_20261012_125555.log`. Current-input predictor and compound-disturbance exports were also refreshed.

An independent numerical audit retains all 144 primary histories, their sample counts, model choices and analysis designs. Fifteen individual histories change, all PISCO sheephead biomass. The other 129 individual contrasts and SEs are unchanged. All nine pooled-response direction and FDR classifications are unchanged; the BH adjustment was independently reconstructed. Density coefficients are unchanged, although their adjusted P values can shift because the nine-test correction also includes biomass tests.

## Numerical Impact

The exponentiated pooled sheephead biomass contrast changes from **3.3474 to 3.4366**; its new 95% interval is **2.1632 to 5.4595**, adjusted P = **0.0000100077**. These pooled contrasts combine different design-specific estimands; they are not a common physical inside/reference biomass ratio.

All five pooled biomass coefficients shift because they share the joint multilevel fit. The sheephead-biomass/purple-urchin-biomass meta-regression slope changes from 0.2884 to 0.2835 (P = 0.20093 to 0.20206; nine MPAs). Its interpretation remains inconclusive. The other five primary pathway regressions and all six memberships are unchanged.

Modern kelp resistance/recovery summaries and historical Landsat comparisons are unchanged. The correction does not establish a predator-urchin mechanism or resolve the separately documented Naples exclusion, reference-habitat or satellite-unit questions.

## Reproduction

Run raw harmonization into a staging directory using `DONHAM_HARMONIZED_OUTPUT`, then use `code/promote_sheephead_allometry.py` in the source-data workspace against a preserved pre-correction checkpoint. Do not write directly over checked release files or modify the provider's CSV. Refresh the reference registry explicitly before running the full analysis.

```sh
Rscript --vanilla code/R/export_reference_registry.R /path/to/data-workspace
Rscript --vanilla code/R/run_all.R
Rscript --vanilla code/tests/test_sheephead_release.R
python3 -B code/audit_sheephead_allometry_refit.py --before /path/to/preserved/before
python3 -B code/audit_kelp_calibration.py --strict
```

Module 35 validates the adopted source, parameters, units, scope and current fingerprints before fitting. The full resilience runner uses the same checks. Figure-only regeneration also refuses a snapshot from different inputs; its fingerprints are location-independent. The before/after audit requires the dated local checkpoint, while routine analyses use only tracked inputs and release metadata.

Existing sparse-subgroup, failed-alternate-fit and diagnostic warnings remain visible. Successful execution does not establish that every model passed its assumptions. No push or circulation was performed for this correction.
