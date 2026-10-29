# Matched predator-size refits

## Biological result

Restricting predators to large individuals did not consistently strengthen the temporal or predator-purple-urchin relationships. This is not evidence that body size is unimportant. The comparisons have limited spatial replication, reconstructed lobster measurements and uncertainty that depends on the fitted model.

The expanded mean-level comparisons contain ten MPAs for each density relationship and eight for biomass. The large-lobster density slope is -0.502 (nominal 95% CI -0.991 to -0.012; P = 0.0458), but its four-link adjusted P is 0.151. None of the four mean-level large-predator associations survives that correction. Omitting Naples changes the large-lobster density interval to -1.035 to 0.040; omitting Naples and Campus Point gives -1.074 to 0.060. These are influence checks, not a fitted covariance correction.

All calendar-AR1 large-minus-all temporal slope intervals include zero. An expanded sheephead-biomass difference is positive under MPA-cluster CR2 uncertainty (0.0263 log units per year, 95% CI 0.0062 to 0.0465; adjusted P = 0.0325), but not under AR1 (-0.0166 to 0.0646). Thus the size-difference conclusion is not stable across these uncertainty models.

The legacy CI-route alternative gives a negative lobster-biomass association for both all sizes (slope -0.403; adjusted P = 0.0067) and large lobsters (-0.343; adjusted P = 0.0381), and positive sheephead-density associations for both. Those results do not show a stronger large-predator effect. They mix mean levels and fitted changes across MPAs and condition on the original IID route-selection rule; they differ from the mean-only comparisons and are explicitly exploratory.

## What was refitted

- Temporal models use the Figure 4 species-by-time structure, with MPA intercepts/slopes where convergence and covariance checks permit. All-size and large-size versions use exactly the same history/year/time keys. Source is PISCO only, so no unidentifiable one-level programme random effect is added.
- Retained-membership scope preserves primary PISCO taxon-response eligibility but restricts annual predator records to years with usable size information. Its biomass model contains 542 rows across 44 histories and 15 MPAs; lobster biomass contributes only 25 years from three MPAs.
- Expanded scope adds the secondary reconstructed PISCO lobster-biomass records. Its biomass model has 634 rows across 51 histories and 15 MPAs; lobster contributes 117 paired years from ten MPAs. Density models have 558 rows across 45 histories and 15 MPAs in both scopes. These are not the full 831 observations of main Figure 4.
- Both temporal scopes compare all-size and large-size slopes, then model their paired within-year differences. Independent-error intervals use contrast-specific Satterthwaite degrees of freedom. Calendar AR1 checks preserve missing calendar years; CR2/HTZ checks allow arbitrary within-MPA covariance with small-sample degrees of freedom. Negative-Hessian and singular candidates cannot supply the retained independent-error results.
- Four predator-purple-urchin relationships are refitted for all/large sizes and both scopes. Each MPA comparison uses identical paired years. The mean-only family compares post-designation mean log ratios. The legacy alternative chooses a mean versus 11-year linear contrast from the all-size predictor and locks that route across thresholds; resource routing is shared as well. HAC uncertainty is conservative relative to IID variance, but remains conditional on selection. Extrapolated 11-year contrasts are flagged.
- REML meta-regressions use Knapp-Hartung inference and separate four-link BH families within scope, estimand and size group. Fewer than four MPAs yield an explicit not-estimable status. Retained-membership lobster-biomass pathways have only two PISCO MPAs and are not fitted; the four-MPA primary relationship also uses LTER, which lacks comparable size-filtered inputs.

## Measurement and inference limits

Sheephead are filtered at >30 cm total length; lobsters at >=82.55 mm carapace length. The lobster threshold is a legal-size reference, not a universal feeding threshold. All/large thresholds share the same positive-value offset within MPA/response. Surveyed zeros are retained; positive counts with absent size information are not converted to zero. No missing years are imputed.

The expanded lobster reconstruction remains different from primary biomass. No primary data, result tables, fitted models or main Figure 4 were replaced. Refits condition on reconstructed annual values, without propagating all size-pool, count, allometric or detection uncertainty. Sparse size pools can give misleadingly precise reconstructions. Pathway models treat estimated predator contrasts as fixed predictors; their exported variances are not incorporated as errors in variables. Shared reference reefs limit independence between MPAs. Temporal and pathway omission results address influence, not complete shared-reference covariance. None of these models measures feeding rates or identifies a designation-caused trophic cascade.

## Reproduction and outputs

Run modules 38 and 39, then `Rscript code/R/41_predator_size_refits.R`. All three are wired into `code/R/run_all.R`. Module 41 checks upstream source fingerprints, records unavailable inputs explicitly and never presents an unavailable run as completed. `renv.lock` now records clubSandwich 0.7.0; exact fitted package versions are exported.

- `tables/table_s_size_refit_temporal_annual.csv`: every matched annual model row.
- `tables/table_s_size_refit_temporal_slopes.csv`: independent-error, AR1 and CR2 slopes, including paired size differences.
- `tables/table_s_size_refit_temporal_candidates.csv` and `table_s_size_refit_ar1_candidates.csv`: convergence, singularity/Hessian checks and selected structures.
- `tables/table_s_size_refit_temporal_predictions.csv`: fitted-mean intervals explicitly labelled by model; not future-observation prediction intervals.
- `tables/table_s_size_refit_shared_reference_slopes.csv`: temporal CR2 results omitting Naples and Campus Point together.
- `tables/table_s_size_refit_pathways.csv`, `table_s_size_refit_pathway_effects.csv`, `table_s_size_refit_pathway_annual.csv` and `table_s_size_refit_pathway_omissions.csv`: estimates, MPA-level inputs, calendar membership and single/shared-reference omissions.
- `tables/table_s_size_refit_inputs.csv`, `table_s_size_refit_versions.csv` and `table_s_size_refit_status.csv`: fingerprints, versions and completion state.
- `plots/fig_s_size_refit_temporal.*` and `plots/fig_s_size_refit_pathways.*`: repository-only displays. The temporal display uses AR1 bands; the pathway display uses the mean-only family. Legends are in `docs/size_refit_figure_legends.md`.

Run `Rscript code/tests/test_predator_size_refits.R`. It checks exact matched cohorts, offsets, source hashes, model-labelled uncertainty, accepted covariance structures and guarded small-MPA fits. The full primary pipeline was not rerun for this secondary addition; module 41 was run directly and primary output hashes were checked unchanged.
