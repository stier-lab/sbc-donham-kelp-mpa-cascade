# Fixed-station PISCO sheephead comparisons

The pooled increase in sheephead biomass persists when Campus Point, Naples and Point Vicente are reconstructed from consistent stations and refitted on matched years. Neither sheephead–purple-urchin relationship is resolved. This establishes that some PISCO fish comparisons can support before/after models; it does not establish that their reference reefs are equivalent or settle the designation-year boundary.

## Design and checks

Module 59 reconstructs sheephead counts and mass from the original PISCO fish file, using non-canopy levels, the provider species-attribute eligibility, the source transect–zone–site averaging, the 60 m² footprint and W(g) = 0.0144 × TL(cm)^3.04. Complete survey events supply observed zeros. Its full-mixture annual density and biomass match every retained raw value to <1e-7, and all six retained log-ratio histories match to <1e-7. Source files remain unchanged. Invalid count/length site-years are excluded from fixed-station eligibility; no unsized positive fish are silently treated as absence in the selected histories.

Station selection uses coverage alone. Enumerate every nonempty protected/reference station group, require at least four common post-designation years and three or four common baseline years, maximize the number of stations, then the number of common years, with alphabetical ties. Retain only years observed at every selected station. All feasible groups and the chosen membership are exported. This favors spatial replication; it is one explicit selection rule, not an optimization of biological results. It does not test every possible station subset's response. Shared reference stations and overlapping station pairs are not added as independent MPA replicates.

Inclusive baseline ends in 2012; strict baseline ends in 2011 and treats the observed 2012 survey as after. These are alternative boundaries, not evidence about actual survey timing relative to enforcement. Calendar time remains max(year − 2012, 0). Trends are evaluated between the first observed after year and the last observed year no more than 11 years later. No unobserved endpoint is substituted.

For each selected station group, compare the original mixture and fixed stations on exactly the same years. Within each resulting history, compare CI and eligible pBACIPS fits on exactly the same annual observations. Original-model-family fits are also exported so membership effects can be examined without model reselection. The adaptive zero correction is fixed at half the smallest positive value in the complete original MPA/response history; changing station or year membership cannot change that correction.

CI selects a mean log ratio unless its calendar-time slope is detectable. pBACIPS compares Step, Linear, Asymptotic and Sigmoid by AICc after the retained baseline-trend screen. Failed nonlinear fits are logged, never silently substituted. Mean CI estimates describe protected/reference differences; steps and trends describe changes. Their magnitudes are not estimates of the same biological quantity. A nonsignificant short-baseline trend does not establish parallel counterfactual trends. IID and calendar-gap HAC uncertainty are both exported; all intervals remain conditional on annual measurements, membership and model selection.

## Coverage

| Rule | Campus Point | Naples | Point Vicente |
|---|---|---|---|
| 3 baseline years including 2012 | 2 protected + 2 reference, 13 years | 3 + 2, 14 years | 2 + 3, 14 years |
| 4 baseline years including 2012 | 2 + 2, 13 years | 3 + 2, 14 years | 2 + 2, 15 years |
| 3 strictly prior years | 2 + 2, 13 years | 3 + 2, 14 years | 2 + 2, 15 years |
| 4 strictly prior years | Ineligible | Ineligible | 2 + 1, 16 years |

Campus Point uses IV_REEF_E/W and ARROYO_QUEMADO_E/W, with 2009–2013 and 2016–2023. Naples uses NAPLES_CEN/E/W and the same two reference stations, with 2009–2013 and 2015–2023. Point Vicente uses LONG_POINT_E and POINT_VICENTE_W. Its reference group is HAWTHORNE_REEF/RIDGES_N/ROCKY_POINT_N for inclusive minimum three (2009, 2011–2023), HAWTHORNE_REEF/RIDGES_N for inclusive minimum four and strict minimum three (2009–2023), and RIDGES_N for strict minimum four (2008–2023). These groups are constant within each alternative, but differ across alternatives; compare models within alternatives to isolate model choice.

## Findings

For the four-year inclusive rule, Campus Point biomass has a mean CI log ratio of −0.118 (SE 0.306), versus a before/after step of 0.367 (SE 0.683). Neither resolves a difference or a change. Naples has a mean CI log ratio of 1.919 (SE 0.176), versus a step of 0.747 (SE 0.345). Thus much of its protected/reference biomass difference predates designation. Point Vicente biomass selects a linear contrast of 1.131 (SE 0.393) under either model route. Point Vicente density fails the inclusive baseline-trend screen narrowly (p = 0.04974); its CI estimate is retained in the pooled sensitivity rather than forcing a before/after model. Under four strictly prior years, its selected density step is negative (−0.456, SE 0.313). Short baselines and boundary choices limit site-specific interpretations.

Replacing only these reconstructed histories preserves all nine pooled response directions. Across the four fixed-station pBACIPS alternatives, the sheephead biomass log contrast is 0.866–0.905 with IID adjusted p = 0.0016–0.0027; calendar-HAC checks likewise retain the positive result. Sheephead density remains unresolved. The sheephead–purple-urchin biomass and density intervals include zero under every alternative and uncertainty check. The lobster pathways are numerically unchanged because no lobster or urchin history is replaced.

The joint model shares variance across taxa. The inclusive minimum-three, inclusive minimum-four and strict minimum-three pBACIPS alternatives change adjusted support for lobster density and red-urchin biomass (IID adjusted p ≈ 0.039–0.042), despite no change to those histories. This is model-dependent precision, not new lobster or urchin observations. Four strictly prior years retain the original adjusted classifications. The full pooled/pathway exports include both uncertainty choices.

## Implication for the paper

Keep the present primary estimates while the baseline rule is resolved. Describe PISCO CI routing as an analytical choice, not a universal lack of before-data. The supplement should state that fixed stations preserve the sheephead biomass result, that strong protected/reference differences can predate protection, and that uncertainty in marginal cross-taxon results depends on model assignment. Do not choose the rule by significance. VERIFY WITH EMILY: historical three-versus-four-year criterion and whether designation-year surveys should count as baseline remain separate from this completed reconstruction.

Rebuild: `Rscript code/R/59_pisco_fixed_stations.R`, then `Rscript code/tests/test_pisco_fixed_stations.R`. Provider inputs default to the companion data repository; set `MPA_DATA_REPO` for another location. Outputs: `tables/table_s_pisco_fixed_{coverage,feasible_groups,station_values,annual,reconstruction,fits,candidates,decisions,pooled,pathways,fit_log,inputs}.csv`. Module 59 is included in both review and full-analysis runners. These are sensitivity products; no primary table or annual input is overwritten.
