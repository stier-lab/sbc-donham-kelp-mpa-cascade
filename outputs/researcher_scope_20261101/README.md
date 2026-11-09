# Researcher exclusions retained; baseline-trend exclusions removed

November 1, 2026. Review only: no manuscript/SI, canonical annual input, or canonical result changes. This implements Adrian’s latest direction and supersedes the earlier provisional decision to retain Gull KFM. Compare against outputs/upstream_sites_20261101 (the expanded review), not the old mixed-estimand manuscript table.

## Selection decisions

- Remove Anacapa Island SMCA and Painted Cave SMCA: preserve the documented ecological-scope exclusion for permitted lobster take. Keep the distinct Anacapa Island SMR.
- Remove Point Conception SMR: preserve the original ecological-scope decision concerning rare focal predators, even though its satellite measurements are valid.
- Remove Gull Island KFM: preserve researcher concerns about its reference; retain the separate PISCO comparison.
- Retain Scott Hamilton’s established excluded PISCO reference stations; do not reintroduce Crook/Tyler/Yellowbanks or substitute new controls. Preserve Jenn’s revised Naples pairing.
- Keep histories formerly removed for baseline trends/significance and related-response trend screens. Keep Point Dume sheephead and other reconstructed responses with adequate coverage and no recovered ecological/advice exclusion.
- Arrow Point, Cat Harbor, Dana Point, and San Miguel access closure remain excluded. Unreconstructed/unresolved site assignments remain pending.
- Judith Rock stays in this main review: its documented overlap is with excluded San Miguel, and no independent ecological/advice exclusion was recovered. A complete companion rerun omits Judith. It does not change the interpretation.
- The preserved site-specific decisions include historical retained exceptions such as Naples and Swamis. This is not a newly applied universal all-species-no-take rule; do not describe the cohort that way.

## Coverage and comparison

Candidate response histories:206 to175. Primary before-after fits:62 to54. Post-protection fits:191 to162. These are response histories, not independent MPAs. The removed candidate roster includes two Painted Cave biomass histories that already failed primary coverage, hence31 removed candidates versus29 removed post fits.

| Area/program removed | Before-after fits | Post-protection fits |
|---|---:|---:|
| Anacapa Island SMCA / KFM | 0 | 7 |
| Anacapa Island SMCA / PISCO | 0 | 8 |
| Gull Island SMR / KFM | 7 | 7 |
| Painted Cave SMCA / PISCO | 0 | 6 |
| Point Conception SMR / Landsat | 1 | 1 |

No retained annual measurement or individual estimate changed. Differences in pooled estimates come from membership and the consequent weights/uncertainty, not recalibration or a different time window. The same conservative dependence clusters are retained; removing an intermediary does not create new independent replication.

## Results

Estimates and95%robust intervals are in log-ratio units. Adjusted p uses the same9response family.

| Comparison | Response | Expanded review | Researcher-scope review | 95% interval | Adjusted p |
|---|---|---:|---:|---|---|
| Change | Red urchin biomass | 1.175 | 1.140 | Insufficient independent information | Unavailable |
| Change | Red urchin density | 0.054 | 0.238 | -0.754 to 1.231 | 1.000 |
| Change | Kelp biomass | 0.328 | 0.341 | -0.147 to 0.829 | 0.865 |
| Change | Lobster density | 0.978 | 1.472 | Insufficient independent information | Unavailable |
| Change | Sheephead biomass | 0.237 | 0.237 | Insufficient independent information | Unavailable |
| Change | Sheephead density | 0.277 | 0.333 | -0.817 to 1.482 | 1.000 |
| Change | Purple urchin biomass | -0.597 | -0.519 | Insufficient independent information | Unavailable |
| Change | Purple urchin density | -1.137 | -0.916 | -2.698 to 0.866 | 0.865 |
| Post_difference | Red urchin biomass | 0.231 | 0.201 | -0.579 to 0.981 | 0.730 |
| Post_difference | Red urchin density | 0.062 | 0.008 | -0.742 to 0.758 | 0.980 |
| Post_difference | Kelp biomass | 0.152 | 0.258 | -0.314 to 0.829 | 0.548 |
| Post_difference | Lobster biomass | 0.199 | 0.199 | -1.315 to 1.713 | 0.766 |
| Post_difference | Lobster density | 0.658 | 0.737 | -0.233 to 1.707 | 0.354 |
| Post_difference | Sheephead biomass | 1.209 | 1.216 | 0.540 to 1.892 | 0.017 |
| Post_difference | Sheephead density | 0.655 | 0.710 | 0.253 to 1.168 | 0.023 |
| Post_difference | Purple urchin biomass | -0.436 | -0.395 | -1.066 to 0.277 | 0.476 |
| Post_difference | Purple urchin density | -0.274 | -0.262 | -0.895 to 0.372 | 0.548 |

Sheephead biomass and density remain higher inside protected areas after designation. Kelp, lobster and urchin post-protection differences remain uncertain. The kelp before-after estimate changes little. Removing Gull KFM leaves insufficient independent information for robust intervals for red/purple urchin biomass change; lobster density change now has only two contributing MPAs. Withheld intervals represent limited information, not no biological effect.

The flat-baseline kelp change is0.341 (−0.147,0.829). The same-row trend sensitivity is now0.156 (−0.528,0.840), compared with−0.382 in the expanded review; removing Gull KFM removes a strongly extrapolated baseline contribution. Both leave regional kelp change unresolved. Full follow-up produces0.659 (0.160,1.158), adjusted p=.122; report the duration sensitivity without substituting it for the common first-nine-year window.

Omitting Judith Rock gives kelp change0.422 (−0.068,0.913) and post difference0.317 (−0.263,0.897). These remain uncertain; the other taxa are unchanged because Judith contributes satellite kelp only.

Biological take-home: protected sites support more sheephead on average, but the regional before-after comparisons do not establish a consistent predator–urchin–kelp cascade or a clear kelp increase attributable to protection. The revised exclusion policy does not reverse that conclusion, but reduces precision for several responses. These are observational comparisons; before-after data alone do not remove all siting confounding.

## Figure

pooled_changes.png/pdf: gray upper points are the preceding expanded cohort; blue lower points apply the researcher/ecological exclusions. Lines show95%CR2 intervals; hollow points lack enough independent information for reliable inference. No embedded legend, title or source prose. Separate panels distinguish before-after change and post-protection association.

## Reproduce and verify

```sh
python3 code/tools/build_researcher_scope_review.py
Rscript code/tools/refit_common_estimand.R outputs/researcher_scope_20261101
Rscript code/tools/refit_common_estimand.R outputs/researcher_scope_20261101/omit_judith
python3 code/tools/check_common_estimand.py outputs/researcher_scope_20261101
python3 code/tools/check_common_estimand.py outputs/researcher_scope_20261101/omit_judith
Rscript code/tools/compare_diver_exclusions.R outputs/researcher_scope_20261101 outputs/upstream_sites_20261101
python3 code/tools/summarize_researcher_scope.py
```

Main review:1521 individual fits and252 pooled rows independently checked; Judith omission:1507 individual fits and252 pooled rows checked. All487 protected artifacts unchanged. All1521 retained individual fits reproduce the prior review; maximum estimate/SE difference below1e-9. Model fitting logged no warnings. Comparison figure visually inspected.

This reruns the common-estimand synthesis and its existing sensitivities. Resilience, pathway and other downstream manuscript analyses have not been promoted or rebuilt. The requested paper-update hold remains in place.
