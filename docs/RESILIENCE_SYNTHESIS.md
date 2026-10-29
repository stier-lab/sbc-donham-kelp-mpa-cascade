# Resilience analysis synthesis

Interpretation and numerical outputs updated 2026-10-29 using SBC127.5 monthly model-1 dry-foliar coefficients.
Five main figures are retained. Figure 5 uses diver-derived giant-kelp biomass;
the essential supplement uses S2.3 for design checks and S4.2/Table S6 for the historical-event comparison. Module 45 non-target comparators remain repository-only following Adrian's 21 October 2026 decision to remove them from V27 pending discussion with Bob. Longer legacy supplements have different numbering and must not be mixed with the current reading package.

## What the data support

Kelp had stronger baseline-relative resistance and recovery inside the sampled
MPAs than at their paired reference reefs during and after 2014-2016. This is an
association with fishing restrictions, not evidence that designation insulated
reefs from all climate or disease stressors. The same-pair Landsat comparison
with the 1997-1998 event did not detect a designation-associated increment in
relative disturbance response. Both primary historical intervals include zero
and are broad; neither improved causal resilience nor no benefit is established.

| Question | Module | Interpretation |
|---|---|---|
| Modern-event diver resistance/recovery | 22 | Ten kelp pairs; inside/reference median baseline ratios are 2.27/0.87 during the event, 2.34/0.80 in early recovery and 2.15/0.63 recently. Three focal signed-rank tests retain detection after correction. |
| Baseline and test dependence | 22 | Directions persist across baseline and deletion checks. The recent-period sign test does not detect a difference; do not call the result test-invariant. |
| Historical event-response increment | 28 | Resistance: 0.422 [-0.575, 1.419], 14 pairs; recovery: 0.505 [-0.387, 1.398], 15 pairs. Both q=0.377. |
| Same-history change versus post-only level | 27 | The 17-MPA arithmetic before/after kelp estimate is positive; 11 MPAs contribute satellite data only. The six-MPA diver-only and descriptive AR1 GLS intervals include zero. GLS re-estimates period differences rather than only their uncertainty. A before/after change does not establish parallel untreated trajectories. |
| Non-target reef comparability | 45 | Fixed legacy KFM reefs support three kelp-matched/four urchin-matched MPAs. Sparse non-target taxa and window-sensitive kelp pretrends leave comparability unresolved; no siting correction or validated negative control is supplied. |
| Geography and source coverage | 27 | Programme, designation timing and geography remain confounded. PISCO common-calendar post-only checks are narrower than a geographic BACIPS test. |
| Shared controls and influence | 26/27 | Known shared references reduce independent replication. Group uncertainty and deletion checks retain all observations. Sparse subsets remain descriptive. |
| Heatwave-window model and method bridge | 14/15 | Contextualize the modern-event contrast and analytical choices; neither establishes a causal designation effect. |
| Second heatwave | 19 | Kelp contrasts have consistent direction, but the second event is not detectably distinct from the intervening period. It is not an independent causal replication. |
| Disease context | 20 | No detected moderation by pre-disease sunflower-star abundance. Low power and spatially variable disease preclude ruling out disease contributions. |
| Temporal variability | 21 | No clear kelp stability difference; variability is distinct from resistance or recovery. |
| Reserve-level consistency | 23 | Repeated-response correlations are exploratory and do not establish ecological memory or priming. |
| Environmental moderators / prediction | 16/18 | No screened environmental moderator survives correction; small cohorts and method-dependent predictions limit explanation. |
| Cross-study context | 17 | Paired-design signs provide context; predator-size mechanism attribution is specification-sensitive and is not directly tested by this paper. |

## Historical sensitivities

The uncorrected Landsat index retains true annual zeros. Log contrasts require
positive window means. Cabrillo is retained in the primary comparison but cannot
enter a 2013-only baseline sensitivity because that baseline is zero. The
common-cohort table separates membership changes from baseline changes. Its
resistance contrast is already larger under the primary windows when Cabrillo
is omitted. Report that influence, rather than select a nominally significant
baseline or deletion result. Primary endpoints, alternative windows and regional
subsets have separate declared correction families. Sign flipping assumes
symmetry, which historical resistance does not clearly satisfy.

## Limits

MPA placement was not randomized. Existing monitoring stations are imperfect
reference matches, some controls are shared, fishing effort is unavailable and
warming/disease overlap designation and reserve maturation. Historical events
differ in intensity, duration and oceanography. The 2010-2013 baseline overlaps
mainland designation. Satellite canopy and diver biomass sample different
footprints and are not interchangeable physical measurements. Predation rates,
grazer consumption and predator-size mechanisms were not tested. Consistency
with parts of a trophic cascade is not direct evidence that the mechanism operated.

Kumagai et al. (2024) and Ortiz-Villa et al. (2025) supply California context, not
a claim that this paper recreates or independently confirms their complete
analyses. Our local pairs and historical comparison answer distinct questions about population patterns, cross-taxa links and designation-related disturbance responses.

## Rebuild and paper connection

Run `Rscript code/R/run_all.R` for the main pipeline. The shared resilience list
includes 28; scripts 27/29/30 add matched checks, the atlas and assertions.
`Rscript code/R/run_review_concerns.R` reruns the new checks alone from current
primary outputs and tracked inputs. `run_resilience.R` also offers optional
external-input modules 17/18/20. See `RESILIENCE_PIPELINE_INTEGRATION.md` and
`COAUTHOR_ANALYSIS_CHECKS.md` for the claim/output map. Essential supplement S1-S4 covers measurements, population contrasts, trophic relationships and disturbance responses. Supporting Atlas A supplies every retained pair history; it is not a numbered S15 figure in the current package.

Recovery rates/return times and formal hysteresis tests remain possible future
analyses, not results of the current paper or requirements silently added here.
