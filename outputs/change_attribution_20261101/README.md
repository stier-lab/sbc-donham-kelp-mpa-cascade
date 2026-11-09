# What changes the urchin and lobster results?

November 1, 2026. Review only; the manuscript, supplement and canonical outputs remain unchanged. No eligibility tolerance selected.

## Main finding

The changes have different causes for different responses. Purple-urchin density and lobster density lose support when the legacy mixture of means, steps and fitted trends is replaced with an average post-protection difference and a common nine-year window. They do not lose support merely because cluster-robust or annual-correlation uncertainty is added. Red-urchin biomass is especially sensitive to the pooling structure and the added histories. Purple-urchin biomass also loses support under a consistent post-protection target before new sites are added. Lobster biomass and red-urchin density remain uncertain throughout the comparisons below.

These are controlled, sequential and factorial comparisons, not a universal percentage decomposition. Changes in pooling weights and estimated heterogeneity interact with membership and covariance. BH-adjusted p values depend on the full nine-response family, so q can change even when a response's estimate and raw p are unchanged. Use the point estimates and raw intervals alongside q.

All tables below report log-scale estimates, 95% intervals and BH q; zero denotes no difference/change. Legacy contrasts and post-protection differences answer different questions despite sharing a scale.

## 1. Pooling and uncertainty, exactly the same 143 legacy histories

The original joint model shares heterogeneity across taxa within biomass/density models and has MPA and program random effects. Separate-response pooling estimates heterogeneity for each response, with MPA and repeated series components. The next step changes only inference to CR2 with shared-reference clusters. The last changes only individual covariance to calendar HAC3; point estimates remain fixed individually, but their weights can change the pooled estimate.

| Response | Current manuscript | Separate-response pooling | Cluster-robust intervals | Annual HAC covariance |
|---|---|---|---|---|
| Purple urchin density | -1.560 [-2.498, -0.621]; q=0.005 | -1.325 [-1.830, -0.819]; q=0.000 | -1.325 [-1.780, -0.870]; q=0.001 | -1.379 [-1.846, -0.913]; q=0.001 |
| Purple urchin biomass | -0.584 [-1.126, -0.043]; q=0.052 | -0.826 [-1.328, -0.323]; q=0.009 | -0.826 [-1.289, -0.362]; q=0.011 | -0.828 [-1.266, -0.389]; q=0.010 |
| Red urchin density | -0.465 [-1.392, 0.462]; q=0.360 | -0.601 [-1.299, 0.096]; q=0.111 | -0.601 [-1.353, 0.150]; q=0.130 | -0.623 [-1.390, 0.145]; q=0.125 |
| Red urchin biomass | 0.571 [0.062, 1.081]; q=0.051 | 0.087 [-0.756, 0.930]; q=0.827 | 0.087 [-0.879, 1.054]; q=0.834 | 0.066 [-0.898, 1.030]; q=0.874 |
| Lobster density | 1.048 [0.119, 1.978]; q=0.051 | 1.124 [0.403, 1.845]; q=0.009 | 1.124 [0.318, 1.930]; q=0.019 | 1.125 [0.340, 1.911]; q=0.020 |
| Lobster biomass | 0.416 [-0.209, 1.041]; q=0.242 | 0.426 [-0.696, 1.547]; q=0.395 | 0.426 [-1.241, 2.092]; q=0.511 | 0.472 [-0.926, 1.870]; q=0.386 |

Red-urchin biomass changes substantially at the pooling-structure step, before any histories or annual observations are changed. Purple-urchin density remains supported after both uncertainty changes. Lobster density and purple-urchin biomass actually receive stronger adjusted support under these alternative pooling specifications. Robust uncertainty is not uniformly wider or more conservative.

`A0b_current_pool_HAC_floor.csv` separately reproduces the manuscript's annual-correlation convention: original joint pooling and max(original SE, HAC SE). This differs from the proposed raw-HAC convention; do not conflate them.

## 2. Effect definition and follow-up window, fixed intersecting roster

First retain only legacy histories estimable under the proposed target, then compare models without any further membership changes. For post differences this is 137 histories. The first two columns use identical annual observations for every history. Common models use a flat Before level and a post-period slope (or a post-period line when no Before observations exist). The target is the average fitted After ratio. Designation-year removal/calendar classification, nine-year restriction, and removal of unused Before observations then occur in distinct steps. All steps use the same separate-response pooling, HAC3 covariance and CR2 inference.

| Response | Legacy contrast | Mean After ratio | Strict calendar periods | After years 1–9 | After observations only |
|---|---|---|---|---|---|
| Purple urchin density | -1.379 [-1.846, -0.913]; q=0.001 | -0.902 [-1.614, -0.189]; q=0.058 | -0.902 [-1.615, -0.189]; q=0.058 | -0.433 [-1.093, 0.227]; q=0.255 | -0.430 [-1.089, 0.229]; q=0.258 |
| Purple urchin biomass | -0.828 [-1.266, -0.389]; q=0.010 | -0.698 [-1.491, 0.094]; q=0.112 | -0.699 [-1.492, 0.094]; q=0.112 | -0.345 [-1.081, 0.392]; q=0.383 | -0.343 [-1.078, 0.391]; q=0.384 |
| Red urchin density | -0.637 [-1.438, 0.164]; q=0.132 | -0.235 [-1.300, 0.830]; q=0.621 | -0.234 [-1.298, 0.830]; q=0.621 | 0.106 [-0.765, 0.978]; q=0.783 | 0.103 [-0.767, 0.972]; q=0.790 |
| Red urchin biomass | -0.022 [-0.988, 0.944]; q=0.957 | 0.272 [-0.516, 1.060]; q=0.554 | 0.272 [-0.516, 1.061]; q=0.559 | 0.621 [0.137, 1.105]; q=0.060 | 0.607 [0.119, 1.095]; q=0.068 |
| Lobster density | 1.161 [0.332, 1.990]; q=0.022 | 1.010 [0.133, 1.887]; q=0.066 | 1.011 [0.134, 1.887]; q=0.065 | 0.737 [-0.233, 1.707]; q=0.255 | 0.737 [-0.233, 1.707]; q=0.255 |
| Lobster biomass | 0.472 [-0.926, 1.870]; q=0.386 | 0.279 [-0.967, 1.525]; q=0.554 | 0.276 [-0.995, 1.546]; q=0.567 | 0.197 [-1.330, 1.725]; q=0.770 | 0.199 [-1.315, 1.713]; q=0.766 |

For purple-urchin density, changing the target on the same observations shifts the log contrast from −1.379 to −0.902; restricting follow-up shifts it further to −0.433. Adding annual-correlation uncertainty alone did not erase the old finding. For lobster density the corresponding sequence is 1.161 → 1.010 → 0.737. Both intervals include zero after standardizing the follow-up window.

## 3. Uncertainty alone on the same new-model observations

All cells below use the same 137-history roster, first-nine-year observations and effect definitions. IID versus HAC changes only individual covariance; model versus CR2 changes only pooled inference. CR2 leaves the fitted pooled point estimate exactly unchanged for any fixed covariance input.

| Response | IID / model interval | IID / CR2 | HAC / model interval | HAC / CR2 |
|---|---|---|---|---|
| Purple urchin density | -0.426 [-1.032, 0.180]; q=0.228 | -0.426 [-1.066, 0.215]; q=0.246 | -0.430 [-1.056, 0.196]; q=0.240 | -0.430 [-1.089, 0.229]; q=0.258 |
| Purple urchin biomass | -0.373 [-1.035, 0.289]; q=0.302 | -0.373 [-1.077, 0.332]; q=0.314 | -0.343 [-1.034, 0.348]; q=0.373 | -0.343 [-1.078, 0.391]; q=0.384 |
| Red urchin density | 0.145 [-0.605, 0.895]; q=0.685 | 0.145 [-0.672, 0.962]; q=0.689 | 0.103 [-0.692, 0.897]; q=0.786 | 0.103 [-0.767, 0.972]; q=0.790 |
| Red urchin biomass | 0.633 [0.183, 1.084]; q=0.029 | 0.633 [0.155, 1.112]; q=0.053 | 0.607 [0.147, 1.067]; q=0.042 | 0.607 [0.119, 1.095]; q=0.068 |
| Lobster density | 0.751 [-0.138, 1.639]; q=0.207 | 0.751 [-0.217, 1.718]; q=0.230 | 0.737 [-0.155, 1.628]; q=0.222 | 0.737 [-0.233, 1.707]; q=0.255 |
| Lobster biomass | 0.237 [-0.783, 1.257]; q=0.623 | 0.237 [-1.377, 1.852]; q=0.689 | 0.199 [-0.819, 1.217]; q=0.693 | 0.199 [-1.315, 1.713]; q=0.766 |

Purple-urchin density and lobster density are already uncertain using conventional IID/model intervals under the standardized target and window. Thus loss of support cannot be attributed simply to the robust uncertainty method. Red-urchin biomass is borderline across uncertainty choices on this intersecting roster.

## 4. Membership alone: removals and additions separated

Every retained history has exactly the same observations, individual effect, covariance, pooling specification and inference across these four columns. The old roster contains 141 post histories after the fixed nine-year coverage requirement removes the two Farnsworth sheephead histories. Researcher removals then remove four Gull KFM histories. Additions supply 25 histories. The final roster has 162 post histories. No baseline-trend tolerance screen is applied.

| Response | Old eligible roster | Removals only | Additions only | Final roster |
|---|---|---|---|---|
| Purple urchin density | -0.430 [-1.089, 0.229]; q=0.258 | -0.430 [-1.089, 0.229]; q=0.258 | -0.262 [-0.895, 0.372]; q=0.548 | -0.262 [-0.895, 0.372]; q=0.548 |
| Purple urchin biomass | -0.343 [-1.078, 0.391]; q=0.384 | -0.343 [-1.078, 0.391]; q=0.384 | -0.395 [-1.066, 0.277]; q=0.476 | -0.395 [-1.066, 0.277]; q=0.476 |
| Red urchin density | 0.127 [-0.724, 0.978]; q=0.737 | 0.103 [-0.767, 0.972]; q=0.790 | 0.031 [-0.704, 0.767]; q=0.925 | 0.008 [-0.742, 0.758]; q=0.980 |
| Red urchin biomass | 0.708 [0.259, 1.158]; q=0.025 | 0.607 [0.119, 1.095]; q=0.068 | 0.245 [-0.549, 1.038]; q=0.639 | 0.201 [-0.579, 0.981]; q=0.730 |
| Lobster density | 0.684 [-0.255, 1.623]; q=0.237 | 0.737 [-0.233, 1.707]; q=0.213 | 0.684 [-0.255, 1.623]; q=0.394 | 0.737 [-0.233, 1.707]; q=0.354 |
| Lobster biomass | 0.199 [-1.315, 1.713]; q=0.737 | 0.199 [-1.315, 1.713]; q=0.766 | 0.199 [-1.315, 1.713]; q=0.766 | 0.199 [-1.315, 1.713]; q=0.766 |

Added histories substantially reduce the red-urchin biomass contrast (0.708 → 0.245 on the old roster plus additions). Removing Gull KFM alone changes it less (0.708 → 0.607). For purple-urchin density, additions move the estimate toward zero (−0.430 → −0.262), but it was already uncertain without them. No new lobster-density histories are added: its numerical change (0.684 → 0.737) comes from removing Gull KFM, not additions. Changes in q across other membership arms can reflect the nine-response BH adjustment.

## 5. Before–after changes: membership alone

These are distinct from the post-protection comparisons above. The old eligible roster supports 41 before–after histories; researcher removals leave37, additions alone yield58, and both together yield54. Several responses have too few independent clusters/effective degrees of freedom for reliable intervals. No before–after lobster biomass history is eligible.

| Response | Old eligible roster | Removals only | Additions only | Final roster |
|---|---|---|---|---|
| Purple urchin density | -1.870; no reliable pooled interval | -1.870; no reliable pooled interval | -0.916 [-2.698, 0.866]; q=0.577 | -0.916 [-2.698, 0.866]; q=0.865 |
| Purple urchin biomass | Not estimable | Not estimable | -0.519; no reliable pooled interval | -0.519; no reliable pooled interval |
| Red urchin density | -0.142 [-1.964, 1.680]; q=1.000 | 0.080; no reliable pooled interval | 0.054 [-0.782, 0.889]; q=1.000 | 0.238 [-0.754, 1.231]; q=1.000 |
| Red urchin biomass | 1.181; no reliable pooled interval | 1.149; no reliable pooled interval | 1.175 [0.203, 2.146]; q=0.305 | 1.140; no reliable pooled interval |
| Lobster density | 0.978; no reliable pooled interval | 1.472; no reliable pooled interval | 0.978; no reliable pooled interval | 1.472; no reliable pooled interval |
| Lobster biomass | Not estimable | Not estimable | Not estimable | Not estimable |

Purple-urchin density had only one eligible legacy before–after history. Adding histories provides broader coverage but not a precise pooled decline. Purple-urchin biomass had no eligible legacy before–after histories; its four-history estimate is new, with inadequate independent information for reliable pooled inference. Removing Gull KFM reduces already sparse red-urchin biomass and lobster-density before–after replication. Missing robust intervals are not evidence of no biological effect.

## Source reconstruction check

C1–C3 in pooled.csv separate common-year input values from reconstructed survey support. For the six focal urchin/lobster responses, reconstructed source values/support do not change individual contrasts on the retained intersecting roster. KFM kelp reconstruction affects the kelp result and can indirectly change family-adjusted q. Source changes were not silently attributed to site additions.

## Checks, limits and reproduction

Current Table2 and the proposed endpoint were independently matched to saved outputs. The R reconstruction checks all143 legacy effects and SEs and requires convergence of accepted nonlinear fits. Three warnings arose from unsuccessful nonlinear starting fits; accepted final fits converged and reproduced the saved estimates. New common estimands are full-rank linear fits. This is conditional on existing legacy model selection, not a new model-selection exercise or proof of causal validity.

Independent Python checks reconstruct common individual contrasts/HAC covariances, pooled GLS point estimates, interval arithmetic and nine-response BH correction. They assert unchanged cohort membership through the B model steps, unchanged per-history values across membership arms, and original observation counts for the same-observation estimator comparison. Existing protected-file hashes are checked. See verification.json and endpoint_verification.csv. Sampling/size uncertainty and causal identification limitations remain; this comparison does not resolve them.

```sh
Rscript code/tools/attribute_analysis_changes.R
python3 code/tools/check_change_attribution.py
python3 code/tools/report_change_attribution.py
```

## Figure captions

- pooling_uncertainty.png/pdf: same143 legacy histories; points and95%intervals after each pooling/uncertainty step. The current joint model is followed by separate-response pooling, CR2 inference and calendarHAC3 individual covariance. The axis shows the legacy mixed contrast; zero is no contrast.
- definition_window.png/pdf: fixed intersecting legacy roster. Effect definition, designation-year treatment, follow-up window and removal of Before observations change sequentially. Points/intervals use the same pooling and robust covariance convention; their biological interpretation changes from the legacy mixture to mean post-protection differences.
- membership.png/pdf: fixed new estimators and individual inputs, with old eligible roster, removals alone, additions alone and both. Points are log protected/reference ratios; bars are95%CR2 intervals. Open points would indicate unsupported pooled intervals. No embedded legends, titles or caption prose.

Machine-readable effects.csv, observations.csv and pooled.csv retain all nine responses and all comparison steps. The paper has not been updated.
