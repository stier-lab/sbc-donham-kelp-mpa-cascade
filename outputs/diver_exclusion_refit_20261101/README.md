# Historical diver exclusions: review-only rerun

The paper and canonical analysis outputs have not been updated. This comparison keeps the flat-baseline primary model, trend sensitivity, calendar-year window, dependence groups and uncertainty method from the preceding common-estimand review. Membership does not depend on slope, sign or statistical significance.

## What was recovered

Audited all 22 omitted focal diver biomass/density histories already in the harmonized files. Twenty-one contribute to the post-protection comparison; 17 also contribute to the before–after comparison. Harris Point lobster density is zero on both sides throughout the primary period and cannot supply a sampling variance for weighted pooling. It stays in the audit; no arbitrary variance is assigned. Matlahuayl red/purple urchin biomass and Santa Barbara Island red/purple urchin biomass contribute post-only because baseline coverage is insufficient.

The roster now has 169 candidate histories (151 diver and 18 satellite). The primary fits contain 62 change histories and 166 post-protection histories, versus 45 and 145 previously. These are response histories, not independent MPAs.

## Source consistency

Restoring Harris Point kelp revealed that saved KFM kelp ratios and the current raw biomass export use different survey support. The raw export explicitly retains the modern quadrat protocol; the ratio table retains earlier/different-area observations. Rebuilding from current raw annual biomass also changes the adaptive zero correction. Rather than repair only the newly restored history, the consistent_kfm companion applies exactly the same reconstruction to all six KFM kelp biomass histories. Five previously retained histories change; all non-KFM-kelp annual data are untouched. This is a source/protocol correction, separately identifiable from adding excluded responses. No new biomass coefficients were fitted.

The parent directory preserves the additions-only comparison: all 1,333 previous individual fits reproduce unchanged. The consistent_kfm directory is the preferred internally consistent review result. Its kfm_protocol_audit.csv gives removed years, remaining baseline support and correction constants. Modern-protocol coverage still meets the three-year baseline rule where previously eligible, but Santa Barbara Island has only three prior years.

## Changes in the primary pooled results

Estimates are in log-ratio units. Before–after change measures a change in the protected/reference difference; post-protection differences alone do not establish an effect of protection. Intervals are 95% CR2 intervals. Missing intervals indicate inadequate independent information, not zero uncertainty.

| Comparison | Response | Previous | Added histories | Consistent KFM | 95% interval | Adjusted p |
|---|---|---:|---:|---:|---|---:|
| Change | Red urchin biomass | 1.181 | 1.175 | 1.175 | 0.203 to 2.146 | 0.305 |
| Change | Sheephead biomass | 0.237 | 0.237 | 0.237 | Insufficient information | Unavailable |
| Change | Kelp biomass | 0.309 | 0.305 | 0.328 | -0.129 to 0.785 | 0.330 |
| Change | Purple urchin biomass | Unavailable | -0.597 | -0.597 | -1.514 to 0.320 | 0.330 |
| Change | Red urchin density | -0.142 | 0.054 | 0.054 | -0.782 to 0.889 | 1.000 |
| Change | Sheephead density | 0.521 | 0.277 | 0.277 | -0.561 to 1.116 | 0.717 |
| Change | Purple urchin density | -1.870 | -1.137 | -1.137 | -2.542 to 0.268 | 0.330 |
| Change | Lobster density | 0.978 | 0.978 | 0.978 | Insufficient information | Unavailable |
| Post-protection | Red urchin biomass | 0.708 | 0.245 | 0.245 | -0.549 to 1.038 | 0.639 |
| Post-protection | Sheephead biomass | 1.166 | 1.166 | 1.166 | 0.493 to 1.840 | 0.024 |
| Post-protection | Kelp biomass | 0.326 | 0.261 | 0.285 | -0.255 to 0.825 | 0.426 |
| Post-protection | Purple urchin biomass | -0.343 | -0.422 | -0.422 | -1.068 to 0.225 | 0.383 |
| Post-protection | Lobster biomass | 0.199 | 0.199 | 0.199 | -1.315 to 1.713 | 0.766 |
| Post-protection | Red urchin density | 0.127 | 0.031 | 0.031 | -0.704 to 0.767 | 0.925 |
| Post-protection | Sheephead density | 0.621 | 0.630 | 0.630 | 0.228 to 1.032 | 0.024 |
| Post-protection | Purple urchin density | -0.430 | -0.304 | -0.304 | -0.916 to 0.309 | 0.426 |
| Post-protection | Lobster density | 0.684 | 0.684 | 0.684 | -0.255 to 1.623 | 0.383 |

Sheephead biomass and density remain higher inside protected sites after designation. Restoring the omitted red-urchin biomass histories removes the earlier clear positive post-protection difference. Purple-urchin before–after estimates can now be pooled across several sites, but their uncertainty still includes no change. The before–after results do not establish a complete predator–urchin–kelp cascade. A positive red-urchin biomass change has few effective degrees of freedom and does not pass the multiple-comparison adjustment.

## Sensitivity checks

The model suite includes the full observed post period, a ten-year baseline, original timing, HAC lags 1/5, trend-adjusted baseline, exclusion of Gull Island KFM, and omission of Harris Point lobster. All are exported in pooled.csv. Harris lobster omission is identical to primary because its constant primary series is already non-estimable. Gull omission reduces already limited independent information for some urchin estimates. With consistent KFM inputs, the trend-adjusted pooled kelp estimate is -0.382 (95% interval -1.652 to 0.887), versus +0.328 (-0.129 to 0.785) under the flat baseline. Sparse modern-protocol baselines at Gull Island and Santa Barbara Island produce large trend extrapolations; this sensitivity does not establish a negative kelp response. Both models leave the protection-related kelp change unresolved. The full-period kelp result must not be substituted for the first-nine-year primary based on its significance.

## What remains unresolved

See UPSTREAM_AUDIT.md for the broader candidate-site audit. Anacapa Island SMCA, Painted Cave and the Point Dume sites are not defensibly ruled out by the inherited broad filters alone; they need response-level reconstruction from source surveys before inclusion. Blue Cavern Onshore/Long Point also need review beyond their curated sheephead branch. Other sites lack verified reference assignments or matched coverage. These unresolved upstream cases prevent calling every historical exclusion resolved.

## Figure captions

pooled_changes.png compares the preceding common-estimand review (gray, upper point) with the current review (blue, lower point). Filled points have 95% robust intervals; hollow points lack sufficient independent information for an interval. The old purple-urchin density change was a single-site descriptive estimate. No before–after lobster biomass estimate is available. Left: change in protected/reference difference; right: observed post-protection difference. In consistent_kfm/pooled_changes.png, blue includes the uniform KFM reconstruction. Points are not necessarily causal effects.

## Verification

Independent Python calculations reproduce individual least-squares contrasts, calendar-gap HAC errors, pooled generalized least-squares coefficients, robust interval/p arithmetic and multiple-comparison adjustments. Protected-file hashes verify that the 487 canonical/reference artifacts remain unchanged. See verification.json in each output directory. Figures visually inspected separately. No sends, commits or pushes.
