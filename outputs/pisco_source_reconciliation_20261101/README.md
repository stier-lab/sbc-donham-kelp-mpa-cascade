# Three PISCO histories: source reconciliation

November 1, 2026. Review outputs only; manuscript, supplement, canonical results and the preceding review are unchanged.

## Finding

All three discrepancies are explained and the old annual values reproduce exactly. No missing measurement was recovered: the archived copies are byte-identical, the alternate size fields are blank, and the notes provide no numerical replacements. Explicit unknown measurements had been converted to zero by the inherited processing. The preceding audit excluded an entire site-year when it encountered a missing measurement, which explains why its reconstruction differed from those legacy values.

| History | Source record | What is unknown | Legacy treatment |
|---|---|---|---|
| Point Vicente kelp biomass | UNDERWATER_ARCH, July 13, 2023, INNER transect 2 | Two plants counted; fronds per plant missing | Missing count × size product became zero fronds |
| Scorpion red-urchin density | SCI_CAVERN_POINT_W, September 16, 2013, INNER transect 2 | Urchin count blank; notes blank | Blank count became zero |
| South Point sheephead biomass | SRI_CLUSTER_POINT_S, November 1, 2005, INNER transect 3, BOT | Three fish counted; all length fields blank; female recorded | Missing biomass became zero |

A blank count is not verified absence. The positive kelp and fish counts establish that organisms were present; their unknown sizes cannot justify zero biomass. Neighboring transects are different samples and cannot supply the missing measurements. Recovering those values requires original field records or contributor confirmation.

## Treatments tested

1. Reproduce the previous review. This preserves the legacy values in these three histories, which had been blocked from the preceding fixed-footprint correction.
2. Omit the response-specific incomplete zone in the affected year, retaining other measured zones and stations. No replacement count or size is imputed. This is a complete-zone sensitivity: habitat coverage can still vary, and the missing observations may not be random.
3. Test the intersection of valid station-zone combinations observed in every modeled year. This is a mechanical diagnostic, not automatic evidence of ecological comparability. Point Vicente retains 18 combinations at six stations; Scorpion retains 12 combinations at four stations. South Point has only two continuously sampled combinations at one station, so no protected/reference comparison is possible. The third scenario therefore retains the complete-zone treatment for South Point and explicitly flags its changing footprint; it is not a fully fixed-footprint solution for all three.

Every stage retains the same fitted years, MPA histories, species, designation clocks, researcher exclusions, additive correction constants and uncertainty method. Selection uses measurement validity and survey coverage, not outcome direction. The full set of 169 post-protection histories and the earlier proposed 140-history subset are both rerun; the 55 before–after histories are unchanged.

## Supported corrections, with the DEEP-only result withheld

These remain partial sensitivities. Apply the complete-zone treatment to all three and the fixed multi-zone treatment to Scorpion only; do not substitute Point Vicente’s DEEP-only series. South Point still lacks continuous station-zone coverage. Full 169-history pooled annual changes and 95% intervals:

| Response | Previous review | Supported corrections |
|---|---:|---:|
| Kelp biomass | +8.89% (+1.14, +17.23) | +8.89% (+1.14, +17.23) |
| Red-urchin density | -6.84% (-12.99, -0.25) | -7.02% (-13.21, -0.39) |
| Sheephead biomass | -0.06% (-6.49, +6.82) | -0.07% (-6.51, +6.80) |

## Individual trajectories

Annual percentage change in the protected/reference ratio. These are relative trajectories, not absolute population growth or effects uniquely attributable to protection. They should not be read as statistically supported individual effects without their uncertainty; log slopes, HAC standard errors and sample sizes are in history_comparison.csv.

| History | Previous review | Incomplete zones omitted | Mechanical fixed-zone diagnostic |
|---|---:|---:|---:|
| Point Vicente · kelp biomass | +27.45% | +27.45% | -27.58% |
| Scorpion · red-urchin density | -4.07% | -4.07% | -7.28% |
| South Point · sheephead biomass | -1.63% | -1.81% | -1.81% |

The small missing-value repair and the large coverage change are distinct. Omitting the incomplete Point Vicente zone in 2023 leaves the numerical estimate unchanged because the remaining observed zones at that station also have zero kelp biomass. This does not establish that the unmeasured plants had zero biomass. Restricting the entire series to continuously observed station-zone combinations produces the large slope change. The only continuously retained protected zone is POINT_VICENTE_W DEEP. In 2013, INNER/MID/OUTER were sampled in April and filtered out; only June DEEP surveys at 21–22 m remained. DEEP has zero kelp in every year from 2017 onward. LONG_POINT_E also lacks 2013 coverage. The mechanical intersection therefore changes the target to deep habitat. This reversal is not an acceptable corrected estimate of MPA-wide kelp recovery. See point_vicente_2013_depth_audit.csv.

## Pooled results: all 169 histories

Annual percentage change with unadjusted 95% CR2/Satterthwaite intervals:

| Response | Previous review | Incomplete zones omitted | Mechanical fixed-zone diagnostic |
|---|---:|---:|---:|
| Kelp biomass | +8.89% (+1.14, +17.23) | +8.89% (+1.14, +17.23) | +8.30% (+0.10, +17.16) |
| Red-urchin density | -6.84% (-12.99, -0.25) | -6.84% (-12.99, -0.25) | -7.02% (-13.21, -0.39) |
| Sheephead biomass | -0.06% (-6.49, +6.82) | -0.07% (-6.51, +6.80) | -0.07% (-6.51, +6.80) |

Lobster and purple-urchin estimates are unchanged in this reconciliation; adjusted q-values can still move when other tests change. After fixed-footprint treatment where feasible, full-cohort kelp and red-urchin q-values are both 0.142. In the earlier 140-history subset, kelp changes from +9.07% to +8.36% per year; its interval changes from +1.07–17.70% to −0.15–17.60%, and q from 0.064 to 0.121. Sheephead remains uncertain. These findings do not justify selecting the footprint according to significance.

## Companion-response check

The missing Scorpion count also entered the original red-urchin biomass cache: the INNER-zone count total is 120 divided across two transects, although one transect count is unknown. The preceding fixed-footprint biomass sensitivity already excludes SCI_CAVERN_POINT_W entirely, so its reviewed biomass estimate requires no further patch. This issue still needs propagation to canonical density and biomass processing before paper integration. South Point sheephead counts remain valid when lengths are unknown; this audit does not remove those fish from density. There is no separate modeled PISCO kelp-density response.

## Recommendation and remaining confirmation

Keep the explicit distinction between absent taxa and unknown measurements. The source discrepancy is resolved, but the original measurements remain unavailable. Retain the complete-zone and Scorpion fixed-footprint comparisons as sensitivities pending a decision about the spatial population being estimated. Reject Point Vicente’s DEEP-only intersection as a general kelp-recovery correction. The Supported_corrections scenario retains complete-zone Point Vicente and South Point estimates and applies the stable, multi-zone Scorpion comparison. None of these recover the missing measurements. South Point still lacks a continuous footprint; do not describe it as corrected for changing stations/zones. Do not replace the three unknown fields with zeros or infer their values from other transects.

If measurement recovery is pursued, the precise requests are the three source records listed above: a red-urchin count, fronds for two plants, and lengths for three female sheephead. No clarification emails were sent in this task.

## Figure caption

three_histories.png/pdf: annual protected/reference ratios for the three affected histories. Gray: previous review; blue: omit incomplete zones in affected years; orange: mechanical fixed-zone diagnostic, including the rejected DEEP-only Point Vicente subset, with the complete-zone estimate retained at South Point. Blue and orange overlap at South Point. Ratios retain original additive constants. Lines connect observed years without imputing intervening surveys; these are descriptive ratios, not fitted curves or confidence bands. All records shown are post-protection.

## Reproduce and verification

Run code/tools/reconcile_pisco_sources.py, Rscript code/tools/refit_reconciled_pisco.R, code/tools/check_reconciled_pisco.py and code/tools/report_reconciled_pisco.py. The checker independently verifies 100 reconstructed role-year values, 2,012 individual contrasts/HAC variances, 236 pooled points, BH/interval arithmetic, fixed-footprint membership, unchanged unrelated histories and 55 before–after fits, source/archive hashes, the preceding review manifest and 487 protected files. Source hashes, exact anomalies, coverage decisions, annual values and model rows are saved alongside this report.
