# Reference selection and station continuity

This audit separates the ecological reasons for choosing references from the stations that actually entered annual averages. Confirmed selection explanations improve transparency; they do not establish that protected and reference reefs had equivalent habitat or would have followed parallel trajectories.

## Confirmed explanations

- PISCO Harris Point: Cuyler is nearby on the north side of San Miguel with similar habitat. Crook Point and Tyler Bight are on the south side and differ in habitat. Their exclusion is not evidence that their records were poor quality.
- PISCO short histories: supplied fish records contain Prince Island in 2003, Bay Point in 2005 and Ford Point in 2007 only. Supplied swath records contain Prince Island in 2003 and no Bay Point or Ford Point records. These statements describe the supplied files, not all surveys ever collected. Hamilton's recollection of Ford Point's ecological selection remains tentative.
- PISCO Yellowbanks: distinctive habitat and distance make it a weak local match to Gull; Scorpion has closer alternatives. An island-wide comparison would ask a different question.
- KFM Gull: Yellowbanks is the best available reference according to the program investigator, approximately 15 miles away. Keep Gull in the primary analysis and retain its omission sensitivity. Distance remains a limitation, and a detailed map review is pending. PISCO and KFM use different survey networks, so these explanations are not contradictory.
- LTER: reference choices followed consultation with local program investigators. The individual retained-site rationale has not been recovered. Arroyo Hondo's documented exclusion for essentially absent kelp does not validate the habitat equivalence of the remaining sites.

The private manuscript repository retains the dated correspondence and durable source-message identifiers. No correspondence or restricted raw KFM records is copied here.

## Intended continuity rule versus implemented membership

The intended before/after rule was to retain stations observed in both periods, separately where fish and benthic coverage differed, while allowing later stations in after-only comparisons. Historical comments in `pBACIPS_PISCO_V10.R` (lines 18–22, 1124–1128 and 1459 onward) support that intent. The current benthic and sheephead branches do not apply identical exclusions.

`code/audit_kfm_site_continuity.py` checks 37 KFM lobster/sheephead site histories within primary density histories. It uses observed zeros and the independently reconstructed site means from the size-feasibility audit. It finds six after-only stations in the two sheephead pBACIPS histories: Miracle Mile at Harris Point; Potato Pasture, Cavern Point, Devil's Peak Member, Little Scorpion and Pedro Reef at Scorpion. No after-only station occurs in the audited lobster pBACIPS histories. This does not establish continuity for kelp or urchin biomass, LTER or PISCO.

## Same-year sheephead sensitivity

Module 53 keeps only stations observed in both periods and compares them with the original membership on identical years, keeping the original offsets, models and time coding. All 74 original annual log ratios are reconstructed within 1e-10. Both comparisons retain 37 paired years. Scorpion 2023 is omitted from both arms because its longstanding reference was not sampled. No missing observation is replaced by zero.

| MPA | Original membership on matched years | Persistent stations on the same years |
|---|---:|---:|
| Harris Point | −1.172 (95% interval −2.210 to −0.134) | −1.116 (−2.209 to −0.023) |
| Scorpion | 1.574 (0.888 to 2.259) | 0.901 (0.033 to 1.770) |

Values are the existing selected log-contrast estimands, not abundance or a common fold effect. Intervals use the fitted model's residual degrees of freedom. Harris is an 11-unit fitted trend contrast; Scorpion is a before/after step. With calendar-gap HAC uncertainty, the persistent-site Scorpion interval includes zero. These conditional fits do not repeat model selection or validate habitat matching.

Station additions explain part of Scorpion's estimated sheephead increase. Harris changes little. Replacing these two contrasts in the existing synthesis leaves all nine pooled directions and all five multiple-test-supported responses unchanged. Pooled sheephead density remains uncertain (estimate −0.339, 95% interval −1.457 to 0.779). This is a sensitivity check, not a promotion of revised primary estimates.

Reproduce with `python3 code/audit_kfm_site_continuity.py` then `Rscript code/R/53_kfm_site_continuity.R`. Module 53 is included in both analysis runners. Source hashes, membership, annual comparisons, model intervals and pooled results are in `tables/*kfm_continuity*`. The count-site reconstruction is documented by `docs/KFM_PREDATOR_SIZE_FEASIBILITY_2026-10-29.md`.

## Still unresolved

The response-specific versus all-response three-before-sample rule, Campus Point routing, per-site LTER and remaining KFM habitat rationales, and suitability of separate KFM predator sizes remain distinct questions. KFM data-use permission and redistribution conditions are not yet confirmed; no raw data upload is authorized by these replies. Do not describe investigator recollection as a completed matching procedure.
