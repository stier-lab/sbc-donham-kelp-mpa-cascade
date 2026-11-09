# Remaining site limitations

Reviewed November 3, 2026: all 25 protected areas and 46 area–program groups in the current decision table. This is a readiness review, not a change to inclusion, models or paper text. [Site-by-program actions](site_program_actions.csv) cover every group; [named-station contributions](named_station_contributions.csv) identify the stations implicated in the remaining reference questions.

Most flags are work we can finish ourselves or limitations to disclose. The main unresolved contributor questions concern KFM reference suitability and permission. Unknown metadata do not automatically make a comparison unusable, but known ecological concerns require a response-specific decision before finalizing the affected results.

## Work supported by existing files

- **KFM reference influence.** Wyckoff Ledge supplies supported values to seven Harris Point histories; Miracle Mile also contributes to Harris Point sheephead density. Webster’s Arch contributes to Santa Barbara Island sheephead density in 17 years, 2005–2022. Rodes Reef contributes to no currently eligible fitted history. The [completed reference tests](../kfm_reference_sensitivity_20261103/index.html) now separate common-year changes from lost coverage. Combined omission loses seven Harris histories and 217 supported comparison-years. These tests cannot prove ecological equivalence; Scott’s ecological recommendation remains pending.
- **PISCO depth, timing and station coverage.** Original fish/swath deliveries contain measured depths and dates. We can summarize the retained observations, inspect residual depth differences and review the fitted constant-station/zone alternatives. Keep the approved all-month rule; no new seasonal cutoff is implied. At Abalone Cove, Point Vicente and Malibu, Dan’s runoff, sedimentation and turbidity concerns remain possible alternatives that depth matching alone cannot eliminate.
- **Point Dume/LECHUZA.** Both conflicting coordinate candidates are recorded. We can assess whether choosing either would affect a distance or boundary claim using available mapping evidence. Source confirmation is needed if the correct identity cannot be recovered and the distinction affects treatment assignment. Do not pick a coordinate by row order. Existing station labels do not themselves resolve the conflict.
- **Satellite timing and reference weighting.** The cached NetCDF is present and contains 152 quarterly observations, including year, quarter, biomass and environmental fields. Saved corrected masks and source deliveries are present. We can audit quarterly support and sensitivity to annual aggregation/reference weighting without reconstructing stochastic historical clustering. This has not been rerun in this review.
- **Uncertainty and release work.** Review the fitted boundary/influence and sparse-size results, run the outstanding ratio-offset check, then regenerate pathway/resilience analyses and final figures/tables on the adopted cohort. Existing files support this work; new coauthor data are not a prerequisite. Missing detection, coefficient covariance and within-group size variation cannot be manufactured by another bootstrap.

## Contributor input that remains useful or necessary

| Issue | Person | Remaining question | When it matters |
|---|---|---|---|
| Harris Point KFM | Scott Gabara | Given the retained Hare Rock–Wyckoff/Miracle Mile comparisons and their response-specific influence, does he recommend excluding the comparison, retaining it with a qualification, or using another established reference? | Settle the treatment of an explicit ecological concern before finalizing affected results. |
| Santa Barbara Island KFM fish | Scott Gabara | Does his recommendation concerning offshore/pinnacle Webster’s Arch apply to sheephead density as well as benthic responses? | Webster’s Arch is already omitted from benthic responses but still contributes to fish. |
| KFM use and redistribution | Scott Gabara / park data manager | Is the analysis covered by a permit, and which source or derived records may be archived? | Confirm before submission/public redistribution. This does not prevent private analysis work. |
| LECHUZA | PISCO station-metadata custodian, routed through the project collaborators | Which coordinate belongs to the surveyed station, if source mapping cannot settle it? | Required for a claim or selection depending on precise distance or boundary membership. |
| KFM/LTER survey depths and habitat | Program metadata custodians | Are station/transect depth or habitat records available beyond the count deliveries? | Needed to establish stronger physical matching; otherwise disclose missing metadata and avoid that claim. |
| KFM sheephead size footprint | Scott Gabara | Do the roving size surveys represent the same habitats and fish size ranges as the transect counts? | Relevant if promoting restricted biomass beyond a qualified supporting analysis. He has already confirmed one roving diver and no lobster lengths; do not repeat those questions. |
| LTER provider biomass formula | Li | Why do the adult values imply the July coefficient while the linked protocol says August? | Only if replacing our raw adult-stipe, actual-month calculation with the provider aggregate. Current adult source observations through 2023 match. |

Scott’s November 2 reference email was reread. A targeted current mailbox search found no newer reference reply, and the latest permission reply found was October 29, still awaiting park confirmation. No messages were sent. These searches do not certify that no relevant message exists under another sender or subject.

## Disclose or retain missing; do not hold all comparisons open

- **Naples PISCO:** Jenn has already explained diver expertise, community ordination and the lack of a close match. The underlying ordination is optional unless we claim independent quantitative validation. Do not ask her to repeat the rationale.
- **Campus Point/Naples LTER:** Emily has described consultation with Bob/Li. Exact station-specific reasoning remains incomplete. Retain the documented assignments and qualifications; ask for existing station metadata only if needed for a stronger claim.
- **Missing measurements:** Underwater Arch 2023 kelp stipes and South Point 2005 sheephead lengths remain missing. Little Dume 2011 has no station value in the current annual estimates, so recovering its missing stipes is not a current numerical blocker. Preserve observed counts and response-specific missingness; recovery would require another original record.
- **KFM predator biomass:** there are no lobster lengths in the checked KFM sources. The restricted sheephead review remains 52 supported comparison-years; five omitted years were searched without recovering the necessary lengths. Full-history density can proceed. The possible 55-year visit-subset calculation is a different sensitivity, not recovered data.
- **LTER group means and sparse samples:** individual size distributions, unobserved count–size association and detection uncertainty are not recoverable from these files. Keep the conditional uncertainty language.
- **Historical dates and satellite selection:** network dates and present-day restrictions are not complete historical station-boundary records. Saved satellite masks reproduce delivered measurements but do not prove unique historical clustering or ecological equivalence. A pre-2012 selection window is not entirely before island protection. Disclose these limits; seek additional evidence only for claims that require it.
- **Out-of-scope records:** Trailer/Willows, unused KFM 2023 groups, new 2024–2026 records and already excluded histories do not need resolution to complete the current cohort. Existing exclusion decisions remain intact.

## Recommended order

- Review the completed KFM reference tests and the unsent named-station note for Scott.
- Finish local depth/timing descriptions, the offset check and diagnostic review; keep previously approved rules fixed.
- Resolve KFM reference scope before finalizing those results. Keep permission confirmation as a submission/archive requirement.
- Adopt the reviewed cohort, then regenerate dependent analyses, figures and manuscript concordance. Carry unresolved measurement/habitat limits into the methods and discussion without presenting unsupported biomass or stronger causal claims.

Evidence: current integrated station/annual/model tables; `outputs/notes_email_audit_20261102/`; `outputs/kfm_missing_size_recovery_20261103/`; `outputs/fish_lter_reconciliation_20261102/`; `outputs/depth_reference_policy_20261101/coordinate_conflicts.csv`; current uncertainty and freeze tables. Older reports supply source evidence only; their obsolete proposed coverage rules are not adopted. Run `python3 code/tools/review_current_site_limitations.py` to reproduce the membership/contribution tables. Input hashes are in `verification.json`.
