# Approved PISCO fish measurement rules: review rerun

November1,2026. Adrian approved retaining CRANE, averaging repeats and omitting incomplete passes from the combined response. The rules are implemented in this reproducible review branch; the manuscript, supplement and canonical result files are unchanged. frozen_rules.json records the decision.

## Scope and implementation

Rebuilt sheephead counts and biomass directly from the supplied2024 fish file. Completed BOT plus MID or CNMD units are retained, including supported CRANE2004. Standalone CAN remains outside this response. Observed NO_ORG passes remain zeros; unrecorded opposite passes are not filled with zeros. Repeated dates are averaged within site/calendar-year/zone/transect number, followed by the existing equal-transect, equal-zone and equal-site hierarchy. Transect numbers are sampling labels, not claims of permanent locations. Missing positive-count lengths keep biomass unknown. No source rows are deleted.

Researcher/ecological exclusions, response-specific station/year membership, allometry, ratio offsets, actual calendar years, observed follow-up, baseline eligibility and model uncertainty are held fixed. The other proposed depth/date/duration rules have not been adopted by this instruction. The analysis retains175 candidate histories. These are not a declaration that every history is ready for primary-model eligibility.

## Coverage

The full annual review input changes from3,228 to3,225 rows. All175 histories remain. Of420 PISCO sheephead history-year rows,417 remain; the missing rows are SBI2005 density/biomass (incomplete passes) and South Point2005 biomass (previously identified missing lengths). CRANE2004 remains in Anacapa2003, SBI and South Point. The response-specific missingness rules retain South Point2005 density. No MPA is removed.

## What changed

Six PISCO sheephead recovery comparisons change: Abalone Cove, Gull Island, Point Vicente, Scorpion, South Point and Santa Barbara Island. Recovery signs do not reverse. Point Vicente before–after biomass change rises from0.3721 to0.4714 and density change from0.0342 to0.1496 log units. Abalone biomass recovery slope changes0.3514→0.3423 and density0.1781→0.1744; Point Vicente slopes change0.1222→0.1253 and0.0674→0.0700. All other taxa's fitted effects are unchanged.

The full pooled sheephead recovery slope is−0.0094 per year for biomass (95% interval−0.0870 to0.0681), compared with−0.0097 previously. Density is−0.0231 (−0.0628 to0.0165), compared with−0.0230. The intervals continue to include zero. These are regional recovery trends; they do not answer whether protected sites have greater average biomass.

Pooled before–after biomass has five histories but insufficient effective independent information for a robust interval; do not treat its point estimate as a new significance result. Two pooled fits required documented BFGS optimization fallbacks. Study-level and pooled tables preserve unavailable intervals rather than supplying a misleading precise answer.

## Separate processing and coverage

Total compares previous full data against the approved rules. Same_years applies both treatments to precisely the approved years. Coverage_only retains previous values and removes only the three unavailable history-year rows. This separates processing changes from the effect of losing early observations. Station-role membership is held fixed; no affected year is silently rescued by dropping an unavailable station and averaging the remainder.

## Model and validation boundaries

This continues the ongoing review analysis: post-protection recovery is a linear trend over actual after-years; before–after change uses a flat baseline plus after-period trend, with at least three verified prior years. Designation-year observations are retained in the ledger but omitted from these fits. Temporal uncertainty uses calendar-gap HAC3; pooled models use the existing MPA/shared-reference clustering and CR2. No sign/significance rule determines inclusion. These estimates are not a direct replacement for the paper's older mixture of model-specific effects until the remaining model/cohort decisions are finalized.

The preceding full-follow-up review's224 corresponding effects and SEs reproduce exactly. Across the total/same-year/coverage comparisons,1,344 effects and SEs were independently recomputed in Python;168 pooled fits completed. Non-PISCO-fish observations and effects are unchanged. All487 protected manuscript/result files and both preceding audits remain unchanged. Pathway/resilience models are not refitted in this measurement-policy rerun; those downstream checks are required before paper integration.

## Reproduction

Run code/tools/build_pisco_fish_rules.py, Rscript code/tools/refit_pisco_fish_rules.R, code/tools/report_pisco_fish_rules.py and code/tools/verify_pisco_fish_rules.py. The raw-source rebuild independently matches the prior policy proposal to floating-point precision. Inputs and outputs are hash-recorded; source records are untouched.
