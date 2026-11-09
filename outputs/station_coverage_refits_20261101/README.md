# Station coverage and reference influence

November 1, 2026. PISCO sheephead diagnostic only; paper and canonical outputs unchanged. This follows the approved CRANE/repeat/complete-pass rebuild and the user's coverage decisions.

## What was compared

All 34 sheephead response histories were refitted over their actual post-protection years. First, replace the available reference average with stations present in every approved year; hold protected means, offsets and fitted years constant. This is possible for 30 histories; Abalone Cove and Santa Barbara Island have no persistent reference station over the full approved record. The approved record includes before-years where available; the fitted recovery slopes use only post-years. No station was chosen for its result.

Second, omit each reference station in turn and compare both treatments on exactly the same remaining years. Eighty-two such comparisons show influence, not justification for removing a station.

Third, fit weighted log-station observations with year intercepts and a protected-by-calendar-time interaction, with and without station intercepts on identical rows. Each role/year has equal total weight. This checks whether persistent station differences explain an apparent trend when membership changes. It is a separate estimand: an average of log station values is not the log of an arithmetic mean. Differences between these models and the approved ratio model therefore cannot be attributed solely to station adjustment. Only the adjusted-minus-unadjusted contrast within this model isolates station intercept adjustment. Constant station intercepts do not remove habitat-specific trends or seasonal/depth changes.

## Findings

Point Dume SMCA is most sensitive to station-composition adjustment. Using Nicholas Canyon West throughout the same five observed post-protection years changes sheephead biomass slope from −0.0502 to −0.1733 log units/year and density from +0.0250 to −0.0870. All four 95% intervals span zero. On the separate log-station scale, adding station intercepts changes biomass +0.0784→−0.1442 and density +0.1076→−0.0641. These diagnostics support treating its recovery estimate as sensitive to station coverage; they do not establish a protection effect or justify exclusion by direction.

The persistent reference subset changes Anacapa biomass −0.0023→−0.0350 and density −0.0603→−0.0640; Scorpion biomass −0.0065→−0.0037 and density −0.0066→−0.0080; South Point biomass −0.0188→−0.0156 and density −0.0075→+0.0035. These use identical years within each response. Review the intervals in effects.csv, not just signs.

Reference identity also matters even when membership is constant. Omitting Cuyler West at Harris Point changes biomass slope −0.1349→+0.0898 and density −0.0564→+0.0499 on the same 17 post-years. Omitting Bunker Point reduces Abalone's biomass slope +0.3423→+0.2126 and density +0.1744→+0.0878 on the same 11 post-years. These are influence checks; no researcher-approved reference is removed. Reference-omission intervals and exact years are preserved in reference_omission.csv. They are many exploratory checks, not multiplicity-adjusted claims.

## Shared zones: correction to the broader screen

The preceding role-level screen found shared zones in every fish year. Checking every contributing station is stricter: Abalone has no zone shared across all stations in 2010 or 2014, so a full-history constant common zone is impossible under unchanged membership. Those are coverage findings, not bad-data labels. Other assignments have at least one zone shared at every station throughout, but often a narrower set: Anacapa, Point Vicente and South Point reduce to OUTER; Santa Barbara Island to INNER. Other station-zone sets are in station_zone_summary.csv. This does not endorse an OUTER-only primary analysis: zone definitions and measured depth must be verified, and narrowing zones changes the sampled habitat. No zone-restricted outcome refit is included here.

## Recommendation

Retain the established reference assignments and report these sensitivities. Do not require every historical reference every year. Prioritize a common station/zone prediction population and explicit station adjustment, especially at Point Dume, while checking that seasonal and habitat changes can be separated from recovery. Keep same-year persistent-reference and leave-one-reference checks as diagnostics. Do not select the version with the preferred sign. Before any primary replacement, resolve zone comparability, baseline station continuity and identification, then refit other taxa and pooled/pathway/resilience analyses together.

## Validation and limits

All 34 approved recovery slopes and HAC3 SEs reproduce to <1e-8. R independently reproduces all 68 weighted log-station coefficients to <2e-14. The 162 reported fits include the approved, persistent-reference, matched-year control and two log-station models. All 487 protected manuscript/result files remain unchanged. Year-summed Bartlett HAC3 scores preserve calendar gaps in the exploratory station models; their year-based t reference and finite-sample correction are approximate, not CR2 or validated small-cluster inference. The existing ratio fits retain their original HAC3/t treatment. No regional pooling, before–after refit or claim of causal identification is made. KFM/LTER, satellite, other PISCO taxa and expanded size support remain outside this numerical pass. No emails, commit or push.
