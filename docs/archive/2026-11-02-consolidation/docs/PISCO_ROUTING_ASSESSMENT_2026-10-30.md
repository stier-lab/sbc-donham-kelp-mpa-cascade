# PISCO model assignment assessment

30 October 2026

## Recommendation

Retain the current primary analysis while explicitly reporting the alternative-rule results. Treat PISCO control-impact assignment as a program-wide analytical choice, not as proof that all PISCO baselines are inadequate. This assessment does not establish a final, uniform before/after eligibility rule. The fish exceptions have now been reconstructed with coverage-selected fixed stations (module 59; see [the fixed-station assessment](PISCO_FIXED_STATIONS_2026-10-30.md)). Their pooled sheephead result persists, but historical baseline thresholds and designation-year treatment still require a defensible decision before changing the primary rule.

## Coverage and the remaining exception

The source audit separates fish and benthic surveys and intersects years at the same protected and reference stations. No pair has three baseline years across all retained benthic responses. Fish pairs at Campus Point, Naples and Point Vicente have four baseline and four later years if designation-year surveys count toward baseline. Only Point Vicente has four strictly prior years. Removing fish station-years with invalid counts or missing positive-count lengths does not change these thresholds. Station pairs share sites and are not independent replications.

The initial source audit provided coverage bounds. Module 59 now reconstructs fixed-station biomass and density at the three fish exceptions. Switching model labels on existing annual means leaves their changing station composition intact. A history-level minimum-year rule therefore does not implement the station-continuity criterion. A uniform rule requiring all retained species at each MPA to qualify instead reroutes no PISCO histories, but conflates fish and benthic coverage and is not supported as the final historical criterion.

## Current fixed-input model comparison

All scenarios retain the 145 histories and their annual inputs. Models use the current calendar-year policy and observed contrast endpoints. Inclusive baseline includes the designation year; strict baseline excludes it. Candidate screening and fitting are in module 52.

| Baseline | Minimum years | Histories reassigned | Lobster density adjusted p | Lobster–purple-urchin biomass slope and nominal 95% interval |
|---|---:|---:|---:|---|
| Inclusive | 3 | 20 | 0.013 | -0.560 (-1.581, 0.462) |
| Inclusive | 4 | 11 | 0.042 | -0.560 (-1.581, 0.462) |
| Inclusive | 5 | 5 | 0.095 | -0.555 (-1.032, -0.078) |
| Strict | 3 | 11 | 0.048 | -0.558 (-1.534, 0.419) |
| Strict | 4 | 5 | 0.095 | -0.555 (-1.032, -0.078) |
| Strict | 5 | 2 | 0.090 | -0.555 (-1.032, -0.078) |

All nine pooled directions persist across these six alternatives. Only lobster-density adjusted support changes, in the inclusive three- and four-year and strict three-year alternatives. Those same alternatives widen the four-MPA lobster–urchin biomass interval across zero. These are primary-variance comparisons; the exported calendar-HAC variants provide a separate uncertainty check. Neither stable directions nor altered significance validates a design choice.

## What remains before changing the primary rule

The fixed-station reconstruction and matched-year model comparisons are complete in module 59. See the [design, coverage and findings](PISCO_FIXED_STATIONS_2026-10-30.md). Fixed stations preserve the pooled sheephead biomass increase and leave sheephead–urchin relationships unresolved. Some marginal pooled lobster/urchin classifications change through the joint variance model. Neither stability nor changed significance determines the eligibility rule.

Resolve three versus four baseline years and designation-year treatment, then choose a uniform prospective rule from design evidence. Passing a baseline-trend screen does not establish habitat equivalence or parallel counterfactual trends. Historical thresholds remain imperfectly documented; substantive historical decisions require contributor confirmation. Primary routing has not been changed.

## Reporting and verification

The current main text already discloses the program-wide choice and points to S2.3. Added the missing alternative-model results to S2.3, cross-referenced existing coverage details in S1.1, and shortened the adjacent KFM station check without changing its findings and indexed the existing module and result tables. Primary estimates and the main manuscript are unchanged. The supplement is generated from current exported results with assertions for the stated directions and classification changes.

Checks: code/tests/check_routing_clock_sensitivity.py reproduces 142 linear/mean/step fits, checks all 14 fixed-cohort scenarios and current hashes. code/tests/check_selection_criteria.py verifies coverage inputs and claims. Private correspondence evidence is retained only in the manuscript repository.
