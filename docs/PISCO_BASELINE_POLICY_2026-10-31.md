# PISCO baseline rule

Adopted October 31, 2026: three matched survey years strictly before implementation for the primary analysis; four as a sensitivity. This is Adrian's explicit choice, not a claim that the ambiguous historical three/four-year wording established one threshold.

## Eligibility and annual inputs

Fish and benthic surveys are assessed separately, requiring shared coverage across the retained responses within each group. At least four later survey years are required. Stations are selected using coverage alone: maximize the number of stations, then common years, then alphabetical station identifiers. Keep the same stations and common years throughout each comparison. Exclude fish site-years with invalid counts or missing lengths for observed sheephead. The mainland MPAs took effect January 1, 2012; every checked 2012 fish survey occurred later that year, so 2012 is after implementation.

Six sheephead histories qualify: biomass and density at Campus Point, Naples and Point Vicente. Their baseline is 2009–2011. Campus Point retains 13 years (drops 2014–2015), Naples 14 (drops 2014), and Point Vicente 15 (drops 2004, 2007–2008). The previous annual data included those years under changing station mixtures; they were not omitted baseline observations. Six MPA-years correspond to 12 response-year rows. The retained 145 histories contain 2,815 annual rows, with 42 pBACIPS and 103 CI contrasts.

No benthic group qualifies under the shared-response criterion. Point Vicente's limiting response is lobster biomass, which begins in 2011. Several other benthic responses could support a response-specific alternative, but adopting that would change the eligibility unit and is not part of this decision. Raw 2009 lobster counts at Hawthorne Reef lack corresponding sizes; zero counts in 2010 cannot supply the missing third baseline year.

## Results and limits

The selected sheephead models are steps at Campus Point and Naples and linear changes at Point Vicente. Explicit baseline initialization and multiple nonlinear starts leave these winners unchanged. Three baseline years leave only one residual degree of freedom for the pretrend screen; passing it does not establish parallel trends or habitat equivalence.

The joint primary synthesis supports higher sheephead biomass, kelp, red-urchin biomass and lobster density, and lower purple-urchin biomass and density. The four-year alternative preserves all nine directions and the sheephead, kelp and purple-urchin findings. Red-urchin biomass and lobster density have adjusted p = 0.075 and 0.091, versus 0.041 each under the primary rule. Their data were unchanged; their uncertainty changes through jointly estimated variance components.

Only Point Vicente qualifies for the four-year before/after alternative. Campus Point and Naples revert to their original CI comparisons. This sensitivity changes membership and estimands as well as the minimum baseline length; it is not a same-cohort test of threshold alone. The module also retains fixed-station versus mixture and CI versus pBACIPS comparisons on identical years to distinguish these components. Neither sheephead–purple-urchin relationship is resolved. The four-MPA lobster-biomass pathway is unchanged.

## Reproduction and audit

- `code/R/59_pisco_fixed_stations.R` rebuilds station means from raw fish records, verifies the archived original mixture, selects coverage groups, and fits alternatives.
- `data/external/pisco_policy_20261031/` preserves explicit pre-policy comparator inputs. These are comparison evidence, not primary inputs.
- `code/promote_pisco_baseline.py` verifies the primary annual transformation; `--apply` rebuilds the two scoped harmonized files. Regenerate the reference registry and downstream results after applying.
- `data/harmonized/pisco_baseline_policy.json` records the rule, source date and chained input/output fingerprints.
- `code/R/pisco_baseline_policy.R` applies the primary routing; `pisco_nonlinear_candidates.R` records multiple starts and rejects invalid nonlinear candidates.
- `code/tools/check_pisco_primary_policy.py` verifies all other 139 primary histories are unchanged, exact annual values, strict periods, routes and row counts.
- `code/tests/test_pisco_fixed_stations.R` independently reconstructs linear/step fits, confirms complete station/year balance and coverage-only choices, and checks sensitivity claims.
- `tables/table_s_pisco_primary_*` records candidate models and optimization starts; `table_s_pisco_fixed_*` contains station membership and all sensitivity results.

Official implementation date: https://wildlife.ca.gov/Conservation/Marine/MPAs/Statistics. Historical correspondence is held in the private manuscript repository; it is not bundled with analysis data.
