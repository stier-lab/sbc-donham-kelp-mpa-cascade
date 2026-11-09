# PISCO repeat surveys and fish levels

November 1, 2026. Review analysis only. No manuscript, supplement, canonical result or researcher exclusion changed. This audit covers34 sheephead density/biomass histories at17 retained candidate PISCO comparisons; current-paper membership is flagged separately in comparison_summary.csv. The raw181 repeat-unit and85 unequal-level flags from the inventory are not181 affected comparisons.

## What the flags mean

Some sites were surveyed across several days, with different depth zones completed on different days. At other sites, two provider campaigns fall within one calendar year. Keep actual calendar years, but do not add repeated effort as though it were additional fish in one survey. Retain campaign and full date in the audit trail. Reused transect numbers are sampling labels, not proof of resurveying an identical permanent location. The proposed correction averages repeated observations within a transect/zone/year, then averages transects, zones and sites in the existing hierarchy. A day sampling one zone does not receive the same site-wide weight as a day covering several zones.

The retained comparison windows contain repeat-unit flags at Abalone Cove (2015/2019 reference stations), Point Vicente (2010/2015 protected stations), Point Dume SMR (2008) and Blue Cavern (2012). The latter two affect annual records but not the fitted recovery comparisons: Dume2008 precedes protection without a qualifying before baseline; Blue Cavern2012 is the designation year omitted from these review fits.

## What changes when repeats are averaged

All34 legacy annual histories reconstruct to numerical precision (maximum log-ratio difference below7e-14). Paired contrasts use identical years, roster, offsets, model formulas and uncertainty. The independently missing South Point2005 biomass value is withheld from both arms, not credited as a repeat-survey effect.

| MPA | Response | Quantity | Before | After |
| --- | --- | --- | --- | --- |
| Abalone Cove SMCA | Bio | Post_slope | 0.3514 | 0.3412 |
| Abalone Cove SMCA | Den | Post_slope | 0.1781 | 0.174 |
| Point Vicente SMCA | Bio | Change | 0.3721 | 0.4167 |
| Point Vicente SMCA | Bio | Post_slope | 0.1222 | 0.1252 |
| Point Vicente SMCA | Den | Change | 0.0342 | 0.0894 |
| Point Vicente SMCA | Den | Post_slope | 0.0674 | 0.0698 |

Recovery quantities are changes in log protected/reference ratio per calendar year. Change is the mean after-minus-before difference, with a flat baseline and post-period trend; it is not the manuscript's mixture of model-specific contrasts. These targeted refits use all observed after years and preserve before/after versus post-trend separation. The correction leaves Abalone Cove and Point Vicente recovery estimates positive. It does not change the remaining15 PISCO recovery point estimates.

Across the retained sheephead histories from all programs, the pooled recovery slope changes from-0.00991 to-0.00927 for biomass and from-0.02297 to-0.02309 for density. Both pooled intervals still span zero. This concerns regional recovery slopes, not the separate question of whether protected sites have higher average biomass. Five-history pooled before/after biomass has too few effective independent degrees of freedom for a robust confidence interval; no significance claim is made.

## Unequal levels are not all missing effort

The inspected units comprise11,482 BOT+MID,2,249 BOT+CNMD,164 BOT-only,42 MID-only and10 CNMD-only units. These are candidate-source records before restricting every comparison to its paired years. All2,249 BOT+CNMD units are INNER-zone records. CNMD means combined canopy/midwater in the integrated-data documentation, and the published PISCO methods describe one or two additional upper-water passes depending on depth and canopy. Literal MID-only selection therefore removes a recognized shallow-water protocol, changing habitat coverage at all17 comparisons. It is not a neutral cleanup.

[BOEM integrated-data definitions](https://www.boem.gov/sites/default/files/documents/regions/pacific-ocs-region/environmental-science/BOEM-2019-063.pdf); [Rassweiler et al. 2020 methods](https://www.frontiersin.org/journals/marine-science/articles/10.3389/fmars.2020.00272/full).

The complete-level sensitivity accepts BOT plus MID or CNMD. Relative to date-normalized standard surveys on common years, it changes recovery estimates at five comparisons: Abalone Cove, Gull Island, Point Vicente, Scorpion and South Point. It also affects Point Vicente before/after change. A missing level in this species-record delivery is a missing record of that pass, not proof that no diver surveyed it. Preserve that distinction until survey-completion metadata establish the reason.

## Coverage and method sensitivities

A restriction to SBTL_FISH_PISCO/VRG removes2004 at Anacapa2003, Santa Barbara Island and South Point because those records use SBTL_FISH_CRANE. This audit does not establish that CRANE was biologically invalid. Requiring complete BOT+upper records also removes Santa Barbara Island2005 under the fixed-role station rule. South Point2005 biomass remains separately unresolved because lengths are missing. Across34 histories there are420 original history-year rows; repeat normalization retains419 valid response rows, standard-method restriction413, and complete-level restriction411. These counts include before/designation years; they are not counts of independent surveys. Density remains available where biomass alone is missing.

Before/after arms of every numerical sensitivity share exactly the same available years and original station-role membership. We do not silently average over fewer stations when a station lacks a valid response. Separate complete-level and strict-MID contrasts distinguish incomplete-pass handling from removing INNER-zone protocols. Effects of losing years are listed, not conflated with the matched-year treatment effects.

## Recommendation

1. Correct repeat effort by averaging within transect number, zone and calendar year before existing hierarchical means. Keep all actual years, source campaigns and dates.
2. Preserve BOT+CNMD as a recognized protocol category alongside BOT+MID. Do not add standalone CAN observations or claim these two upper-water categories have identical habitat/detection. Keep protocol-aware and bottom-only sensitivities.
3. Flag single-pass units and seek survey-completion evidence before calling them true zero upper/bottom counts. Use complete BOT+(MID or CNMD) as a transparent sensitivity until that meaning is established.
4. Resolve what SBTL_FISH_CRANE denotes before automatically excluding its2004 surveys. Method labels alone do not establish incompatible sampling.
5. Keep researcher exclusions, same-year protected/reference support, response-specific missing measurements and existing offsets fixed during these checks. Freeze final depth/season rules separately; this audit has not applied the proposed shallow-water/date-gap policy.

## Outputs and checks

comparison_summary.csv identifies every affected comparison. affected_site_years.csv and campaign_date_ledger.csv retain exact station/year/source flags. effect_changes.csv and pooled.csv provide all six separate contrasts; coverage_changes.csv lists lost years. Source rows are never deduplicated by equal values. Models are diagnostic before/after linear contrasts and observed-year recovery trends with calendar-gap HAC3 and clustered pooled uncertainty. Four pooled biomass-change fits required a recorded BFGS fallback; those pooled confidence intervals remain unavailable from insufficient effective degrees of freedom. No new thresholds were chosen based on effect direction. Downstream pathway/resilience models were not rerun. Reproduce from the analysis root with code/tools/audit_pisco_fish_effort.py, code/tools/refit_pisco_fish_effort.R (Rscript), code/tools/report_pisco_fish_effort.py and code/tools/verify_pisco_fish_effort.py, in that order.
