# Required source years, Pelican Bay, area and missing data

November 2, 2026. Review only; paper and frozen results unchanged.

## Required-year check: fixed

The four held KFM groups from 2023 remain unavailable but no longer block the audited current matched histories, which use none of them. The source build now reports one required hold (Pelican Bay red urchins, 1995) and four deferred holds. The explicit group policy is hash-checked. A request to validate biomass support for a deferred group fails; broad preparation can retain its count with biomass marked missing, rather than halt unrelated histories or borrow around the source gap. A later extension needs a reviewed scope change and verified sources. This is not permission to add 2023 observations.

## Pelican: five original measurements remain unexplained

All 140 original red-urchin rows are maximum test diameters from Pelican Bay on August 7, 1995. The combined file contains 135 animals and is an exact frequency subset of the original. It omits one measurement each at 51, 63, 105, 134 and 184 mm. This is not a general unit conversion or a single upper-size cutoff. Equal sizes are legitimate repeated observations; the data do not identify individual animals, so arbitrary deduplication would be unjustified.

No recorded exclusion reason was recovered for these five measurements. Recommendation: prefer the original 140-row source under the source-first policy, with the 135-animal version as a sensitivity. The required source hold is retained pending that source choice; this audit does not claim the original version has been independently adjudicated.

The existing same-data/same-model comparison changes the Scorpion red-urchin before–after estimate from 1.253 to 1.270 log units. Conditional 95% intervals are 0.704–1.803 and 0.719–1.822, respectively. Both remain positive; post-protection trends are unchanged because the affected observation is before protection. These are previously verified conditional fits, not a new joint-uncertainty analysis. The affected 1995 reference-pool mean mass changes by 39.4%, so the small pooled-effect change should not be mistaken for a negligible input discrepancy.

## LTER area and missing-size handling: checked and clarified

- All 23,236 source LTER urchin observations have an area of 1 m² and a recorded nonnegative count. The suspected missing area conversion therefore causes no numerical difference in this dataset. The production calculation now explicitly divides total biomass by quadrat effort and area; invalid effort or missing counts cause an error.
- The existing 1,156-row bootstrap cache contains 72 missing biomass estimates. They affect 30 physical-site/year/species combinations. Previously the aggregation temporarily filled those values with zero before excluding unsupported biomass. The revised branch preserves missing biomass directly and produces the same retained numerical rows.
- Explicit area normalization changes the checked cached values by exactly zero. This is a source/aggregation comparison, not a fresh bootstrap or model refit.
- No repeated dates occur within LTER physical-transect/species/year groups, and no physical stations belong to both of these LTER comparison rosters. This resolves those particular concerns for this branch, not for every program.

Density remains supported when sizes are unavailable; the source count is not changed to zero or removed because biomass is missing. Group-average size limitations and complete uncertainty remain separate issues.

## Broader source check: further missing-data cases

The inventory preserves program, survey method and actual area. KFM areas vary by method/year and must not be replaced by one blanket constant. No missing or nonpositive areas were found among the focal MBON rows checked.

PISCO source records needing targeted handling are:

| Site | Year | Recorded value missing | Why it matters |
|---|---:|---|---|
| Cavern Point West | 2013 | Red-urchin count | Scorpion protected station; unknown count must not become absence |
| Little Dume West | 2011 | Kelp size/stipe count, with two plants recorded | Point Dume; plant count does not determine stipes |
| Underwater Arch | 2023 | Kelp size/stipe count, with two plants recorded | Point Vicente reference; check use in final matched coverage |
| Point Loma Central | 2023 | Kelp size/stipe count, with 18 plants recorded | South La Jolla is outside the retained ecological scope |

These are documented, not numerically corrected here. The next targeted check should trace each retained record through PISCO's completed survey grid and zero-filling step, retaining explicit unknown measurements separately from inferred species absences. VRG lobster/urchin swath sizes being blank is expected: their size measurements are in the separate accepted frequency table.

Five MBON LTER kelp counts are also missing (Naples 2000 and four Santa Cruz Diablo 2004 rows). They are listed separately. The primary LTER kelp calculation uses the separate annual kelp delivery, so these are not automatically five errors in that primary response.

## Tables and reproduction

- [Area and missingness verification](area_missing_verification.json)
- [LTER source areas](lter_source_area.csv)
- [Missing biomass rows](lter_missing_biomass.csv)
- [Explicit versus legacy area calculation](lter_area_comparison.csv)
- [Pelican size frequencies](pelican_size_frequencies.csv)
- [Pelican source inspection](pelican_review.json)
- [Cross-program source-area inventory](mbon_area_missing_inventory.csv)
- [PISCO records with missing focal values](pisco_missing_focal_values.csv)

Reproduce the area/missingness comparison with `code/tools/audit_area_missingness.py`; its cached-input CSV is an export of the existing LTER RDS, not a new fit. The data repository `make check` tests source gating, missing values and effort/area conversion. No manuscript, supplement or frozen primary output was changed.
