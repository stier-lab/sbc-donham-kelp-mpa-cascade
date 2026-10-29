# Fish coverage recovery

## Why pairs were lost

The strict full-roster rule required every site ever retained at an MPA to occur
in every usable year. Among missing site/year/MPA links outside designation years,
447 have no recorded fish survey, 12 lack a standard-method survey and 4 lack a
matched BOT+MID visit/zone. These are missing links, not unique surveys or MPAs;
shared references may appear more than once. No complete site-year was lost
solely to the count/size quality filter in this audit.

Examples: CARP_REEF and LEAD_BETTER_BEACH occur only in 2009 in the Campus Point
matched sampling frame. Several Scorpion and South Point subsites ended in 2013;
South Point's BEE_ROCK and JOLLA_VIEJA subsites occur only in 2007-2009. Requiring
all these sites simultaneously discards long records at continuing sites.

No site-name merger was justified. In particular SCI_SCORPION_E and
SCI_SCORPION_W share site_new=SCI_SCORPION but have opposite protection status.
They must remain distinct. The CRANE method and CNMD strata were not relabelled
as standard PISCO/MID surveys.

## Recovery without changing membership through time

The main recovery diagnostic selects sites sampled in >=80% of the relevant
primary-history years, freezes that core, and then requires all selected sites
in every retained year. Before/after selection requires >=80% coverage in each
period separately. These are two endpoint-specific rosters, not one universal
MPA roster. No effects or P values enter selection.

Post-designation trend coverage increases from five to ten MPAs. Recovered MPAs:
Anacapa Island SMR 2003, Campus Point SMCA, Santa Barbara Island SMR, Scorpion
SMR and South Point SMR. The 70% and 90% thresholds also retain ten trend MPAs.
Before/after coverage increases from Naples alone to Naples, Campus Point and
Point Vicente; 90% retains two. Bottom-only recovery does not increase these
main-core cohort counts. Extending the window to all raw years does not increase
trend coverage and changes before/after eligibility because its coverage
denominator includes earlier sparsely sampled years. It is not an automatic
improvement over the primary-year diagnostic.

## What changes in the evidence

Eight of ten targeted-minus-warm-nontargeted trend estimates are positive.
Anacapa and South Point estimates are negative. The equal-MPA mean is 0.0287
log units/year (approximate shared-site-group interval -0.0043 to 0.0617;
ten MPAs, nine connected groups). The original five estimates are exactly
unchanged on common years; broader coverage, not recalculation of those pairs,
weakens the apparent consistency of the response.

Before/after differences are 0.279 at Naples, -0.117 at Campus Point and -0.621
at Point Vicente. Their lag-2 HAC 95% intervals are respectively [-0.296, 0.854],
[-0.348, 0.114] and [-1.125, -0.118]. Campus Point and Naples share reference
sites, leaving only two independent identity groups for these three comparisons;
no pooled confidence interval is reported. Point Vicente's before/after core
has two sites, whereas its post-trend core has eleven, so these are not conflicting
answers to an identical question.

This is defensible recovery of comparisons among regularly sampled sites, not
restoration of the complete MPA-wide sampling design. Site retention may be
related to reef conditions, and sampled zones can still vary within sites.
The broader results offer weaker evidence against siting than the five-pair
benchmark. Neither diagnostic establishes a protection effect or a cascade.

## Reproduction

Run `Rscript code/R/48_fish_coverage_recovery.R`, then
`Rscript code/tests/test_fish_coverage_recovery.R`. Module 48 follows module 47
in `run_all.R`. The standalone module and tests were run; the complete primary
pipeline was not rerun. Outputs use the `table_s_fish_recovery_` prefix. They
include the missing-link audit, site histories, identity conflicts, selection
scores, fixed rosters, exact years, densities, group contrasts, fits, pooled
results, same-year benchmark comparison and input/protected-file hashes.

All module 47 tables and primary harmonized inputs remain unchanged.
