# Fish coverage recovery audit

23 October 2026. The strict module 47 comparison remains unchanged. Inspection
of site/year footprints, not fish effect estimates, identifies one-year sites,
retired subsites and nonoverlapping site histories as causes of exclusions.

Module 48 will classify each missing raw-site/year as no recorded fish survey,
nonstandard/unusable method, missing bottom/midwater coverage, or removal in the
quality filter. It will audit site_new mappings without merging names: some
subsites share a parent name but have opposite protection status. No inferred
alias, CRANE-to-standard conversion or CNMD-to-MID conversion is allowed.

Recovery comparisons use the same nine-species panel and raw-site identities:

1. Primary history years versus all available standard raw survey years through
   2023 at the same declared sites. Additional years are a separate diagnostic,
   not changes to the manuscript cohort. Denominators for raw-year coverage are
   years with any standard fish survey at a mapped site, not years selected for
   fish abundance or response.
2. Existing matched BOT+MID sampling versus independently reconstructed BOT-only
   sampling, without requiring MID. Match all species and size versions to the
   same cells and omit unknown count/size cells as in module 47.
3. Full historical roster versus a stable core present in at least 80% of anchor
   years; 70% and 90% thresholds are sensitivities. For before/after comparisons,
   each site must meet the coverage threshold separately before and after.
   For post-designation trends, select on post-designation coverage only. Then
   freeze that roster and use only years when every selected site is observed.
   Never replace a missing site with another reef or a missing survey with zero.

Require >=1 inside and >=1 reference site, >=3 pre and >=5 post years for
before/after, or >=8 post years spanning >=8 calendar years for trends. Retain
module 47's fish counts, species weights, size cutoffs, zero-offset sensitivities
and calendar-gap HAC inference (lag 2). Pool equal-MPA effects with shared-site
groups; do not use species or site cross-products as independent MPA replicates.
The main recovery scenario is primary-history years, BOT+MID and 80% core.
Alternative thresholds, raw years and bottom-only estimates are explicitly
separate. Selection is based on survey coverage, never effect sign or P value.

Report both recovered cohorts and effects restricted to original five trend
pairs on common calendar years so cohort changes are not mistaken for a change
in the result at the same reefs. Core-roster results concern regularly sampled
sites, not all sites in each MPA. Changing zone coverage and nonrandom survey
retention remain limitations. Added preperiod records do not certify reference
equivalence or parallel trajectories. Keep old estimates and outputs intact.
