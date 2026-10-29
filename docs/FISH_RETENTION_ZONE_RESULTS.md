# Retention and within-site sampling results

23 October 2026. Optional module 49 implements the recorded plan in
`FISH_RETENTION_ZONE_PLAN.md`. Earlier results and primary inputs are unchanged.

## Retained versus excluded sites

Within the ten eligible core80 trend pairs, five MPAs have excluded sites:
Anacapa, Campus Point, Santa Barbara Island, Scorpion and South Point. All 34
excluded-site links have overlapping years with a fixed retained comparator
roster of the same protection status. This is a descriptive check, not 34
independent experimental replicates.

In the first three available overlapping years (or fewer for short histories),
28/34 excluded sites have lower mean species log targeted-fish counts: 10/12
protected sites and 18/22 references. Warm-affinity non-targeted fishes are lower
at 22/34. Using all overlapping years, targeted fishes are lower at 30/34.
Bottom-depth differences vary in sign; retained and excluded reefs are not
interchangeable on the measured variables.

Only the two Campus Point reference comparisons contain a pre-designation
observation, each from one year. The comparisons therefore do not establish
untreated baseline habitat differences, the reason sampling stopped/started,
or a selection bias specifically favoring protected reefs. Canopy cover is
entirely missing in this fish frame, and substrate/exposure are not supplied.
Retention results for the other five screened MPAs and for the before/after
selection rule are also exported, but are not pooled into the ten-pair diagnostic.

## Fixed survey zones and bottom-depth ranges

All ten trend pairs retain zones sampled every original retained year at every
core site. Holding those zones and their weights fixed changes the mean trend
from 0.028668 [-0.004326, 0.061662] to 0.018920 [-0.018943, 0.056783] log units/year.
These comparisons use identical sites and years.

Requiring BOT transects within 2 m of the site-zone median bottom depth retains
eight trend pairs, with 0.024732 [-0.022580, 0.072044]. The original calculation
on exactly those eight pairs and retained years is 0.035602 [-0.008088, 0.079293].
Point Vicente and Santa Barbara Island lose trend eligibility. The 1 m and 3 m
checks retain three and nine trend pairs and give -0.003709 [-0.151263, 0.143845]
and 0.019163 [-0.020638, 0.058964], respectively. Cohorts differ across ranges.
Intervals are approximate equal-MPA delete-shared-site-group intervals, not
evidence that all protection effects are zero.

All three before/after pairs remain under fixed zones and the 2 m restriction.
Point Vicente is particularly sensitive: -0.621420 [-1.125333, -0.117507] in the
previous calculation becomes 0.069961 [-0.489866, 0.629789] with fixed zones.
The 2 m result is -0.112086 [-0.860575, 0.636403]. All fixed-zone
before/after intervals include zero. Shared sites leave only two independent
identity groups, so no pooled before/after interval is available.

Depth restrictions are within-site checks, not matching habitat between
protected and reference reefs. Almost all VRG MID depths are missing; only
bottom depths are constrained. Both vertical survey strata remain required,
missing bottom depths are not imputed, and actual depths can vary inside the
fixed ranges. Neither check resolves siting or establishes a cascade.

## Reproduction and checks

Run `Rscript code/R/49_fish_retention_zones.R`, then
`Rscript code/tests/test_fish_retention_zones.R`. Module 49 follows module 48
in `run_all.R`. Tables have prefix `table_s_fish_sampling_`. They include
retention values and comparator rosters, cell-depth metadata, fixed-zone
rosters, depth centers and restriction audits, retained years, observed annual
contrasts, fits, pooled estimates, and input/protected hashes. The unrestricted
reconstruction reproduces module 47 densities to numerical precision.

Bob's condensed R Markdown note reads these outputs for Sections 3-4 and
Figures 3-4. The fuller earlier note is preserved in its dated backups.
No main manuscript or supplement changes were made. The standalone diagnostic
and its tests were run, not the full primary pipeline.
