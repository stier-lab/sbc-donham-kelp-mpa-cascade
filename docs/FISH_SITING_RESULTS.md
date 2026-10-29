# Matched fish siting diagnostic

Module 47 implements `FISH_SITING_PLAN.md` without changing the primary analyses.
It reads raw PISCO fish observations, the supplied species attributes and site
map, and the retained sheephead-density reference registry and analysis years.

## Interpretation

The fixed-site, matched bottom/midwater screen retains five post-designation
comparisons: Abalone Cove, Gull Island, Harris Point, Naples and Point Vicente.
All five have a positive targeted-minus-warm-nontargeted trend difference.
The equal-MPA mean is 0.057 log units/year (approximate shared-site-group 95%
interval 0.002 to 0.112). Removing any one group leaves a positive mean.
The signs also remain positive across size, vertical-frame, species-panel and
zero-offset sensitivities, but several pooled intervals include zero. The
conclusion must not depend on declaring the main interval significant.

Only Naples meets the before/after rule: three pre-designation years and ten
post-designation years. The direct targeted-minus-warm-nontargeted step is
0.279 log units (calendar-gap HAC lag-2 interval -0.296 to 0.854). No network
before/after inference is available. No retained pair has six pre-designation
years for detrended fish tracking.

These are standardized density diagnostics for a fixed panel of fishes, not
biomass refits or reproductions of the primary all-noncanopy sheephead index.
They weaken an explanation of uniform improvement for all fishes, but do not
exclude species-specific habitat effects, recruitment or indirect protection
effects. The post-designation windows differ among pairs. Even with fixed raw
site IDs, sampled zones can vary across years and are exported for inspection.
Zeros are reconstructed only on observed standard fish survey cells. Entirely
unrecorded transects cannot be recovered from this encounter file.

## Reproduction

```
Rscript code/R/47_fish_siting.R
Rscript code/tests/test_fish_siting.R
```

The module is called by `run_all.R` in an isolated R process. It explicitly
reports missing raw inputs when those linked files are unavailable. A standalone
module run and synthetic/integration tests were completed; the full primary
pipeline was not rerun for this diagnostic.

Outputs are `tables/table_s_fish_siting_*.csv`: species panel and thresholds,
site/year/zone coverage, invalid cells, fixed rosters, exclusions, paired
densities, species positivity, annual group contrasts, preperiod tracking,
HAC fits and predictions, pooled effects, leave-group-out estimates, hashes,
protected-file checks and completion status. `fish_siting_sessionInfo.txt`
records the runtime. Tests cover fixed rosters, eligibility, designation-year
omission, direct group differences, size-filter monotonicity, shared identities,
consistent species years and unchanged primary inputs/figures.

Bob's note generates Figures 7-9 directly from these tables: coverage, observed
group trajectories with fitted trends, and Naples species before/after contrasts.
The manuscript and supplement are unchanged. Nothing was sent.
