# Malakhoff and Miller (2021) reconstruction

This is a standalone, **incomplete** reconstruction workspace. See
`docs/MALAKHOFF_2021_REPLICATION.md` for the source audit, specification
differences and missing-data requirements. All model outputs are quarantined
from manuscript use, even when fitting succeeds.

The best-available extension in `common_support.R` is called by `run.R`.
It freezes an original-window before-after-supported cohort across species
and compares 2019 with extensions through 2022 and partial 2023. See
`docs/MALAKHOFF_2021_AVAILABLE_DATA_RESULTS.md` for findings and limits.
The 75 additional fits include both period boundaries, independent-error
versus calendar-AR1 fits, reserve influence checks and before-period trends.
These are explicit modifications, not exact replication.

`method_bridge.R` now tests why the reconstruction and manuscript appear to
differ: matched observations for density/biomass, exact size accounting,
all-four-reserve density support, reference/aggregation algebra, primary
estimand reconstruction, calendar-AR1 and joint-synthesis sensitivities,
and strict-calendar geographic coverage. Read
`docs/MALAKHOFF_2021_METHOD_BRIDGE.md` for the findings and limits. It runs
automatically after `common_support.R`; `test.R` includes `test_bridge.R`.
The bridge also requires metafor. Outputs remain separate from the main paper.

Run from the analysis repository root:

```sh
Rscript code/replications/malakhoff2021/run.R
Rscript code/replications/malakhoff2021/test.R
```

Dependencies: data.table, nlme, ggplot2. Set `MALAKHOFF_DATA_ROOT` to override
the adjacent data repository's input directory. No raw input or current
manuscript-pipeline file is modified. Outputs are written to
`outputs/malakhoff2021/`, including an input-hash manifest, annual membership,
model results, serialized models and diagnostics.

The pinned author code is at
https://github.com/kmalakhoff/urchin-baci/tree/f892f7ace06e680c60175b1428ee8fc2463f4c71.
It does not include the required analysis CSVs. Latest-source discovery found
a 1982-2025 NPS package whose files are restricted; the available local files
do not provide a complete latest-data or original-period replication.
