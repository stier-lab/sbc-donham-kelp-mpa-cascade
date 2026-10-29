# Script 32 verification and visual QA, 2026-10-10

## Coverage

The corrected production primary table has MD5
`72582d436d7d070f36965cefb458d601`. The parent's completed pipeline log is
`logs/pipeline_20261010_121815.log` (exit 0; finished 12:22:26). It records
script 32 completing successfully through the full-run reconstruction
wrapper. Script 32 was also explicitly reconstructed and cache-checked
standalone after that run and after the parent's clock-comment updates.
The final cache MD5 is `361a36acc4c08619df2519b31c5ce7ac`. Its current
script-08 MD5 is `da9e31f11df5da453cd430ee1f1e4c5e` and constants MD5 is
`bee769d4e7d60a7c54e64d05a9f16ecd`. Final core comparison PDF MD5:
`6650250f7fb69a33429b078e7eeaa02e`; nine-page fit PDF MD5:
`3f6dcd3f6e9a3b7452aa7a16b6651df7`.

Individually inspected 17 current PNGs with `view_image` requesting original
resolution, not contact sheets: Atlas B plot pages 001-009 (all 36 panels),
the three-panel nonlinear comparison, all five refreshed main figures,
the design comparison, and the historical heatwave comparison. The viewer
resizes large images for display. This review does not cover Atlas A's 144
raw-history panels, the parent's Atlas B title/index packaging, or a new
rendered manuscript/Word file.

The first pass found long MPA labels crowding panel letters on Atlas B plot
pages 3 and 7. MPA labels now wrap without reducing font size; all nine
final plot pages and the compact comparison were inspected again. Panel
letters, points, fit labels, axes and units are legible and unclipped.
No overall titles or subtitles were added. Both new PDFs embed their
Helvetica regular/bold fonts; the retained-fit PDF contains nine pages.

## Checks

- All 36 retained pBACIPS fits verified: 27 Linear, 6 Step, 3 Sigmoid;
  no failures or dropped keys. Maximum absolute effect discrepancy is
  `1.110223e-15`, against the unchanged `1e-7` acceptance tolerance.
- The six Step effects independently equal the annual After-minus-Before
  means. The display's jump separator is designation year + 0.5, not an
  exact designation date.
- All 144 membership keys and finite annual sample counts match primary
  inputs; 108 CI histories have an explicit non-display reason.
- All three retained Sigmoid histories appear in core panels A-C. Their
  MPA/taxon/response/source membership matches the legend and ledger.
- Native clocks retained: 19 observation-index histories and 17 supplied
  calendar-time histories. Fourteen progressive contrast horizons extend
  beyond observed model time and are flagged. Step horizons are NA.
- All 16,794 prediction rows are finite, bounded to observed calendar
  support, and checked against each native clock and retained label.
- Focused synthetic regression tests pass, including sparse calendars,
  Step direction, cache corruption/staleness rejection, cache-only render
  with fitting functions disabled, shared-state preservation in reconstruction
  and render modes, and each
  runner's mode restoration after both success and error.

Scored checks (1 = defective, 5 = no detected defect): numerical concordance
5/5; cohort/panel traceability 5/5; clock/legend agreement 5/5; layout 5/5;
readability 4/5 (dense but readable atlas axis labels). These scores do not
grade the adequacy or causal interpretation of the underlying primary fits.

## Scientific limitations

No reconstruction or layout blocker remains. The retained nonlinear shapes
are not three demonstrated biological plateaus: Gull Island lobster is
still increasing at the end of the record and has `M = 50` at a fitting
bound; Harris Point red-urchin biomass decreases and records an optimizer
warning; Scorpion kelp has a steep fitted transition (`K = 40.0033`). The
legend explicitly distinguishes the displayed shapes from biological
saturation, new model selection, uncertainty intervals or causal inference.

The additional script-11 clock lookup confirmed that its labels at the
Figure 4/biomass-density legend belong to calendar-time temporal models,
not the S1 processing illustration. All 2,080 eligible After rows satisfy
`time = year - designation year`. The correct "t = 11 years" labels remain;
only misleading comments equating them with every primary pBACIPS route
were corrected. S1's simulated trend uses `lnRR ~ year`, with no 11-unit
label, so no processing-figure rerender was needed.

## Changed files

Owned code changes are `code/R/32_pbacips_shape_display.R`, the script-32
wrappers in `run_all.R`, `run_review_concerns.R`, and `run_figures_only.R`,
and the two misleading clock comments in `code/R/11_figures.R`.
The focused test is `code/tests/test_pbacips_shape_display.R`. Documentation
is `docs/pbacips_shape_display.md`, `docs/pbacips_shape_legend.md`, and this
verification note. Generated artifacts are the shape cache, six
`table_s_pbacips_shape_*.csv` files, the core PNG/PDF, nine individual
PNG/PDF plot pages, and the combined nine-page fit PDF. No primary analysis
rules, manuscript files or main figures were changed by this task.

## Inspected PNG fingerprints

Paths below are relative to the analysis repository. PDF file fingerprints
can change on rerender because Cairo includes creation metadata; PNG hashes
identify the visually checked content.

| File | MD5 |
| --- | --- |
| plots/fig_s_pbacips_nonlinear_comparison.png | f34cfe6e3245c00b6af2982853a671de |
| plots/pbacips_shape_histories/page_001.png | f7271642d0778c6c512df20bcb85e39b |
| plots/pbacips_shape_histories/page_002.png | 891f83e3254d915fe69d4ec21664e217 |
| plots/pbacips_shape_histories/page_003.png | afc2bf555b0747ef97f2868d280ebe09 |
| plots/pbacips_shape_histories/page_004.png | bd2de6821183e5026279978d53921744 |
| plots/pbacips_shape_histories/page_005.png | 7f5c86c36d6cd264542622db054a90bf |
| plots/pbacips_shape_histories/page_006.png | 265b0f54031b0d0b4b3ebfe95e545ef7 |
| plots/pbacips_shape_histories/page_007.png | aed22c5e9eecf47003a287872eb1efd9 |
| plots/pbacips_shape_histories/page_008.png | 570c9eb8435a5b9deafc5613ef991765 |
| plots/pbacips_shape_histories/page_009.png | ac776f21aa329b409e35db5e311b2a40 |
| plots/fig_01_mpa_map.png | 9db09e7e45efb947f8a9e39952073d91 |
| plots/fig_02_cascade_case_studies.png | b9160a7218485b818dabfd6414ebcf96 |
| plots/fig_03_mean_effects.png | faca97b343f88ec2397cdad8978400be |
| plots/fig_04_recovery_curves.png | 1da25b042c5730dbfc7d7d6ae8d09086 |
| plots/fig_kelp_resilience.png | cdeaee6a6e45136a164947ff7fdbb4d2 |
| plots/fig_s_design_matched.png | 4035755e91aae552b089b26a1dc1df3b |
| plots/fig_s_historical_heatwave_bacips.png | 78e6d2d2c1032be7bd2f2823c6a833f3 |
