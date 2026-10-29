# Retained pBACIPS shape display (Reed comment 116)

## Scope

Reed's original comment 116 (2026-10-06, in the annotated DOCX delivered
2026-10-08) requests a figure showing the patterns behind the discussion
of 27 linear and three sigmoidal histories among 36 pBACIPS histories.
Main Figure 4 instead displays pooled, response-specific linear mixed-model
predictions. Live SI Figure S8 displays residual diagnostics, and Figure
S11 summarizes retained model counts and residual dependence. The paired
history atlas shows observations, not the retained pBACIPS curves.

Script 32 fills that display gap without choosing illustrative linear
histories by their outcomes, fitting new candidate families, changing
retained labels or changing primary estimates. The compact comparison
includes every retained nonlinear history; the paginated display includes
every retained pBACIPS history, including Step models. CI histories remain
in the membership ledger with an explicit non-display reason. They are
not suitable additional evidence for the before-after model-count claim.

Retained model labels are read from the current primary and selection
tables. This script does not independently establish that these labels
are current AICc winners. The primary effect-size code includes explicit
per-history routing. This display is not a new model-selection analysis.
Neither a sigmoidal functional form nor a linear functional form by
itself demonstrates an observed plateau or indefinite biological change.

## Execution

```sh
Rscript code/R/32_pbacips_shape_display.R --validate-only
# Only after the parent's full pipeline has finished:
Rscript code/R/32_pbacips_shape_display.R --reconstruct
# Subsequently, no fits:
Rscript code/R/32_pbacips_shape_display.R --check-cache
Rscript code/R/32_pbacips_shape_display.R
Rscript code/tests/test_pbacips_shape_display.R
```

Ordinary execution never reconstructs a missing or stale cache silently.
Only `--reconstruct` fits the retained model for each pBACIPS history and
a same-data linear comparator for each retained nonlinear history. It
uses the existing standalone nonlinear definitions, not the candidate
selection routine. Definitions are read as parsed R expressions without
sourcing the effect-size pipeline. A non-NLS substitute, nonconvergence,
nonfinite predictions or a reproduced effect differing from the primary
table by more than 1e-7 causes a recorded failure and stops figure/cache
replacement. No primary cohort member is silently dropped.

The script is isolated with `local()`, uses qualified package calls, and
restores options and RNG state. It is safe to source in a shared pipeline
environment. `run_all.R` and `run_review_concerns.R` explicitly set
reconstruction mode after primary outputs are finalized, then restore
the previous environment setting. `run_figures_only.R` sets render mode
and never reconstructs; its selective key is `pbacips_shape_display`,
independent of SI numbering or legacy figure stems. No raw/upstream files
are read.

The focused regression suite uses sparse annual histories in a temporary
fixture. It checks clock units, Step direction and separator, exact sample
membership, horizon extrapolation flags, corrupt/stale-cache rejection,
no-fitting cached render, and caller-state preservation. Each runner's
scoped mode wrapper is also evaluated with a stub source, including error
paths, to verify restoration of both set and previously unset modes.

## Clock and uncertainty

The primary `add_time_columns()` function uses an observation-index clock:
all Before observations are zero; After observations are numbered
0,...,n_after-1. Thus the first After observation also has model time zero.
Missing calendar years are not elapsed time in that clock. Script 32
reuses the exact definition rather than replacing it with calendar time.
The explicit script-08 exceptions are preserved: Landsat linear histories
and Campus Point LTER lobster biomass use the supplied `time` column.
The ledger distinguishes `supplied_time` from `observation_index`;
observations retain the supplied time, original index and actual fit time.
For index-clock fits, calendar years are plotted by interpolating the index
between observed After years; Before predictions stay at zero model time.
Supplied-calendar fits instead use `pmax(year - MPA_Start, 0)` directly,
including between observations. `clock_units` identifies
`calendar_years_since_designation` or `after_observation_index`. The audit and
observation exports retain both clocks and list missing years.

Step contrasts are After minus Before. Progressive contrasts are evaluated
at model times zero and the existing common horizon (currently 11), even
when that horizon extends beyond a particular observed clock. The audit
flags such cases. Curves themselves are drawn only within observed
calendar-year support. There are no confidence or prediction bands: the
lines are descriptive retained shapes, not new inference accounting for
serial dependence, shared references or model-selection uncertainty.

In the 2026-10-10 cohort, 19 pBACIPS histories use the index clock and
17 use supplied time. The progressive contrast horizon extends beyond
observed model time for four LTER index-clock fits (maximum 10) and ten
2012-designated Landsat fits (maximum 9). A Step model has a period
contrast, not a progressive horizon. For the three sigmoidal KFM cases,
index 11 corresponds to the twelfth After observation, in 2015 rather
than eleven calendar years after 2003 designation. Describe the progressive
estimand as "fitted change from model time 0 to 11, using each
source-specific input clock", not uniformly "year 0 to year 11".

Exact reconstruction is a numerical concordance check, not a validation of
asymptote estimates. In this cohort the Gull Island lobster amplitude is
at the standalone fitter's upper bound (`M = 50`), and its curve still
rises at the last observation. Harris Point red-urchin biomass decreases;
its reconstruction records a singular-convergence warning from the existing
fitting sequence. Scorpion kelp has a steep transition (`K = 40.0033`).
Parameter values and warnings remain visible in the audit. These retained
shapes must not be presented as three independently validated biological
plateaus. Script 32 does not change their primary fits or fitting rules.

## Files and schemas

- `plots/fig_s_pbacips_nonlinear_comparison.pdf` and `.png`: all retained
  nonlinear histories, observations plus retained fit and linear comparator.
- `plots/fig_s_pbacips_selected_histories.pdf`: paginated retained-fit
  display for all pBACIPS histories; four stacked panels per page.
- `plots/pbacips_shape_histories/page_NNN.pdf` and `.png`: individual pages.
- `data/cache/pbacips_shape_display.rds`: numerical plotting cache, input
  and code MD5 fingerprints, and per-history reconstruction audit.
- `tables/table_s_pbacips_shape_membership.csv`: all primary keys,
  `Model`, `AnalysisType`, `N`, `MPA_Start`, `display_scope`,
  `nonlinear_comparison`, `nonlinear_panel`, `page`, `panel`,
  `primary_input_md5`, `clock`, `clock_units`, `n_finite_annual_inputs`,
  `annual_N_matches_primary`. All 144 annual sample counts must match
  the primary table; only pBACIPS histories have reconstruction checks.
- `tables/table_s_pbacips_shape_observations.csv`: source-specific annual
  `year`, `BA`, `lnRR`, fitted `time.model`, `input_time`, `index_time`,
  `clock`, `clock_units`, `time.true`, designation
  year, retained label and missing-year string.
- `tables/table_s_pbacips_shape_predictions.csv`: key, retained `Model`,
  calendar `year`, native-clock `time.model`, `clock`, `clock_units`, `BA`,
  `curve`, `fitted_lnRR`.
  `curve` is `Retained fit` or `Linear comparator`; it is not a new selection.
- `tables/table_s_pbacips_shape_fit_audit.csv`: all displayed keys,
  `N`, `primary_N`, `primary_effect`, `reconstructed_effect`, `difference`,
  `horizon`, `horizon_within_observed_clock`, `contrast_definition`,
  `calendar_year_at_horizon`, `first_after_model_time`, `clock`, `clock_units`, coverage, coefficients,
  `fit_class`, `fallback`, `status`, `error`, `warnings`.
- `tables/table_s_pbacips_shape_counts.csv`: retained labels, `n_histories`,
  `n_pbacips`, `n_primary`, `primary_input_md5`.
- `tables/table_s_pbacips_shape_provenance.csv`: tracked input/code paths
  and their MD5 fingerprints. A changed input or fitting code invalidates
  the cache and requires explicit reconstruction after the main run.

The parent assigns the final SI number and Results/Discussion callouts.
No manuscript or SI source is changed by script 32.
