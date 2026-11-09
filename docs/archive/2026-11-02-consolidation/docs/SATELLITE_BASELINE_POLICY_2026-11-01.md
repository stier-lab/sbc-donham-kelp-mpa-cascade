# Satellite baseline eligibility

**Later November 1 update:** the p-value exclusion policy below is historical. Adrian authorized an outcome-independent20-candidate audit and flat-baseline primary rerun with trend sensitivity. The frozen18-comparison revised satellite analysis is complete in `outputs/satellite_primary_20261101/README.md`; `CRITERIA.md` and its hash manifest control that analysis. Canonical143-history/manuscript-facing outputs still retain the earlier policy pending coordinated integration. Do not apply the old screen to the revised satellite cohort or silently insert its new mean-change estimand into the old pooled synthesis.

Adopted November 1, 2026: apply the same baseline rule to every satellite candidate while retaining the existing ecological, regulatory and whole-MPA exclusions. This decision concerns satellite canopy histories only. Naples LTER sheephead biomass remains included.

Fit annual protected–reference log response ratio against calendar year within each supplied Before period. Require at least three unique finite baseline years and a finite slope test with p >= 0.05. Histories failing this rule do not enter the primary synthesis as either pBACIPS or CI. Keep their observed data available for descriptive and explicitly separate analyses. Passing the test does not establish parallel trajectories or causal attribution.

| Satellite history | Baseline years | Calendar slope | p | Primary decision |
|---|---:|---:|---:|---|
| Naples | 28 | 0.0509313 | 0.057406 | Retain |
| Santa Barbara Island | 19 | -0.225798 | 0.0001209 | Exclude |
| Swamis | 28 | -0.0713588 | 0.040344 | Exclude |

Naples is close to the threshold; absence of a detected trend should not be described as proof of a flat baseline. Across all 20 satellite candidates only Santa Barbara and Swamis fail. Existing whole-MPA exclusions remain separate. The primary cohort becomes 143 histories (103 CI, 40 pBACIPS); diver membership is unchanged. Model families and corrected observed-calendar endpoints are unchanged by this eligibility decision.

Implementation: `code/R/satellite_baseline_policy.R`, called before satellite fitting in `08_effect_sizes.R`. Per-candidate diagnostics: `tables/table_s_satellite_baseline_screen.csv`. Test: `Rscript code/tests/test_satellite_baseline_policy.R`. Independent scipy linear regressions reproduced the three p-values to floating-point precision.

The pre-change tables, figures and scientific inputs are preserved in ignored `local_archive/satellite_screen_20261101/`. See the manuscript repository's integration record for downstream completion status; this policy document alone does not certify all document outputs.

## Primary rerun

The 143 retained individual effects are unchanged; only the two failing satellite histories were removed. Pooled estimates change because the joint models reestimate variation among MPAs and programs.

| Taxon | Response | Estimate | Adjusted p |
|---|---|---:|---:|
| S. purpuratus | Biomass | -0.5842 | 0.05222 |
| M. franciscanus | Biomass | 0.5710 | 0.05144 |
| M. pyrifera | Biomass | 1.1379 | 0.00043 |
| P. interruptus | Biomass | 0.4163 | 0.24246 |
| S. pulcher | Biomass | 0.9298 | 0.00295 |
| S. purpuratus | Density | -1.5595 | 0.00453 |
| M. franciscanus | Density | -0.4648 | 0.36010 |
| P. interruptus | Density | 1.0485 | 0.05144 |
| S. pulcher | Density | -0.2481 | 0.59419 |

## Verification and scope

Core modules 08–23 and downstream modules 24–59 (as listed in run_all.R), including the reference registry, pooled and temporal uncertainty, before/after coverage, programme and size comparisons, fitted-history displays and results summary, were rerun. Execution resumed after updating obsolete 143-versus145 cohort checks and the changed environmental-moderator minimum p-value (0.085); no environmental test survives its multiple-test adjustment. The separate 1000-draw size sensitivity and 250-draw second seed per weighting scheme were also rerun. The two removed histories remain in the harmonized source files; event-based satellite disturbance comparisons retain their separate window-coverage criteria. No raw data or source coefficients changed.

Passed: satellite selection and invalid-year tests; independent scipy trend tests; unchanged 143 individual effect comparison; exact-effect/HAC reconstruction; before/after grouping and shared-reference checks; selection-source/coverage fingerprints; calendar/routing OLS checks; size-weight/domain reconstruction and second-seed Monte Carlo check. Manuscript tables and current figures are checked separately in the private manuscript repository.


## November 1 — 16-pair baseline comparison, awaiting user review

Review-only follow-up: `outputs/satellite_baseline_comparison_20261101/README.md` , with eligibility table, exact annual inputs, full/10/5-year baseline comparisons, calendar-HAC lags1/3/5, baseline holdout predictions and seven figures. Both diagnostic models include a level change plus post slope; changing model adds pretrend. Estimand is average After departure from projected Before trajectory, not the canonical primary contrast. 569 annual observations across16pairs; no slope-sign/significance eligibility filter. Reference ecology remains qualified; Naples/Campus overlap documented.

Full baseline: SBI departure1.97→6.19, Swamis−0.53→0.81 (both uncertain), Naples−0.07→−1.02. Linear baseline wins last3baseline-year prediction in6/16 only; baseline-window sensitivity substantial. Do not automatically replace flat baseline or interpret long linear extrapolations causally. Independent NumPy/HAC reconstruction passes1152contrasts (maxdifference1.5e-11);487protectedfilesunchanged. Figures visually checked and margins repaired. Manuscript/supplement and primary143historycohort unchanged; Kumagai/data source roots unchanged. No send/push/commit, no running jobs. User explicitly requested review of figures/results BEFORE manuscript changes.


## November 1 — common-estimand pooled review refit

AR outputs/pooled_common_estimand_20261101/README.md is authoritative. Existing129diver+audited18satellite=147candidatehistories; no historical omitted diver histories restored. Frozen plan separates BeforeAfterchange (45histories;27diver18satellite) from postdifferences(145; FarnsworthPISCOsheepheadBio/Den only2postyears excluded in9yearwindow). Primarystrictprior baseline and observedyears1–9, designationyearomitted allprograms; originaltiming/fullAfter/baseline10/HAC1&5/trend sensitivities. SixfixedstationPISCOsheepheadhistories retain baselineeligibility; others notautomaticallypromoted.

Per-response REML with MPA/repeatedseries heterogeneity, sharedreferenceMPA CR2clusters anddf gates (no inference<2df or<3groups). Commonkelpchange+0.309 CI−0.159,0.778; postdifference+0.326 CI−0.141,0.792. FullAfterchange+0.608 CI0.142,1.074 butBHq0.128. Postcomparisons supporthigherredurchinbio,sheepheadbio/density(q~.025); no clearpurpledecline. Severalbaseline responses too sparse, not evidenceofnoeffect. OldnineTable2 reproduced; same19kelphistories oldmixed+1.146 versusnewcommonchange+.470; expanded23kelphistories+.309. Exactmatchedyears4MPAs stillNaplesdiver+1.382sat−.537. Do not equate measurementfootprints.

1333individualfits independentlyverified(<7e−13),238pooledrowsGLSreconstructed(<6e−15),CR2 arithmetic/BH, no duplicate scenariohistory keys,487protectedhashes unchanged. Implementation bugs caught andfixed beforeacceptedresults; finalrun clean. Comparativefigure inspected. MRscientificdocs/canonicalARresults/sourceworkspaces/Kumagai unchanged; no send/push/commit/jobs. Next userreviews design/estimand/window; historicaldiveromissions andfullsize/pathway/resiliencerebuild remain beforepromotion. Do not label newpostmean ascausalprotectioneffect.
