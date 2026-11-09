# Primary processing labels

The 145 retained histories comprise 103 control-impact (CI) and 42 progressive before-after (pBACIPS) contrasts. Of the CI histories, 85 have fewer than three paired baseline years and 18 have at least three. Six of the 89 PISCO histories use fixed-station before/after comparisons under the three-strictly-prior-year rule; the other 83 use CI. Fish and benthic coverage are assessed separately across retained responses within each survey group. A CI label does not by itself establish failure of a baseline-trend screen.

Campus Point LTER lobster biomass is a calendar-time linear CI trend, with one observation in the designation year and none strictly before designation. Its source helper exports CI, BA=N and LinearBefore=NA. The fitted-history atlas contains the 42 pBACIPS histories plus this explicitly labeled CI trend. Table S2 and the primary export use the same labels.

Run `python3 code/tools/check_processing_labels.py` after regenerating the primary and matched-coverage outputs. It verifies exact identities and routes. The separate matched before/after analysis has 60 qualifying histories and is distinct from primary model selection. Its coverage table records the actual before/after labels for each series; the six PISCO sheephead comparisons treat 2012 as after implementation.

See `docs/PISCO_BASELINE_POLICY_2026-10-31.md` for the current baseline rule, station selection and four-year sensitivity. Historical dated audits describe the versions they evaluated.
