# Frozen satellite eligibility and analysis specification

Frozen November 1, 2026, before the expanded-cohort refit. This is a transparent revision after earlier results were examined, NOT a preregistered or outcome-blind design. Site decisions below use treatment identity, source verification and coverage, never response direction, trend significance or influence.

## Scope and eligibility

1. Candidate universe: all 20 MPA labels in the verified 1984–2021 harmonized satellite delivery. This is not a census of all California MPAs or all 24 extraction-workbook rows.
2. Protection treatment: a named SMR/SMCA and documented designation year, with ordinary recreational AND commercial take prohibited for at least one focal predator (sheephead or lobster). Both no-take and qualifying limited-take areas are included; classify them separately. Tribal/research/incidental exceptions do not imply ordinary fisheries are open. This scope represents predator protection, not all possible effects of any regulatory restriction. No observed predator abundance threshold is used. Kelp harvest prohibition is not required (Naples exception).
3. Source validity: a distinct supplied MPA mask and nonempty reference mask(s), annual arrays matching the corrected delivery, consistent measurement/aggregation and unique finite pair-year ratios. Existing equal-available-reference averaging is retained. Known placeholders are not observations; actual biological zeros retain the established pipeline handling. Do not modify source data or fill missing years.
4. Reference standard applied to every pair: supplied geometric assignment with verified data correspondence and no demonstrated invalidity. This is NOT independent verification of habitat equivalence, legal-boundary geometry or untreated reference status at every pixel. Document those limitations for all pairs, not only newly restored ones. Overlap with an omitted extraction is not grounds for omission; overlap among included comparisons is dependence to disclose and assess, not independent replication.
5. Temporal coverage: at least three unique paired baseline years, three paired After years, and at least ten years between first/last observations. Use all available valid years within the delivered period. Preserve audited satellite timing: designation year belongs to Before; After starts the following year. Do not borrow PISCO's different annual timing.
6. Baseline slopes, p values, effect magnitude, direction and model fit do not determine inclusion. No Cook's-distance or result-driven removals. Full versus shortened baseline and uncertainty lags are sensitivities, not selection devices.

## Four audited decisions

- Point Conception SMR: INCLUDE. Named no-take reserve, 2012 designation, verified 37-year paired record. The old broad data-quality omission has no demonstrated satellite measurement defect; Tom confirmed its second reference and corrected masks reproduce annual arrays. Low/rare focal predators, if raised as biological scope, are not a satellite data-validity criterion.
- Judith Rk SMR: INCLUDE. Named no-take reserve established 2003; supplied extraction identified separately from San Miguel SC; verified 37-year paired record. Its 364 shared MPA pixels with San Miguel do not duplicate an included comparison because San Miguel SC is excluded. The overlap alone does not invalidate Judith.
- San Miguel Island SC: EXCLUDE from primary. A boating/access special closure is not a verified predator-fishing protection treatment; the harmonized row is misleadingly typed SMR and dated 2003. Do not treat that coding as legal confirmation or silently relabel it. Its valid kelp measurements remain available.
- Arrow Point to Lion Head Point SMCA: EXCLUDE from this predator-protection primary. Recreational finfish and all commercial take remain allowed, so neither focal predator has the complete ordinary recreational/commercial protection required here. Recreational invertebrate take IS prohibited: do not describe it as unprotected. This is a scope decision, not bad data. Its coded 2003 year also cannot be assumed to represent the current SMCA treatment; no alternative-effect fit is conducted using an unverified clock.

The same protection rule keeps Swamis (lobster protected despite shore fishing), Naples (focal predators protected despite kelp harvest), and the retained pelagic-exception SMCAs. All 20 receive explicit audit rows.

## Frozen models, estimand and uncertainty

Primary satellite model: annual log(MPA/reference biomass) ~ post + post_time. This has a flat Before mean and separate After intercept/slope. Primary effect is average predicted After minus the flat Before mean, algebraically equal to observed After mean minus observed Before mean on the same rows.

Trend sensitivity: add centered calendar year, retaining identical rows and allowing Before/After slopes. Effect is average After departure from projection of the fitted Before trend. A constant extrapolated trend is an assumption, not evidence of causality.

Full-baseline and calendar-gap HAC lag3 uncertainty are primary. Retain the previously specified 10-/5-calendar-year baseline windows, HAC lags1/5 and last-three-baseline-year predictive holdout as sensitivities. Use n/(n-p) covariance correction, residual-df t intervals and BH adjustment within each model/window/lag/quantity family (18 comparisons in the expanded cohort). No selecting a model from its p value or AICc.

Report site estimates and descriptive equal-site mean/median; do not treat shared-reference sites as independent or attach an independence-based pooled CI. Compare the frozen18 cohort with old14/16 using the SAME models to separate model changes from inclusion changes. The paper's current pooled synthesis mixes other estimands; do not silently replace its satellite effects with these mean shifts. Diver fits and all canonical manuscript-facing outputs remain unchanged pending review.

## Evidence

- https://wildlife.ca.gov/Conservation/Marine/MPAs/Judith-Rock-San-Miguel-Island (checked November1)
- https://nrm.dfg.ca.gov/FileHandler.ashx?DocumentID=100879 (Judith2003 establishment)
- https://wildlife.ca.gov/Conservation/Marine/MPAs/Point-Conception (checked November1)
- https://wildlife.ca.gov/Conservation/Marine/MPAs/Arrow-Point-to-Lion-Head-Point (checked November1)
- docs/LANDSAT_V2_VERIFICATION_2026-10-21.md; data/harmonized/landsat_v2_verification.json
- tables/table_s_selection_judith_pixel_overlap.csv; docs/SELECTION_CRITERIA_AUDIT_2026-10-29.md
- Original archived pBACIPS_PISCO_V10.R lines6690–6694; private manuscript selection_history/HISTORY.md and email_evidence.json (Tom October20 clarification).

Missing historical ecological matching details do not become confirmed by this audit. Exact legal-polygon overlays and a common estimand for pooled diver/satellite inference remain separate tasks, not reasons to select histories using their results.
