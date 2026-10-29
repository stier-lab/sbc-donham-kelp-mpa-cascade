# Monthly Kelp Calibration Refit

Implemented 12 October 2026 following Adrian's approval of Dan Reed's suggestion.

## Calibration

The pinned source is SBC knb-lter-sbc.112.15, Table 5, model 1, DOI 10.6073/pasta/710844a10e2fe8873a4b49d97999d549. Its twelve zero-intercept slopes estimate dry foliar standing crop, excluding holdfasts and sporophylls. Provider CSV bytes and EML are preserved under data/external/sbc11215/.

PISCO and LTER convert recorded-month observations before temporal averaging. KFM converts each dated stipe-size record, then resamples whole records so counts and months remain paired. Its existing site-year abundance matching is retained; size and abundance surveys are not necessarily month-matched. An unweighted monthly average is not used.

KFM's one-element sample(count, ...) branch previously sampled integers from 1 to count. Sampling record indices corrects one unpaired Anacapa Island 2007 density observation (0.007435 to 0.01). No paired annual density lnRR changed. The correction is separate from the biomass calibration. Legacy RNG consumption is retained for singleton samples to preserve unrelated bootstrap draws.

## Independent Checks

- PASS: harmonized_raw_responses.csv: 6,280 identical identities; 4,984 non-kelp records byte-preserved; 532 kelp records changed.
- PASS: harmonized_response_ratios.csv: 3,341 identical identities; 2,599 non-kelp records byte-preserved; 354 kelp records changed.
- PASS: Satellite inputs and site metadata are byte-identical to the checkpoint.
- PASS: All 144 primary histories retain their membership, N, model selection, analysis design and designation date.
- PASS: Non-kelp and Landsat individual primary contrasts and SEs are unchanged.
- PASS: Explicit registry re-export retained all 857 reference rows and 144 primary keys; only source-definition fingerprint changed within reference rows.
- PASS: Independent nine-test BH calculation reproduces all adjusted P values; all nine direction/detection classifications are unchanged.
- PASS: Historical Landsat pair-level and summary results are unchanged.

## Results

The pooled kelp log contrast changed from 0.802004 to 0.806823; exponentiated values are 2.230 and 2.241. These combine distinct contrast types and are not one physical inside/reference biomass ratio.

| Window | Paired geometric ratio | 95% CI | Signed-rank P |
|---|---:|---|---:|
| resistance | 4.22 | 1.37 to 13.01 | 0.00977 |
| recovery | 5.3 | 1.39 to 20.26 | 0.0137 |
| recovery_recent | 5.02 | 1.59 to 15.88 | 0.0137 |

All three modern windows retain ten paired MPAs. The recent-period sign test remains P = 0.109; the result is not test-invariant. The historical comparison still does not detect an increment after designation. Prior reef differences, shared controls and observational attribution limits remain.

## Reproduction

```sh
python3 -B code/audit_kelp_calibration.py --strict
Rscript --vanilla code/tests/test_monthly_kelp_refit.R
python3 -B code/audit_monthly_kelp_refit.py --before /path/to/preserved/before
```

The before/after audit requires the local checkpoint and is not a normal pipeline dependency. Routine input validation uses the pinned source bytes and kelp_monthly_calibration.json. All pre-change dirty files were preserved in local_archive/monthly_kelp_2026-10-12/before/ in each repository.

## Completed Workflow Checks

The final full pipeline exited successfully on 12 October 2026; its log is
`logs/pipeline_20261012_094719.log`. Current-input predictor and disease-moderation
exports were also regenerated. Existing warnings about sparse histories,
alternate rank-deficient fits and exploratory AR1 convergence remain visible.

The manuscript repository's `supplement/readability/monthly_kelp_qa_2026-10-12.md`
records the independent tests, numerical concordance, dataset citation and
changed-page visual checks. The exact source is cited through Zotero item
`NZF246ZD` / source ID `Rassweiler_2026`. Main-text figures remain five.
Previously circulated and frozen submission exports were not replaced.
No commit, push or circulation was made for this request.
