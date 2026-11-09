# Anacapa kelp: are KFM and PISCO comparable?

Review only, November 1, 2026. No manuscript, supplement, canonical input, previous review output, or pooled model was changed. Reproduce with `python3 code/tools/audit_anacapa_comparability.py` and `python3 code/tools/plot_anacapa_comparability.py` from the analysis root; MPA_DATA_REPO can override the sibling data repository.

## Finding

The two programs measure a comparable biological quantity, but their Anacapa observations are not interchangeable samples of the same habitat. The disagreement is concentrated in the reference surveys. Across the same six years (2007–2012), mean protected biomass is almost identical: KFM 18.02 and PISCO 17.81 g dry foliar mass/m². Reference means are 0.279 versus 170.899. These are averages of annual site means, not uncertainty-weighted pooled effects.

The opposite directions already occur in plant density (protected/reference means 0.04583/0.000556 plants/m² for KFM and 0.02708/0.16539 for PISCO). They persist for stipe density, before applying biomass coefficients. Thus calibration and bootstrap randomness cannot explain the sign disagreement. KFM references contain only two counted plants across all three reference sites over these six years. The resulting log ratios are sensitive in magnitude to zero corrections; the raw observations show the discrepancy without a log transformation.

## Sites and spatial coverage

- KFM protected: Keyhole (a-k-31). References: Admiral’s Reef (a-k-11), East Fish Camp (a-k-32), Lighthouse (a-k-34).
- PISCO protected: West Isle CEN/E/W. References used in this comparison: Lighthouse Reef CEN/E/W every year, with Admirals CEN/E/W in 2007 and 2009. East Fish Camp is not part of this selected PISCO comparison.
- Keyhole’s recorded coordinate is about 104 m from West Isle CEN. KFM Lighthouse is about 93 m from Lighthouse Reef E. Admiral’s Reef and PISCO Admirals W have essentially the same nominal coordinate. Shared coordinates or names do not establish that the permanent KFM line and PISCO transects coincide.
- Holding references to Lighthouse alone preserves the disagreement. Restricting PISCO to the closest nominal protected station and corresponding nearest Admirals/Lighthouse stations also preserves it. Reference turnover alone does not explain the pattern. These are diagnostic restrictions, not new eligibility rules.

## Methods and dates

KFM’s modern protocol counts kelp plants taller than 1 m along a fixed 100-m line in forty 5-m² quadrats. The integrated file stores paired sides as twenty 10-m² records: 200 m² per site per year here. Adult and subadult plants are separate categories; both exceed 1 m. A separate size survey supplies stipes per plant. It can extend farther from the line than the abundance quadrats. Biomass is estimated from plant density and the sampled plant sizes. The original source-code comments describing quadrat counts as stipes are misleading; the arithmetic correctly uses plants times stipes per plant. Summing the two modern categories does not mean summing stipes twice. [NPS monitoring handbook, printed pp. 11–12 and 26–28](https://irma.nps.gov/DataStore/DownloadFile/485444).

PISCO records plants taller than 1 m and their stipe counts directly on 30 × 2 m swaths distributed among depth zones and stations. Here the hierarchy averages visits, transects, zones, then station means. Both conversions use the adopted month-specific SBC127.5 coefficients. [PISCO methods](https://opc.dataone.org/metacat/d1/mn/v2/object/doi%3A10.25494%2FP6%2FMLPA_kelpforest_methods.4).

KFM sampled between May and September; PISCO mostly sampled July. There is no shared year in which both the protected counterpart and an admitted reference counterpart have overlapping survey months. A strict same-month protected/reference comparison therefore cannot be made. The July 2009 Lighthouse surveys nevertheless provide a useful local check: PISCO sampled July 13 and KFM July 28. PISCO found kelp in its six shallow transects (4.0–5.2 m) but none in its six outer transects (7.6–9.8 m). KFM also counted none. This is consistent with habitat/depth differences, although it does not isolate depth from all local spatial variation or the 15-day timing difference.

## Depth checks

The raw PISCO records span substantial depth differences. A separate wave-sensor study lists nominal KFM-associated depths near 8.84 m at Keyhole, 14.33 m at Admiral’s Reef, 10.97 m at East Fish Camp, and 7.62 m at Lighthouse. These are later sensor-site measurements, **not verified survey-year depths for the KFM quadrats**. They motivate an exploratory depth comparison, not an exact match. [Lafferty and colleagues, Table 5.1](https://www.tomwbell.net/uploads/5/6/9/7/56976837/lafferty_etal_2019.pdf).

All diagnostics below use 2007–2012 and hold the original source-specific biomass zero correction fixed. They are descriptive mean log ratios without new inferential tests.

| PISCO sampling restriction | Mean protected/reference log ratio |
|---|---:|
| All included zones | −3.191 |
| Inner zones only | −3.769 |
| Outer zones only | −0.803 |
| 5–15 m | −2.812 |
| 6–15 m | −1.099 |
| 7–15 m | +0.580 |
| 8–15 m | +1.775 |
| 7–12 m | −1.779 |

The sign depends on the habitat included. The 7–15 m result does **not** establish that the programs agree after proper matching: narrower restrictions change station composition and available transects, and other plausible limits retain a negative contrast. Depth is represented by the median recorded within each survey-event key; event-specific minimum and maximum values are preserved. No KFM record has been assigned an invented survey depth.

## Recommendation before integration

Retain both full comparisons as program-specific observations and preserve their dependence at the MPA level. Do not interpret them as independent replications of a single Anacapa protection effect. Do not select a program, depth band, or reference subset because it produces a desired direction. For the eventual synthesis, retain the full-coverage comparison and an explicit omission sensitivity; any future depth-standardized primary rule should be defined across all relevant sites before rerunning.

The strongest explanation supported now is different sampled habitat, especially the reference footprint. Dates and separate size sampling can contribute, but their independent contributions are not identifiable here. A map of the actual KFM line and PISCO transects, with survey-year depths and habitat, would be needed to establish a closer match. There is no reason from this audit alone to declare either program wrong or remove its records. No new numerical correction is justified by the opposite signs.

## Figures and verification

- `protected_reference.png/pdf`: blue circles KFM, orange triangles PISCO. Each point is an annual mean across the program’s included stations. Both panels use the same scale. No confidence intervals are implied. KFM uses the exact expectation of the empirical size bootstrap; this removes Monte Carlo noise (mean log ratio 1.99496 versus prior bootstrapped 1.99542).
- `lighthouse_2009.png/pdf`: orange triangles are the individual PISCO transects on July 13; the blue point is the KFM site mean on July 28. KFM’s displayed approximate depth is the later sensor-site value, not a measured 2009 survey depth. Both deeper PISCO transects and KFM sampled zero kelp.
- `diagnostic_contrasts.csv` includes plant/stipe/biomass checks, fixed-reference and nearest-station restrictions, and all tested depth limits. Raw-record coverage, dates, distances, size pools, and zero corrections are separate CSVs.
- Independent PISCO reconstruction reproduces the preceding selected annual biomass values to <3.2e−12. Observed zeros are retained; unsurveyed stations are not zero-filled. No kelp count/size was missing in the selected raw PISCO records. Source hashes and verification results are saved. Exact-date/habitat matching remains unverified, not assumed.
