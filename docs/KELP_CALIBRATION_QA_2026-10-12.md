# Kelp Biomass Calibration: Exact-Source Check

**Historical pre-adoption audit.** Adrian subsequently authorised the monthly
calibration. Current implementation and refit status are in
`MONTHLY_KELP_REFIT_2026-10-12.md`; the observations below describe the preserved
pre-change version, not the current pipeline.

Checked 12 October 2026. **Source verification is complete; adoption of the
monthly calibration is not.** The fixed production slope 0.08 does not
implement SBC `knb-lter-sbc.112.15` Table 5. No production constant, input
value, model, table, figure or manuscript was changed by this audit.

## Verified Source

The package was published on 9 September 2026:
[SBC 112.15](https://doi.org/10.6073/pasta/710844a10e2fe8873a4b49d97999d549).
The original EML and fifth CSV entity are retained in
`data/external/sbc11215/`, with their license and retrieval details. Direct
PASTA access returned 403; DataONE's exact-object resolve service recovered
both objects. The CSV has 108 rows and its MD5 matches the EML:
`50eda5d3664bf70780b3bead9674fac1`.

Table 5 means the fifth data entity, not the rows labelled model 5. **Model
1** converts frond density (>1 m tall, per square metre) to total dry foliar
standing crop in kg dry mass per square metre. Foliar biomass excludes
holdfasts and sporophylls. Other models estimate canopy or water-column
components, subsurface biomass, and production; they are not replacements
for total foliar standing biomass.

Model 1 has zero intercepts and these twelve monthly slopes:

| Month | kg dry mass per frond |
| --- | ---: |
| January | 0.08648 |
| February | 0.09158 |
| March | 0.10540 |
| April | 0.11654 |
| May | 0.11996 |
| June | 0.11355 |
| July | 0.10495 |
| August | 0.10138 |
| September | 0.09043 |
| October | 0.08687 |
| November | 0.07272 |
| December | 0.07965 |

The associated [SBC seasonal biomass protocol](https://sbclter.msi.ucsb.edu/external/Reef/Protocols/Seasonal_Benthic_Survey/Seasonal_Benthic_Survey_Protocol_Biomass_invert_algae_fish_calculation.pdf),
page 5, independently describes month-specific conversion from frond density
to dry biomass. The original calibration paper is Rassweiler et al. (2018),
*Ecology* 99:2132, https://doi.org/10.1002/ecy.2440. The exact numerical check
uses the 2026 package CSV, not an assumed unchanged 2018 supplement.

## Current Calculation And Discrepancy

The data-processing repository's `code/R/01_utils.R` defines
`MACRO_AVE_SLOPE <- 0.08` and `bio_macro(stipe) = stipe * 0.08 * 1000`.
PISCO applies it to stipe density; LTER applies it to frond density before
annual aggregation. KFM bootstraps annual site-level stipe counts per plant,
converts them with the same function, then combines mass with plant abundance.
None of these calls selects a slope by survey month.

Reed et al. (2009), Figure 2b, page 2081 (PDF page 5), was checked visually
and against the local full text. It reports `FSC = 0.08 * frond density +
0.01`, in **dry** kg per square metre. Our zero-intercept fixed conversion
uses its slope, not its full fitted equation. KFM's code comment describing
that conversion as wet mass is incorrect. The main Figure 5 axis currently
says grams per square metre without a mass basis; a future refresh should
specify dry foliar biomass if this response is retained.

Neither the old six-value set (0.10386, 0.10103, 0.09267, 0.09204, 0.08054,
0.08505) nor its mean is the verified Table 5 model-1 calibration. The
unweighted twelve-month mean is 0.09745917 and the May-October mean is
0.10285667. Neither is the survey-month lookup and neither was adopted here.
Changing 0.08 to either mean would remain an approximation.

## Dates And Consequences

The raw files retain months for all three diver programmes. A read-only
inventory found:

| Raw input scope | Kelp rows | Missing/invalid months | Site-years with multiple months |
| --- | ---: | ---: | ---: |
| LTER annual kelp file | 29,526 | 0 | 3 of 250 |
| KFM Macrocystis raw file | 65,215 | 0 | 5 of 612 |
| PISCO swath rows labelled MACPYRAD | 84,581 | 0 | 223 of 2,327 |

These counts describe complete raw-file scopes, not the retained primary
cohort. LTER observations occur in July-October; KFM in May-October; PISCO
kelp records include all twelve months. Raw-field parsing used Latin-1 to
retain non-UTF-8 notes; month and site identifiers used here are ASCII.

A common constant cancels from proportion-normalized lnRR. Month-specific
coefficients need not: months and their contribution can differ among years,
sites and statuses. Applying a slope after annual averaging would discard
that information. KFM additionally needs an explicit date-aware decision for
combining stipe-size and plant-abundance observations, rather than inventing
months for annual bootstrap summaries. Consequently this audit does **not**
claim unchanged annual ratios, primary contrasts, model selection, pooled
results or resilience estimates after a monthly-conversion change.

## Required Decision And Follow-Through

Choose either the dated monthly calibration or a documented fixed-slope
approximation. A monthly implementation should convert dated source
observations before the existing spatial and annual aggregations, preserve
non-kelp records and all inclusion decisions, and check KFM date alignment.
Then regenerate harmonized kelp inputs, refit the full analysis, inspect
membership and uncertainty changes, and refresh affected tables, figures,
main-text/SI wording, units and source citations together. Do not cite Table
5 as the source of the unchanged 0.08 calculation.

The harmonization repository remains archived on GitHub. This check did not
unarchive or push it, and did not silently promote a new calibration to the
public analysis inputs.

## Reproducible Checks

```bash
python3 -B code/audit_kelp_calibration.py
python3 -B code/audit_kelp_calibration.py --strict
python3 -B -m unittest discover -s code/tests -p 'test_*.py' -v
```

The first command reports source **PASS** and production **REVIEW**. The
strict command deliberately exits 2 while that decision is unresolved; it
does not certify the current calculation as the monthly calibration. Six
new synthetic tests plus the 33 existing audit tests pass. The source check
verifies the pinned metadata SHA256, package revision, fifth entity identity,
provider MD5, row count, response/predictor units and labels, all twelve
months, positive finite slopes and zero intercepts. It writes no files.

Fingerprints of eleven production anchors were checked before and after:
the harmonized raw and ratio files, primary effect sizes, Tables 1-3, Figure
5, the live manuscript, both legacy supplement source/Word files and the
separate readability Word supplement. All remain unchanged. No model was
refitted and no email, upload, commit or push was performed.
