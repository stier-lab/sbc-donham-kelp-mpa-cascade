# Remaining KFM source checks and analysis impact

## What is resolved

**Hare Rock, purple urchins, 1987:** the original legacy table contains102measurements. The NPS1982–1989report, PDFpage274/printedAppx.5-26, independently givesN102, matching all20published size bins and the3–56mm range. The combined file has103animals and five differing size-frequency entries. Use the original individual measurements; the published bins corroborate them but are not converted into invented individual sizes. The exact reason for the later discrepancy remains undocumented. [NPS report](https://irma.nps.gov/DataStore/DownloadFile/640174).

The verified identity bundle now replaces this single site/year/species group with its 102 original records. The original raw files are preserved. The guard verifies unique source IDs and retains all researcher exclusions. The replacement reduces estimated mean animal mass for this protected-site sample by 2.19%. It affects a before-protection point in the Harris Point purple-urchin review history, which is not a primary effect in the current manuscript. Every post-protection recovery observation and trend is unchanged.

## What remains uncertain

**Pelican Bay, red urchins, 1995:** the legacy original contains 140 measurements; the combined file contains 135. The latter lacks one animal at each of 51,63,105,134 and 184 mm relative to the legacy distribution. The larger animals explain why the legacy mean mass is 39.39%higher even though the sample totals differ by only five. Do not silently discard the large sizes or choose the version by which produces the preferred effect.

The downloaded1995NPSannual report describes database conversion and changed reporting conventions, but does not explain this particular difference. Its available PDF lacks the natural-habitat size-frequency appendix needed to adjudicate the exact counts. The report lists Pelican survey datesJune20/August9, whereas the legacy size records carryAugust7; that timing difference is not evidence that either frequency table is wrong. No checked coauthor email explains these five measurements. A current database export or documented correction from the data steward is still needed. [1995report](https://irma.nps.gov/DataStore/DownloadFile/485211).

Both versions were tested on the same 37 primary years with the same before/after assignment, model and offset. The current Step model was reproduced exactly before changing 1995. The sensitivity rescales that year's biomass using the two mean individual masses while keeping counts and other measurements fixed; it does not claim a new whole-pipeline bootstrap or model-selection run.

| Scorpion red-urchin biomass effect | Estimate | Conditional 95%interval |
|---|---:|---:|
| Current primary |1.25304|0.70353 to 1.80255|
| Legacy 140-animal sensitivity |1.27010|0.71870 to 1.82150|

Units are the before–after change in the log protected/reference biomass ratio. Difference:+0.01706, about 1.36%of the original log-effect estimate. The positive before–after finding is unchanged. These reproduce the primary model's conditional uncertainty, not total biomass uncertainty. Post-protection recovery slopes are exactly unchanged because 1995 is before protection.

## The four 2023 groups

Scorpion red/purple sizes contain 144/207 animals; Arch Point red/purple contain 151/211. Original physical-site identities and exact measurements remain unverified against an independent 2023 source. Scott's November 2 email explicitly supplies 2012–2022 and confirms those records against the database. It does not validate 2023; we requested the earlier window.

All four 2023 groups have zero uses in the current expanded review size pools. None of their corresponding current primary biomass histories includes 2023. They therefore have **no effect on the analyses currently reported or reviewed**. Keep them as unverified out-of-current-analysis evidence, rather than calling them erroneous or adding them to extend the series. An unrestricted all-years import still fails its guard; this is distinct from their demonstrated lack of use in the current analysis.

## Correspondence and source scope

Read targeted Scott Gabara work/personal and Emily Donham messages about size data, KFM, source filenames and the disputed site/years, including Scott's November2delivery and earlier October discussions. Private email evidence is preserved locally under `private/`. No direct reply resolving the old discrepancies or supplying2023sizes was found in those searches. No clarification email was sent. Original/source copies, dated reports and saved audit inputs were checked; source versioning remains a limitation, not grounds to invent a correction history.

## Implementation and checks

Registry now 1482 groups:1461 verified,21 unverified. Existing exclusions leave five unresolved raw groups (Pelican 1995 plusfour 2023 groups). The first is the only unresolved group identified here with use in a current primary effect. Full raw ingestion remains fail-closed; no threshold or researcher exclusion changed.

- `verify_hare_1987_source.py` verifies102records against20published bins and min/max; source hash stored.
- `build_kfm_identity_bundle.py` incorporates only the corroborated Hare correction and installs the private bundle. `size_identity_utils.R` replaces explicit keyed groups, preserving all other rows.
- `check_remaining_kfm_sources.py` computes mass differences, same-year effect sensitivities and2023use checks. All post-protection annual rows are exactly unchanged.
- `verify_remaining_kfm_effects.R` independently reproduces both Scorpion Step estimates, SEs and intervals.
- `check_kfm_pipeline_guard.R` verifies102Hare source records and correct failure on the five remaining raw groups. Mapping/conflict/source-change tests pass.
- The 487 protected primary files are unchanged. Manuscript, supplement, frozen primary and raw sources remain unchanged; no email or push.

Tables: [mass comparisons](older_mass_comparison.csv), [exact size disagreements](exact_size_disagreements.csv), [published-bin validation](hare_1987_published_bin_check.csv), [primary sensitivity](scorpion_primary_step_sensitivity.csv), [post-protection trends](post_protection_slopes.csv), [2023 use](2023_analysis_use.csv).
