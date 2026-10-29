# KFM predator sizes: matching and biomass feasibility

29 October 2026

**Sheephead sizes can support a limited follow-up comparison, but the retrieved records do not support replacing the full KFM density history with a reliable biomass history. No compatible lobster size records were found in the inspected KFM releases.** No primary result or manuscript has been changed.

## What we recovered

- The [Kushner et al. 2013 archive](https://doi.org/10.1890/13-0562R.1) contains 1,565 sheephead size records from 2007–2011. Records can represent one fish or a group of fish; they are not all independent individual lengths.
- The [NPS 2015 data release](https://irma.nps.gov/DataStore/Reference/Profile/2259652) contains sheephead size-frequency summaries at 33 sites. These are bins with sample counts and percentages, not individual records or dated survey events.
- Neither that release's natural-habitat size table nor the historical archive's 376,511 invertebrate-size rows contains lobster sizes. Species code 8001 identifies lobster in the KFM handbook. This establishes absence from these releases, not absence from all KFM collections.
- The original supplied MBON files contain sheephead visual-transect counts and lobster band-transect counts. The PISCO size file contains lobster measurements, but those records do not establish a verified KFM site/visit match. We did not borrow sizes from another reef or pool years to manufacture coverage.

Public 2013, 2014, 2016 and 2017 NPS appendix downloads were identified but returned server errors. An alternate NPS host returned invalid placeholder files, which were rejected and quarantined. The newer 1982–2025 package is publicly described but its downloads are marked Internal. Those records remain unexamined; unavailable downloads are not biological absences.

## Sheephead: what matches the current paper

Matching retained the paper's four KFM sheephead comparisons, original MPA/reference assignments, and annual observation windows. NPS numeric site identifiers link the historical lengths to MBON counts. The 2015 workbook's natural-habitat table supplies its own site-number/name crosswalk. No fuzzy name match was used.

Across the four comparisons, there are **423 contributing count site-years in 109 paired MPA-years**. Sizes match **118 site-years**, representing **1,766 measured fish**. All 98 historical matched site-years have size-survey dates matching count dates. The remaining 20 matches are from 2015, whose published size summaries provide only the year. Season and exact visit cannot be checked for those summaries.

| Current comparison | Site-years with sizes / count site-years | Median fish measured per matched site-year (range) | Matched samples with fewer than 10 fish |
|---|---:|---:|---:|
| Anacapa Island SMR 2003 | 24 / 68 | 19 (4–105) | 7 |
| Harris Point SMR | 16 / 89 | 3 (1–10) | 13 |
| Scorpion SMR | 42 / 164 | 12 (1–55) | 13 |
| South Point SMR | 36 / 102 | 10 (2–32) | 16 |

All four comparisons have sizes on both sides in **2007–2011 and 2015**. All count sites have sizes in 23 of these 24 MPA-years. The exception is Wyckoff Ledge in the Harris Point comparison in 2010: its transects recorded zero sheephead. A size distribution is not needed to retain its observed biomass point estimate of zero. Keeping that zero-count site preserves the comparison; dropping it would change site weights. This does not establish a true population absence or eliminate sampling uncertainty.

No recovered size observation predates the 2003 designation. These records cannot support a KFM sheephead biomass before/after analysis. They also cannot fill the later years in the current histories. There are 255 positive-count site-years without retrieved sizes. Missing sizes must not become zero biomass.

## How reliable would sheephead biomass be?

**The limiting issues are representativeness and sparse size samples, not the number of bootstrap repetitions.**

1. **Survey footprint differs.** The documented size protocol searches broadly around the permanent transect during or after a roving survey, whereas the paper uses visual-transect counts. Total length is measured in centimeters, so the length definition is compatible with the adopted sheephead equation. Matching site and date does not establish equal depth, habitat or size-selective detectability. The [KFM methods in BOEM 2019-063, printed p. 26](https://www.boem.gov/sites/default/files/documents/regions/pacific-ocs-region/environmental-science/BOEM-2019-063.pdf) describe the search area, prioritization of species and changes in sampling completeness. Sheephead stages are priority 1 in the historical size records.
2. **Many pools are small.** Of the 118 matches, 49 have fewer than 10 measured fish and 88 have fewer than 20. These are descriptive counts, not a universal sample-size cutoff. Harris Point is particularly sparse.
3. **Large fish strongly influence inferred mass.** As a diagnostic, we converted recorded lengths using the paper's equation, W = 0.0144 × TL^3.04. In 38 matched pools, one source record supplies more than half the inferred total mass. A record can contain multiple fish. In 68 pools, a conditional bootstrap interval for mean individual mass has a half-width exceeding 50% of its estimated mean. Six pools have insufficient size variation to estimate a meaningful bootstrap interval; they are flagged rather than reported as precisely known.
4. **Those intervals omit important uncertainty.** Historical diagnostics resample original size records, preserving grouped counts. They do not recover within-group length variation where only a group mean is recorded, and they cannot resample independent dives without dive identifiers. In the matched pools, 187 grouped records lack minimum/maximum lengths. The 2015 diagnostic samples the reconstructed bin frequencies at bin midpoints; explicit lower/upper bin-mass bounds are saved separately. Neither calculation includes visual length error, allometric parameter uncertainty, count-survey variation or a mismatch between the sampled habitats. They are not confidence intervals for biomass density.
5. **Life stages need a compatible interpretation.** The size records distinguish female, male and juvenile sheephead. The MBON species-level counts do not retain those categories. The 2015 size release labels some fish over 10 cm as juvenile, so that label must not be silently equated to a fixed <10-cm or young-of-year class. Seventeen matched site-years have sizes from roving searches despite zero visual-transect counts; that is another reason to preserve the distinction between the two samples.

A defensible next numerical step would be a separately labeled, restricted-period sensitivity comparison after confirming that the size survey represents the fish counted on transects. It would retain site-specific sizes, retain zero-count sites, and propagate both count and size uncertainty. The current evidence does not justify promoting those provisional estimates into the primary biomass analysis.

## Lobster

The five retained KFM lobster comparisons contain **436 contributing count site-years**, including **208 with positive counts**. None has a verified size match in the retrieved KFM size releases. Bootstrapping cannot solve that gap. Converting these counts would require measurements of lobster carapace length, with units, dates, site identifiers and sampling methods that can be matched to the band counts. Borrowing a PISCO or mainland size distribution without establishing compatibility would introduce a new assumption rather than recover missing KFM measurements.

## What remains to obtain or confirm

Scott Gabara is the appropriate contact for:

- The current fish-size export, including dates, site IDs, lengths, counts, life-stage categories and observer/dive identifiers if retained; whether those samples can represent the visual-transect count population across habitat and depth.
- Whether lobster carapace lengths were collected separately, and, if so, their years, sites, units and sampling protocol.

These are data/measurement questions. The existing permission request remains separate. No new email was sent during this audit.

## Reproduction and checks

Run `python3 code/audit_kfm_size_feasibility.py` from the analysis repository (numpy and openpyxl required). The default source cache is the sibling data repository's `data/source_archive/kfm_size_feasibility_2026-10-29`; set `KFM_SIZE_CACHE` to override it. Its download manifest records URLs and SHA256 hashes. Inputs are read-only.

Outputs are `tables/table_s_kfm_size_{coverage,site_year_matches,pools,count_sources,count_reconstruction,inputs}.csv` and `table_s_kfm_size_audit.json`. The calculation reconstructs all 218 retained sheephead MPA/status/year density values from the underlying count rows and reproduces the canonical values to within 1e-10. The 235 identical historical sheephead rows are retained: no individual IDs establish that they are erroneous duplicates rather than repeated same-sized fish. Annual sample counts repeated across histogram bins are counted once per site/stage, not summed over bins. Blank/missing measurements never become zeros. Primary input and results checksums remain unchanged.

Two misleading comments in the data-processing script were corrected: it now distinguishes the count-only input from the existence of separate size surveys, and the visual/RDFC taxon identifiers agree with the actual source file. Executable processing logic is unchanged.
