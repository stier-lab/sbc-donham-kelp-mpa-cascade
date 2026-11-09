# Scott’s fish sizes and Li’s updated LTER biomass

November 2, 2026. The new deliveries are preserved and reconciled. The paper, primary inputs, exclusions and fitted results are unchanged. Scott’s and Emily’s recent emails were checked; no replies were sent.

## Recommendation

Use Scott’s export as the source for a restricted-period sheephead biomass review, with count-visit matching and explicit uncertainty. It does not justify replacing the full KFM density history with biomass. Keep our LTER adult-stipe response and survey-month conversion; the updated provider biomass is a different aggregate and should not be substituted directly.

## KFM sheephead: much better coverage, but still a different sampling method

The received CSV has 36,516 records across species, including 5,658 sheephead records representing 10,073 fish at 33 named sites in 2007–2025. All 33 names match the existing physical-site dictionary without fuzzy matching. Counts and representative lengths are positive; grouped fish records retain their frequency weights. Rows are not all independent fish measurements or replicate surveys.

| Current comparison | Count site-years | Earlier size support | New size support |
|---|---:|---:|---:|
| Anacapa Island SMR 2003 | 68 | 24 | 56 |
| Harris Point | 89 | 16 | 40 |
| Scorpion | 164 | 42 | 102 |
| South Point | 102 | 36 | 84 |
| Total | 423 | 118 | 282 |

These are the existing four KFM sheephead histories, not a newly selected cohort. The full export’s later years are inventoried but do not extend our count histories automatically. Existing fish membership is retained; historical benthic exclusions must not be indiscriminately applied to fish. Scott’s new reference concerns remain a separate ecological review.

- Of the 282 matched site-years, 279 have all size-survey dates matching count-survey dates. South Point has additional size visits in 2014 and 2015; Potato Pasture’s 2016 size survey was one day before its count survey. The exact-visit table uses matching dates first and leaves that one-day case explicit rather than silently declaring an exact match.
- Complete comparison-years with some measured size information increase from 24 to 57. Three additional comparison-years have all observed counts zero and no sizes; these are not evidence of a measured size distribution or an estimable uncorrected response ratio.
- There are still 110 positive-count site-years without sizes, and 100 matched site-years have fewer than ten measured fish. These are descriptive support counts, not new sample-size eligibility thresholds.
- No size record predates 2003. A biomass before–after analysis of the 2003 designation is still unsupported.
- Scott confirms that sizes came from one diver’s 30-minute roving survey. Our counts are visual-transect observations. Same site/date does not prove equivalent depth, footprint or detection. He also confirms KFM does not measure lobster carapace lengths; the fish file cannot resolve lobster biomass.

### The historical overlaps are not identical

The 2007–2011 comparison uses site, date, stage, count and representative length as a multiset, preserving repeated rows. There are 269 differing signatures: 220 excess historical rows and 52 excess new-export rows. These differences mostly concern repeated rows in 2011, with additions in 2007 and 2010. This does not prove that identical historical rows were errors, so neither source was deduplicated or overwritten.

Five overlapping early site-years in the current count roster have changed mean-mass estimates. The 20 overlapping 2015 site-years also change because the earlier calculation used histogram-bin representatives, while the new export supplies finer measurements and, at some sites, additional fish. These are changes to provisional size-derived mass estimates; KFM sheephead biomass is not part of the existing primary analysis.

[Coverage](fish_coverage_summary.csv) · [Historical differences](fish_legacy_differences.csv) · [Pool comparison](fish_pool_comparison.csv) · [Exact-visit support](fish_exact_visit_support.csv) · [Count/size matches](fish_count_size_matches.csv)

## LTER: unchanged adult observations, different biomass definitions

Retrieved SBC50.21, published November 2, 2026, with the file `Annual_All_Species_Biomass_at_transect_20261012.csv`. Its metadata points to SBC127.5 coefficients and SBC18.31 adult kelp observations. The filename’s October date is not the release date.

The 29,526 adult kelp records through 2023 in SBC18.31 exactly match the current March 2024 delivery on year, month, date, site, transect, quadrat, side, area and frond count, including repeated rows. The updated biomass table covers 1,100 kelp transect-years; 968 match the current raw-data reconstruction and 132 occur in 2024–2026. Those extra years are not automatically added.

The provider’s `MAPY` density is **adult stipe density plus juvenile kelp density**, not adult stipes alone. This sum reproduces the finite matched density values to floating-point tolerance. Juveniles are recorded separately as `MPJ` in the quadrat data. Mixing this provider density with adult-only diver stipe counts would change what the response means.

Our code converts adult stipes using the observed month’s fixed SBC127.5 Model1 coefficient. Where juvenile density is zero and adult density is positive, the provider dry mass instead implies **0.10234 kg per adult stipe**, the July coefficient in that same table. Its linked protocol describes August. The discrepancy needs clarification; this audit does not infer the provider’s internal code from the coefficient alone.

Adding a constant 1.16 g per juvenile to that fixed adult calculation reproduces 738 of 966 finite matched biomass values, but leaves 228 unmatched, with a maximum residual of 22.33 g/m². Juvenile size information varies through time. We have isolated the remaining component but have not independently reproduced its full size-to-mass treatment. Do not call the entire provider formula verified.

Two matched provider kelp entries are unavailable (`-99999`): Naples transect6 in2000 and Diablo transect3 in2004. Adult stipe records exist for those entries. The reconciliation keeps the provider missing values unknown and does not convert them to zero or remove their valid adult observations. Any provider-based role-year calculation remains unavailable when a required constituent is missing.

### What would change if we substituted the provider column?

On matching observations, provider versus current biomass differs by more than 0.000001 g/m² in 494 of 966 finite transect comparisons. This is largely a response-definition/calibration comparison, not evidence that 494 raw surveys were wrong.

Using the same stations and years, positive protected/reference ratios would differ by a median absolute log difference of 0.0043 at Campus Point and 0.0024 at Naples; maximum absolute differences are 0.2465 and 0.1792. The largest calibration-only difference is 0.2458 at Campus Point and 0.0424 at Naples. These diagnostics use positive values only, with no invented zero offset, and are not fitted protection effects or updated confidence intervals.

The [focal-taxon release comparison](lter_provider_focal_taxa_revisions.csv) separately records changes from provider release20 to21. Provider transect biomass must not be equated to this paper’s raw-based animal estimates, especially where size pooling, mass units or survey coverage differ.

[Transect comparison](lter_transect_comparison.csv) · [Formula decomposition](lter_formula_decomposition.csv) · [Role-years](lter_role_year_comparison.csv) · [Positive pair ratios](lter_positive_pair_ratio_comparison.csv) · [Summary](lter_ratio_change_summary.csv)

## Emails and remaining decisions

Scott’s November2 email supplies the fish file, confirms the roving protocol and absence of lobster lengths, and raises response-specific reference questions. Emily’s latest email reiterates three before observations at the same physical stations, assessed separately for fish and benthic surveys, with the implementation year counted as before. No baseline convention or researcher exclusion was changed here. No newer relevant reply appeared in the targeted searches.

Remaining work is to resolve KFM count/size representativeness, test the restricted-period biomass with joint uncertainty, and clarify the LTER aggregate’s fixed-month description and juvenile calculation if that product is to be used. Existing raw-based monthly adult kelp estimates remain the appropriate like-for-like response for the current review. KFM redistribution permission remains unresolved.

## Reproduction

Run `python3 code/tools/reconcile_scott_fish.py`, `python3 code/tools/reconcile_lter_annual_biomass.py`, then `python3 code/tools/decompose_lter_reconciliation.py` from the analysis repository. Review inputs are pinned in the sibling data repository’s ignored `data/source_archive/fish_lter_reconciliation_20261102/`, with an [input manifest](input_manifest.json). The older SBC50.20 comparison uses the existing NAS cache; this dependency is explicit. The juvenile-source subset is an extraction of SBC19.33 `MPJ`, not an independent delivery. Full source records remain private; no EDI credentials or email bodies are in the report.
