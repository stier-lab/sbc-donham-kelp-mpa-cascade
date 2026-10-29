# Satellite siting checks

## Finding

The five-year comparison shows a larger average kelp increase relative to
references at designation than at the two earlier placebo dates. This pattern
does not hold in the three-year comparison. Pre-designation satellite histories
usually track positively, but some pairs track poorly and trend estimates depend
on the preperiod. These findings neither establish a protection effect nor show
that siting explains the observed kelp advantage.

## Fixed-pair comparisons

Values are equal-pair means of changes in annual log inside/reference contrasts,
using the existing recorded offset. Each row uses the same pairs at all three
anchors. The designation year and each placebo anchor year are omitted.

| Years per side | Pairs | Real designation | Near placebo | Far placebo |
|---|---:|---:|---:|---:|
| 3 | 15 | -0.009 | 0.721 | -0.613 |
| 5 (primary diagnostic) | 15 | 0.427 | -0.110 | -0.228 |
| 7 | 7 | 0.577 | 0.164 | -0.245 |

The three- and five-year cohorts contain the same fifteen pairs; Farnsworth
lacks sufficient earlier coverage. Seven-year comparisons retain only seven
mainland pairs, so they are not a network-wide replication. Only seven of
fifteen pairs have a five-year actual change exceeding both placebo changes.
There is substantial variation hidden by the pooled means.

Changing the response definition does not reverse the five-year ordering:
log ratios of unshifted period means give 0.387 at designation, -0.314 at the
near placebo and -0.076 at the far placebo, on the same fifteen pairs. For
three-year windows, the corresponding values are 0.178, 0.242 and -0.463.
Therefore the near-zero three-year annual-log result is not a universal absence
of change; its magnitude depends on averaging and transformation. Unshifted
index-difference changes are exported separately. Indices are not physical
biomass units, and these local-window changes are not the paper's selected
primary effects or its heatwave-resilience estimates.

No P values are attached to the placebo comparison. Anchors are not randomized,
windows overlap, years are autocorrelated, and shared references may introduce
dependence. Windows are measured relative to each MPA's designation year, so
island and mainland comparisons occur in different calendar years. Mainland
post-designation windows include the later heatwave; this does not isolate
protection from geographically different environmental responses.

## Pre-designation tracking

Across sixteen primary satellite pairs, all full-preperiod detrended index
correlations are positive; their median is 0.645. Medians are 0.796 over the
last ten years and 0.880 over the last five. Short-window correlations have few
observations and are descriptive, not evidence of equivalence or parallel
counterfactual trends. Carrington Point (0.130), Cabrillo (0.147), Abalone Cove
(0.257) and Skunk Point (0.288) track weakly over the full preperiod.

Six of sixteen pairs have positive full-preperiod log-contrast slopes; ten
have positive slopes in each shorter window. The sign balance depends on the
period chosen. This argues against imposing a single long-term pretrend as a
known counterfactual. Strong tracking can coexist with different abundance
levels, changing ratios or different responses to an extreme event.

These satellite results differ from the sparse KFM non-fished-species check.
Do not transfer the non-fished comparison's weak tracking to all satellite
reference pairs, or use satellite tracking to validate the KFM pairs.

## Earlier disturbance

Module 28's existing same-pair comparison is retained in
`table_s_siting_historical_existing.csv`, not refitted or counted as new evidence.
Its network changes in resistance and recovery remain imprecise (95% intervals
include zero). The 1997-1998 event also falls within some placebo windows,
so the historical and placebo checks are not independent corroboration.

## Reproduction and QA

Run `Rscript code/R/46_satellite_siting.R`, then
`Rscript code/tests/test_satellite_siting.R`. Module 46 is also called from
`run_all.R` after module 45. The main pipeline was not rerun for this diagnostic.

The module verifies exact reconstruction of the stored annual log ratios and
guards primary inputs and main figures with before/after MD5 checks. Tests cover
constant and stepped histories, coverage rejection, zero period means,
undefined correlations, pre-designation placebo boundaries and fixed membership.
All passed. Input hashes, exclusions, years, contrasts and raw unshifted index
histories are saved under `tables/table_s_siting_*`. Index and log-history
figures are under `plots/fig_s_siting_*`, without titles or subtitles.

The source-data normalization is retained; the unshifted annual index is not an
absolute canopy-area measurement. No masks, reference assignments, primary
effects or circulated manuscript files were changed. Reference geometry and
cross-pair overlap remain separate unresolved issues. This diagnostic should
not be described as resolving reserve siting.
