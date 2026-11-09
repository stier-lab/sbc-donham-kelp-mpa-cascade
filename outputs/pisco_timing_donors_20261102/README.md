# Survey-date pairing and size-donor checks

This review addresses the six PISCO animal comparisons at Abalone Cove and Point Vicente. It preserves the researcher exclusions, existing station/year membership, all-month eligibility and constant shared INNER/MID zones. No manuscript, primary result, pooled analysis or size imputation was changed.

## Survey timing

The files provide survey dates, not verified campaign identifiers. We therefore tested explicit date matching rather than claiming to recover field campaigns. Within each protected/reference station pair, year and zone, match the largest possible number of visits under each date tolerance, then minimize their total gap. A visit is used once within that station-pair/zone/year. All original station pairs and both zones must remain represented for an annual contrast. Shared stations recur across pairs; they are not independent samples. Pair means receive equal weight within the annual role means.

Tolerances0,14,30,60,90 and366days are a sensitivity grid, not adopted biological cutoffs. The366day version means nearest matching within the observed calendar year; it does not establish seasonal comparability. Dates and all4206matched-pair rows across the scenarios are recorded.

At Abalone Cove, complete coverage retains0,0,1,1,2 and11post-protection years across that grid. Point Vicente retains0,0,0,1,6 and10years. These counts apply to each of the three species in each location. The loss comes from requiring every reference station and both zones, not evidence that every individual protected/reference visit is far apart. We did not improve coverage by selecting whichever references give a preferred result.

Point Vicente's90day check has six years: redurchin slope−0.024, lobster−0.104 and purpleurchin−0.040log-ratio units/year, with all conditional temporal intervals spanning zero. The purpleurchin comparison is especially informative: the all-month ten-year slope is−0.130; restricting the original values to those six years changes it to−0.051, and date pairing then changes it to−0.040. Most of that change comes from the years retained, not pairing itself. Abalone's90day complete record is too short to fit a recovery slope. All scenarios include same-year baseline fits where estimable.

Date-matching figure: A=Abalone Cove; B=Point Vicente. Each line represents the common coverage of the three species. The366day endpoint retains broad timing differences and is not equivalent to close pairing.

## Can other measurements predict local body mass?

We held out624positive-count site/zone/date species samples with actual size measurements. Counts of measured animals are preserved in the mean-mass conversion. Cross-site prediction excludes the entire target site from the same-year donor data; the held-out site cannot predict itself through another zone or date. Donors are selected by species, zone, protection role and date, never by their agreement with the target. The nearest eligible visit within each donor site supplies a frequency-weighted mean; donor sites receive equal weights.

Five strategies are evaluated: another visit at the same site/zone; another zone at the same site/date; another site with the same protection role/zone; another protection role/zone as a diagnostic; and other same-role sites with zones pooled. Cross-role borrowing is not adopted. Date tolerances30,60,90 and366days are diagnostics. There are5860available predictions across all strategies/tolerances, not5860independent validation samples. Like-for-like target comparisons and separate location/protection/sample-size tables prevent comparing accuracy without also checking coverage.

For other sites with the same protection role and zone within90days, predictions are available for158redurchin,204purpleurchin and125lobster samples. Typical multiplicative errors are1.39,1.33 and1.48fold, respectively. Only41%,46% and40%fall within25%of the measured mean mass. These are descriptive validation summaries, not preselected acceptance criteria or calibrated prediction intervals. Target sample sizes are often small: median5redurchins,13.5purpleurchins and3lobsters. Prediction error therefore includes sampling noise in the observed target mean and cannot be attributed entirely to spatial biological differences.

Errors vary by location and protection role. At Point Vicente, lobster predictions for protected stations have median1.79fold errors; at Abalone reference stations, the worst tenth of lobster errors exceed roughly4.3fold. Abalone has only one protected station in this roster, so same-role cross-site prediction cannot be validated for its protected samples. Results at its references do not establish accuracy inside that MPA. Average errors near zero would not make individual predictions precise; some strategies also have reciprocal paired observations.

Same-site other-visit validation is sparse: only28eligible held-out predictions exist, all with date gaps greater than90days. This is inadequate evidence for estimating a reliable within-site seasonal relationship on its own. Same-date other-zone borrowing is also imperfect (typical1.39–1.58fold errors depending on species). Neither nearby dates nor an adjacent depth zone automatically guarantees representative sizes.

Prediction figure: A=redurchin, B=purpleurchin, C=lobster. Same-role, same-zone, other-site donors within90days. The diagonal denotes exact prediction. Both axes are logarithmic; each point is a held-out sample mean, not an individual animal. Allometric coefficients are held fixed.

## What this could recover

Of153positive visits lacking eligible local sizes, same-role/zone donors exist for77within30days,106within60days,120within90days and133somewhere in the same year. Donor availability is not validated biomass recovery. These counts are not additional complete annual histories: one unsupported station/visit can still prevent a complete protected/reference biomass contrast. We have not imputed these visits or rerun biomass using unvalidated predictions.

## Recommendation

Keep the all-month density histories available. A narrow complete-reference date window would remove much of the time series, so do not adopt one simply to create a nominally matched analysis. Retain the date-pair sensitivity alongside the full record, and consider a survey-level model that separates seasonal timing from site and year only after checking whether the sampling calendar supports those terms.

For biomass, retain directly matched measurements and keep borrowed-size results explicit. These validation errors support estimating and propagating borrowing uncertainty rather than substituting donor means as known values. A hierarchical size model is a possible next step, but it needs held-out site/year validation, realistic prediction intervals and explicit limits for locations such as Abalone's single protected site. The current checks do not establish that such a model will recover the missing biomass reliably.

## Verification

Independent R calculations reproduce54timing model coefficients/standard errors and all5860donor predictions. Date pairs have unique within-pair visit use and verified gap bounds; cross-site donors exclude the target site. All487protected manuscript/result hashes are unchanged. Conditional temporal intervals do not include total measurement uncertainty. Scripts: `code/tools/check_pisco_timing_donors.py`, `code/tools/verify_pisco_timing_donors.R`, `code/tools/report_pisco_timing_donors.py`.
