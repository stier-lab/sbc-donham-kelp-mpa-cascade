# Coverage-defined MPA comparisons

Run `Rscript code/R/50_before_after_cohorts.R` from the analysis root. This module is also wired into `run_all.R` and `run_review_concerns.R`. It reads current harmonized annual ratios and the primary 145-history cohort, validates the reference registry, and never changes primary inputs or fitted results.

## Question and grouping

Do protected-reference differences depend on whether an MPA has observations before protection? Define coverage separately for each MPA, taxon, response and source view. All contributing histories require at least three paired post-designation years. The main grouping requires at least three paired pre-designation years in at least one source. The complementary group has fewer than three in every source. These groups are mutually exclusive within each endpoint/view. Short records are not labelled as zero pre-designation coverage: the roster and coverage tables identify zero versus one or two years. A second grouping uses one pre-designation year to implement the literal any-before versus none comparison.

Coverage is not the primary model label. Some primary CI records have before observations; passing a pre-trend screen is not required here. Grouping is based only on recorded coverage, not effect direction or statistical significance.

## Comparisons

- Post-designation difference: mean annual ln(MPA/reference), averaged across all available sources within MPA, then equally across MPAs. This is the same calculation in both coverage groups.
- Before-after change: subtract the pre-designation mean within each eligible source, average eligible sources within MPA, then average MPAs. This uses the 60 histories meeting the three-year rule.
- Post difference, eligible sources: same histories as the before-after estimate, retaining the existing module-27 method comparison. Do not subtract the all-source post estimate from the eligible-source change: source composition can differ.

Source views include all programmes, divers, Landsat, and individual programmes. Group labels are recalculated within each view, so source-stratified MPA counts need not add to all-programme counts. Regional views retain separate island/mainland results. The common calculation does not make calendar periods, reserve ages, geography or habitat identical between groups.

Approximate 95% intervals use the established delete-shared-reference-group jackknife. The between-cohort difference deletes connected groups jointly across cohorts, including any shared references linking the cohorts. No difference interval is reported unless both cohorts contain at least three reference groups. No significance tests or adjusted P values are claimed: these are exploratory descriptive checks. Known reference sharing is accounted for, but satellite polygon overlap and all spatial dependence are not known. Unequal group sizes and few groups limit interval calibration.

## Findings

Using the three-year rule, kelp has 17 MPAs with adequate before coverage and two without. Mean post-designation log ratios are 0.459 and 1.746. Only two MPAs support the latter, so there is no between-cohort uncertainty interval. Both cohorts contain mainland and island reefs, but source composition differs. All 16 Landsat MPAs have adequate pre-designation data. With divers alone, the groups are six versus five MPAs, and means are 0.694 versus 0.790; their difference is -0.095 with approximate interval [-1.714, 1.524]. These are uncertain estimates, not evidence of equivalence.

The same-history kelp change remains 0.693 [0.168, 1.218] across 17 MPAs. Lobster density has a positive before-after change (1.488 [0.807, 2.168], six MPAs), but lobster biomass has no histories meeting the three-before-year threshold. Sheephead biomass has only four eligible MPAs, all mainland. Purple-urchin biomass has only one. The six estimable all-programme between-cohort intervals span zero; the remaining three endpoints lack sufficient cohort replication. The one-before-year sensitivity also leaves every estimable between-cohort interval spanning zero, with several point-estimate directions changing. Neither split establishes that cohort membership or protection causes the differences.

## Outputs and checks

`table_s_coverage_cohort_{series,roster,pairs,summary,difference,coverage,inputs,status}.csv` gives source-level coverage, exclusive membership, MPA estimates, group summaries, with-minus-without comparisons, composition, fingerprints and completion status. `fig_s_coverage_cohorts` shows individual post-designation MPA contrasts, equal-MPA mean bars and group sample sizes; it contains no statistical significance symbols. The result note belongs in the private manuscript repository.

Run `Rscript code/tests/test_before_after_cohorts.R` after generating outputs. It checks complete/unique grouping, source aggregation, agreement with module 27 for eligible histories, missing intervals for sparse cohorts, literal no-before classification, shared-reference deletion, and unchanged input fingerprints. This grouping does not replace primary manuscript estimates or resolve siting bias.
