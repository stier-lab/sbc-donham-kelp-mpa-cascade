coverage <- read.csv("tables/table_figure5_annual_membership.csv")
counts <- read.csv("tables/table_figure5_annual_coverage.csv")
pairs <- read.csv("tables/table_figure5_period_pairs.csv")
summary <- read.csv("tables/table_figure5_period_summary.csv")
primary <- read.csv("tables/table_resistance_recovery_sensitivity.csv")
primary <- primary[primary$baseline == "2010-2013*", ]
stopifnot(nrow(coverage) == 16 * 11, nrow(pairs) == 30,
          all(coverage$mpa_present == coverage$reference_present),
          identical(counts$n_pairs, c(9L, 9L, 9L, 11L, 11L, 10L, 9L,
                                     10L, 10L, 10L, 10L, 11L, 9L, 11L, 11L, 11L)))
for (metric in primary$metric) {
  s <- summary[summary$metric == metric, ]
  p <- primary[primary$metric == metric, ]
  z <- pairs[pairs$metric == metric, ]
  ci <- exp(t.test(log(z$inside / z$outside))$conf.int)
  stopifnot(nrow(z) == 10, !anyDuplicated(z$reserve),
            isTRUE(all.equal(s$ratio, exp(mean(log(z$ratio))))),
            isTRUE(all.equal(c(s$lower, s$upper), as.numeric(ci))),
            round(s$ratio, 2) == p$ratio_in_over_out,
            round(s$lower, 2) == p$ratio_CI_low,
            round(s$upper, 2) == p$ratio_CI_high)
}
stopifnot(identical(as.integer(tapply(pairs$ratio > 1, pairs$period, sum)),
                    c(9L, 9L, 8L)))
cat("PASS: annual membership, period cohorts, ratios and intervals agree.\n")
