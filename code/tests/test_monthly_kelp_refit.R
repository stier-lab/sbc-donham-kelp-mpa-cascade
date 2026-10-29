# Independently reconstruct headline period comparisons from promoted raw data.
raw <- read.csv("data/harmonized/harmonized_raw_responses.csv")
d <- subset(raw, taxon_name == "Macrocystis pyrifera" & resp == "Bio")
d <- aggregate(value ~ CA_MPA_Name_Short + status + year, d, mean)
headline <- subset(read.csv("tables/table_resistance_recovery.csv"),
                   taxon == "Macrocystis pyrifera")
sensitivity <- read.csv("tables/table_resistance_recovery_sensitivity.csv")
counts <- read.csv("tables/table_resistance_recovery_pair_counts.csv")
baselines <- list("2008-2013" = 2008:2013, "2009-2013" = 2009:2013,
                  "2010-2013*" = 2010:2013, "2011-2013" = 2011:2013,
                  "2012-2013" = 2012:2013, "2013" = 2013)
windows <- list(resistance = 2014:2016, recovery = 2017:2019,
                recovery_recent = 2020:2023)
same <- function(a, b) stopifnot(isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 1e-10)))
for (baseline in names(baselines)) for (metric in names(windows)) {
  b <- aggregate(value ~ CA_MPA_Name_Short + status,
                 subset(d, year %in% baselines[[baseline]]), mean)
  m <- aggregate(value ~ CA_MPA_Name_Short + status,
                 subset(d, year %in% windows[[metric]]), mean)
  x <- merge(b, m, by = c("CA_MPA_Name_Short", "status"), suffixes = c("_base", "_period"))
  x$ratio <- x$value_period / x$value_base
  x <- subset(x, is.finite(ratio) & ratio > 0)
  a <- subset(x, status == "mpa", c(CA_MPA_Name_Short, ratio))
  z <- subset(x, status == "reference", c(CA_MPA_Name_Short, ratio))
  paired <- merge(a, z, by = "CA_MPA_Name_Short", suffixes = c("_mpa", "_ref"))
  delta <- log(paired$ratio_mpa) - log(paired$ratio_ref)
  row <- sensitivity[sensitivity$baseline == baseline & sensitivity$metric == metric, ]
  stopifnot(nrow(row) == 1, nrow(paired) == row$n_mpa)
  ci <- t.test(delta)$conf.int
  p <- suppressWarnings(wilcox.test(delta)$p.value)
  same(row$ratio_in_over_out, round(exp(mean(delta)), 2))
  same(row$ratio_CI_low, round(exp(ci[1]), 2))
  same(row$ratio_CI_high, round(exp(ci[2]), 2))
  same(row$p_wilcoxon, signif(p, 3))
  loo_p <- vapply(seq_along(delta), function(i) suppressWarnings(wilcox.test(delta[-i])$p.value), 0.0)
  loo_ratio <- vapply(seq_along(delta), function(i) exp(mean(delta[-i])), 0.0)
  same(row$loo_p_max, signif(max(loo_p), 3))
  same(row$loo_ratio_min, round(min(loo_ratio), 2))
  if (baseline == "2010-2013*") {
    h <- headline[headline$metric == metric, ]
    same(h$inside, median(paired$ratio_mpa))
    same(h$outside, median(paired$ratio_ref))
    same(h$paired_p, signif(p, 3))
    same(counts[counts$metric == metric, ]$n_inside_gt_reference,
         sum(paired$ratio_mpa > paired$ratio_ref))
  }
}
cat("PASS: independently reconstructed all 18 baseline-period comparisons and three headline medians/counts.\n")
