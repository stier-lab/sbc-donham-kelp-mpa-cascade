source("code/replications/malakhoff2021/helpers.R")
stopifnot(all.equal(urchin_mass(50, "red"), .00059 * 50^2.917),
          all.equal(urchin_mass(50, "purple"), .00059 * 50^2.870),
          size_class(35, "code") == "small", size_class(35, "paper") == "medium")
suppressPackageStartupMessages(library(data.table))
d <- fread("outputs/malakhoff2021/baci_input.csv")
stopifnot(!any(d$site %in% c("a-k-12", "a-k-13")),
          all(sort(unique(d$site)) == sprintf("a-k-%02d", c(1:11,14:16))),
          !anyDuplicated(d[, .(site, year, species)]),
          all(d[density == 0]$biomass == 0),
          all(is.na(d[missing_positive_sizes == TRUE]$biomass)))
x <- d[is.finite(mean_mass)]
stopifnot(all(abs(x$biomass - x$density * x$mean_mass) < 1e-8),
          all(x$small_mass_corrected <= x$biomass + 1e-8))
r <- fread("outputs/malakhoff2021/model_results.csv")
stopifnot(!any(r$eligible_for_manuscript),
          all(r[formulation == "original_year_site"]$status %in% c("failed", "not_temporal_AR1")),
          all(r[status == "fitted"]$n_reserve_sites == 4))
cat("PASS: cohort, allometry, zeros, joins, size definitions and quarantine flags.\n")
membership <- fread("outputs/malakhoff2021/common_support_membership.csv")
support <- fread("outputs/malakhoff2021/common_support_counts.csv")
kept <- membership[included == TRUE, site]
stopifnot(all(support[site %in% kept, before] >= 2),
          all(support[site %in% kept, after] >= 2),
          !any(c("a-k-05", "a-k-07", "a-k-14") %in% kept))
common <- fread("outputs/malakhoff2021/common_support_results.csv")
stopifnot(!any(common$eligible_for_manuscript), all(common$status == "fitted"),
          all(common[omitted_reserve == "none", n_sites] == length(kept)),
          all(common[omitted_reserve != "none", n_sites] == length(kept) - 1),
          all(common$lower <= common$estimate & common$upper >= common$estimate),
          all(is.na(common[formulation == "site_independent", phi])))
inputs <- fread("outputs/malakhoff2021/common_support_inputs.csv")
stopifnot(!anyDuplicated(inputs[, .(site, year, endpoint)]))
cat("PASS: fixed before-after support, leave-one-reserve-out membership, uncertainty and independent-error fits.\n")
report <- paste(readLines("docs/MALAKHOFF_2021_AVAILABLE_DATA_RESULTS.md"), collapse = "\n")
for (focal in c("red_biomass", "purple_biomass", "total_biomass")) {
  for (last in c(2019L, 2022L)) {
    z <- common[endpoint == focal & end_year == last & first_post == 2005 &
                  omitted_reserve == "none" & formulation == "site_calendar_AR1"]
    stopifnot(nrow(z) == 1)
    text <- sprintf("%.2f (%.2f, %.2f)", z$estimate, z$lower, z$upper)
    stopifnot(grepl(text, report, fixed = TRUE))
  }
}
cat("PASS: report's main estimate table matches regenerated results.\n")
source("code/replications/malakhoff2021/test_bridge.R")
