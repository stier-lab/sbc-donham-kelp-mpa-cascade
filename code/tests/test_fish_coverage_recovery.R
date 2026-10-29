source("code/R/fish_recovery_utils.R")
source("code/R/fish_siting_utils.R")
source("code/R/primary_uncertainty_utils.R")
source("code/R/review_concern_utils.R")
library(dplyr)

roster <- data.frame(site = c("stable", "late", "short"))
fp <- rbind(data.frame(site = "stable", year = 2000:2015),
  data.frame(site = "late", year = 2008:2015), data.frame(site = "short", year = 2000:2001))
a <- fish_recovery_select(fp, roster, 2000:2015, 2005, "before_after", .8)
b <- fish_recovery_select(fp, roster, 2000:2015, 2005, "post_trend", .8)
stopifnot(identical(a$selected, c(TRUE, FALSE, FALSE)), identical(b$selected, c(TRUE, TRUE, FALSE)),
  !fish_recovery_eligible(2006:2013, 2005, "post_trend"),
  fish_recovery_eligible(c(2006:2012, 2014), 2005, "post_trend"))
read <- function(x) read.csv(paste0("tables/table_s_fish_recovery_", x, ".csv"))
if (file.exists("tables/table_s_fish_recovery_status.csv")) {
  inputs <- read("inputs"); guard <- read("protected")
  stopifnot(identical(inputs$md5, unname(tools::md5sum(inputs$path))),
    identical(guard$before_md5, guard$after_md5),
    identical(guard$after_md5, unname(tools::md5sum(guard$path))))
  cov <- read("coverage"); roster <- read("roster"); site <- read("site_annual")
  old <- read.csv("tables/table_s_fish_siting_coverage.csv")
  benchmark <- filter(cov, scenario == "primary_years|bottom_mid|full|post_trend")
  benchmark <- benchmark[match(old$MPA, benchmark$MPA), ]
  stopifnot(identical(benchmark$eligible, old$post_trend_eligible),
    identical(benchmark$n_after, old$n_after))
  for (i in which(cov$eligible)) {
    row <- cov[i, ]; yy <- as.integer(strsplit(row$years, ";", fixed = TRUE)[[1]])
    rr <- filter(roster, scenario == row$scenario, MPA == row$MPA)
    stopifnot(all(c("mpa", "reference") %in% rr$status), !row$start %in% yy,
      fish_recovery_eligible(yy, row$start, row$endpoint))
    for (s in rr$site) stopifnot(all(yy %in% site$year[site$site == s & site$frame == row$frame]))
  }
  d <- read("density")
  stopifnot(!anyDuplicated(d[c("scenario", "MPA", "species", "size", "year")]),
    all(is.finite(d$mpa)), all(is.finite(d$reference)))
  counts <- d %>% count(scenario, MPA, size, year)
  stopifnot(all(counts$n == 9L))
  annual <- read("annual")
  w <- annual %>% tidyr::pivot_wider(names_from = series, values_from = contrast)
  stopifnot(all(abs(w$Difference_warm - (w$Targeted - w$Nontargeted_warm)) < 1e-10))
  same <- read("common_year_comparison") %>% select(MPA, version, estimate, years) %>%
    tidyr::pivot_wider(names_from = version, values_from = c(estimate, years))
  stopifnot(all(abs(same$estimate_strict - same$estimate_core80) < 1e-10),
    all(same$years_strict == same$years_core80))
  # Original densities and fitted results are byte-identical; no benchmark overwrite.
  ident <- read("identity_audit")
  stopifnot(any(ident$site_new == "SCI_SCORPION" & ident$statuses == "mpa;reference"))
}
cat("PASS: coverage selection, fixed rosters, original benchmark, species/year matching and hashes\n")
