suppressPackageStartupMessages({library(dplyr); library(tidyr)})
source("code/R/coauthor_sensitivity_utils.R")
source("code/R/primary_uncertainty_utils.R")
source("code/R/predator_size_refit_utils.R")
d <- size_refit_annual(here::here())
stopifnot(nrow(d) > 0, all(d$time > 0), all(d$source == "PISCO"),
  setequal(unique(d$size), c("all", "large")), !anyDuplicated(d[c("MPA", "year", "Taxa", "Resp", "size")]))
y <- c(2001, 2002, 2005, 2006, 2010, 2011)
v <- .2 * (y - min(y)) + c(.1, -.1, .05, -.05, -.1, .1)
linear <- size_refit_history(v, y, "Linear")
fit <- lm(v ~ I(y - min(y)))
stopifnot(abs(linear$estimate - coef(fit)[2] * 11) < 1e-12,
  linear$variance >= (11 * summary(fit)$coef[2, 2])^2 - 1e-12,
  linear$extrapolated_11, size_refit_history(v, y)$model == "Mean")
stopifnot(inherits(try(size_refit_history(v, rep(y[1], 6)), silent = TRUE), "try-error"))
cat("PASS: paired size annual keys, zero/offset conventions, calendar time and conservative HAC contrasts.\n")

status <- read.csv("tables/table_s_size_refit_status.csv")
stopifnot(all(status$status == "completed"), !any(status$primary_changed))
inputs <- read.csv("tables/table_s_size_refit_inputs.csv")
stopifnot(identical(unname(tools::md5sum(inputs$path)), as.character(inputs$md5)))
annual <- read.csv("tables/table_s_size_refit_temporal_annual.csv")
for (scope in unique(annual$scope)) for (response in c("Bio", "Den")) {
  x <- filter(annual, scope == .env$scope, Resp == response)
  a <- filter(x, size == "all") %>% arrange(history, year)
  b <- filter(x, size == "large") %>% arrange(history, year)
  rownames(a) <- rownames(b) <- NULL
  stopifnot(identical(a[c("history", "year", "time")], b[c("history", "year", "time")]))
}
candidates <- read.csv("tables/table_s_size_refit_temporal_candidates.csv")
stopifnot(all(candidates$converged[candidates$selected]), !any(candidates$singular[candidates$selected]),
  !any(grepl("negative eigenvalue|failed", candidates$warnings[candidates$selected], ignore.case = TRUE)))
slopes <- read.csv("tables/table_s_size_refit_temporal_slopes.csv")
stopifnot(all(is.finite(slopes$estimate)), all(slopes$ci_low < slopes$ci_high),
  all(c("calendar_AR1", "MPA_cluster_CR2", "independent_errors") %in% slopes$model))
predictions <- read.csv("tables/table_s_size_refit_temporal_predictions.csv")
stopifnot(setequal(unique(predictions$model), c("calendar_AR1", "independent_errors")))
paths <- read.csv("tables/table_s_size_refit_pathways.csv")
stopifnot(all(is.na(paths$slope[paths$status != "fitted"])),
  all(paths$n_mpa[paths$status == "fitted"] >= 4),
  all(grepl("Conditional", paths$inference)))
cat("PASS: fresh input hashes, exact matched cohorts, valid covariance hierarchy, explicit model-labelled intervals and guarded small-MPA pathways.\n")
