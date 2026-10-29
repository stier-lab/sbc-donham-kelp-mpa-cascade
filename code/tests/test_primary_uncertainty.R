suppressPackageStartupMessages(library(sandwich))
source("code/R/primary_uncertainty_utils.R")
set.seed(361012)
d <- data.frame(y = rnorm(25), t = seq_len(25))
m <- lm(y ~ t, d)
v <- uncertainty_hac(m, d$t, 3)
reference <- sandwich::NeweyWest(m, lag = 3, prewhite = FALSE, adjust = TRUE)
stopifnot(max(abs(v - reference)) < 1e-12)
nonlinear <- nls(y ~ a + b * t, d, start = list(a = 0, b = 0))
stopifnot(max(abs(uncertainty_hac(nonlinear, d$t, 3) -
  sandwich::NeweyWest(nonlinear, lag = 3, prewhite = FALSE, adjust = TRUE))) < 1e-12)
# Integer-calendar gaps must change the lag kernel, not compress time.
year <- c(1:12, 20:32)
gap <- uncertainty_hac(m, year, 3)
stopifnot(max(abs(gap - v)) > 1e-6,
  max(abs(uncertainty_hac(m, year + 1980, 3) - gap)) < 1e-12,
  min(eigen(gap, symmetric = TRUE, only.values = TRUE)$values) >= 0)
stopifnot(inherits(try(uncertainty_hac(m, rep(1, 25), 3), silent = TRUE), "try-error"))
primary <- read.csv("data/sumstats_final.csv")
annual <- uncertainty_annual(getwd(), primary)
definitions <- uncertainty_definitions(getwd())
for (i in seq_len(nrow(primary))) {
  x <- uncertainty_reconstruct(primary[i, ], annual, definitions)
  stopifnot(abs(x$estimate - primary$Mean[i]) < 1e-7,
    abs(x$iid_se - primary$SE[i]) < max(1e-7, primary$SE[i] * 1e-5))
}
cat("PASS: all 145 retained contrasts/SEs, library-equivalent HAC, calendar gaps, PSD and duplicate rejection.\n")
