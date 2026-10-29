# Exercise the actual parsed Naples eligibility guard without sourcing the pipeline.
source("code/R/00c_analysis_constants.R")
ast <- as.list(parse("code/R/08_effect_sizes.R"))
env <- new.env(parent = globalenv())
for (expr in ast) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) &&
      as.character(expr[[2]]) %in% c("add_time_columns", "test_linear_before")) eval(expr, env)
}
raw <- read.csv("data/harmonized/harmonized_response_ratios.csv")
d <- subset(raw, CA_MPA_Name_Short == "Naples SMCA" & source == "LTER" &
  y == "Semicossyphus pulcher" & resp == "Bio")
env$LTER.nap.SPUL.bio <- env$add_time_columns(d)
env$naples_biomass_before_p <- env$test_linear_before(env$LTER.nap.SPUL.bio)$pSlope
env$naples_biomass_weights <- data.frame(Step = 17.2, Linear = 51.0, Asymptotic = 19.0, Sigmoid = 12.8)
guard <- Filter(function(expr) is.call(expr) && identical(expr[[1]], as.name("stopifnot")) &&
  grepl("nrow\\(LTER.nap.SPUL.bio\\) == 23L", paste(deparse(expr), collapse = " ")), ast)
stopifnot(length(guard) == 1L)
eval(guard[[1]], env)
stopifnot(abs(env$naples_biomass_before_p - 0.0988993259) < 1e-8,
  sum(d$BA == "Before") == 12L, sum(d$BA == "After") == 11L)
reject <- function(change) {
  trial <- list2env(as.list(env), parent = globalenv())
  change(trial)
  stopifnot(inherits(try(eval(guard[[1]], trial), silent = TRUE), "try-error"))
}
reject(function(e) e$naples_biomass_before_p <- .01)
reject(function(e) e$naples_biomass_before_p <- NA_real_)
reject(function(e) e$naples_biomass_weights$Step <- 99)
reject(function(e) e$LTER.nap.SPUL.bio$year[2] <- e$LTER.nap.SPUL.bio$year[1])
reject(function(e) e$LTER.nap.SPUL.bio$lnDiff[1] <- NA_real_)
reject(function(e) e$LTER.nap.SPUL.bio <- e$LTER.nap.SPUL.bio[-1, ])
cat("PASS: Naples eligibility guard accepts current history and rejects invalid screens, model routes, duplicates, missing values and changed coverage.\n")
