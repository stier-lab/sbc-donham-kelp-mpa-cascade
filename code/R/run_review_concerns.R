# Standalone regeneration from tracked inputs and the current primary effects.
if (!requireNamespace("here", quietly = TRUE)) stop("Package here is required.")
source(here::here("code", "R", "50_before_after_cohorts.R"))
run_before_after_cohorts()
for (module in c("27_design_geography_sensitivity.R", "28_historical_heatwave_bacips.R",
                 "29_pair_history_atlas.R", "30_review_concern_checks.R",
                 "31_primary_pathway_display.R")) {
  source(here::here("code", "R", module))
}

# This regeneration runner explicitly reconstructs shapes from current primary inputs.
(function() {
  previous <- Sys.getenv("MPA_SHAPE_DISPLAY_MODE", unset = NA_character_)
  on.exit(if (is.na(previous)) Sys.unsetenv("MPA_SHAPE_DISPLAY_MODE") else
    Sys.setenv(MPA_SHAPE_DISPLAY_MODE = previous), add = TRUE)
  Sys.setenv(MPA_SHAPE_DISPLAY_MODE = "reconstruct")
  source(here::here("code", "R", "32_pbacips_shape_display.R"))
})()

# Fixed-cohort eligibility and clock sensitivity, preserving current primary effects.
source(here::here("code", "R", "52_routing_clock_sensitivity.R"))
