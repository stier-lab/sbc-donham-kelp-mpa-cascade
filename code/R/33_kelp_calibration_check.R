# Verify that analysis inputs are the adopted monthly-calibration release.
check_monthly_kelp_inputs <- function(root = here::here()) {
  marker <- file.path(root, "data/harmonized/kelp_monthly_calibration.json")
  if (!file.exists(marker)) stop("Monthly kelp calibration release marker is missing.")
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required for the release check.")
  z <- jsonlite::fromJSON(marker, simplifyVector = FALSE)
  calibration <- file.path(root, "data/external/sbc11215/table5.csv")
  if (!identical(z$package, "knb-lter-sbc.112.15") || z$model != 1 ||
      !isTRUE(z$applied) || !identical(z$mass_basis, "dry foliar") ||
      !identical(z$calibration_md5, "50eda5d3664bf70780b3bead9674fac1") ||
      unname(tools::md5sum(calibration)) != z$calibration_md5) {
    stop("Wrong or unapplied monthly kelp calibration release.")
  }
  followup_path <- file.path(root, "data/harmonized/sheephead_allometry.json")
  followup <- if (file.exists(followup_path)) jsonlite::fromJSON(followup_path, simplifyVector = FALSE) else NULL
  source(here::here("code/R/lobster_release_utils.R"),local=TRUE)
  for (name in names(z$files)) {
    path <- file.path(root, "data/harmonized", name)
    expected_md5 <- z$files[[name]]$after_md5
    if (!is.null(followup)) {
      if (!identical(followup$release, "hamilton_caselle_2015_sheephead") ||
          !isTRUE(followup$applied) ||
          !identical(followup$files[[name]]$before_md5, expected_md5)) {
        stop("Sheephead release does not descend from the checked monthly kelp release.")
      }
      expected_md5 <- followup$files[[name]]$after_md5
    }
    expected_md5 <- lobster_release_hash(root,name,expected_md5)
    if (unname(tools::md5sum(path)) != expected_md5 ||
        nrow(utils::read.csv(path)) != z$files[[name]]$rows) {
      stop("Harmonized input differs from the checked monthly release: ", name)
    }
  }
  cat("PASS: monthly kelp calibration release chain and harmonized input fingerprints\n")
  invisible(z)
}

check_monthly_kelp_inputs()
