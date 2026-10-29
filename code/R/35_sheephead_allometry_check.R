# Require the sourced sheephead correction before fitting or using cached data.
check_sheephead_allometry_inputs <- function(root = here::here()) {
  path <- file.path(root, "data/harmonized/sheephead_allometry.json")
  if (!file.exists(path)) stop("Hamilton sheephead allometry release marker is missing.")
  z <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  expected <- list(reference_doi = "10.1098/rspb.2014.1817", a = 0.0144, b = 3.04,
                   length_type = "TL", length_unit = "cm", weight_unit = "g",
                   count_weighted = TRUE, original_lookup_md5 = "6ccfacd18fff07208f510384e914a3b2")
  valid <- vapply(names(expected), function(k) identical(z$source[[k]], expected[[k]]), TRUE)
  if (!all(valid) || !isTRUE(z$applied) ||
      !identical(z$release, "hamilton_caselle_2015_sheephead") ||
      !identical(z$scope$programme, "PISCO") || !identical(z$scope$response, "Bio") ||
      !identical(z$scope$taxon, "Semicossyphus pulcher")) {
    stop("Wrong or unapplied sheephead allometry release.")
  }
  required <- c("harmonized_raw_responses.csv", "harmonized_response_ratios.csv")
  source(here::here("code/R/lobster_release_utils.R"),local=TRUE)
  if (!setequal(names(z$files), required)) stop("Incomplete sheephead release fingerprints.")
  for (name in required) {
    input <- file.path(root, "data/harmonized", name)
    expected_md5 <- lobster_release_hash(root,name,z$files[[name]]$after_md5)
    if (unname(tools::md5sum(input)) != expected_md5 ||
        nrow(utils::read.csv(input)) != z$files[[name]]$rows) {
      stop("Stale or changed sheephead input: ", name)
    }
  }
  cat("PASS: Hamilton sheephead formula, units, scope and input fingerprints\n")
  invisible(z)
}
