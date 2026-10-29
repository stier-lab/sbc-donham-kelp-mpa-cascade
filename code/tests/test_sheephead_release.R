source("code/R/35_sheephead_allometry_check.R")
local({
  root <- tempfile("sheephead-release-")
  input <- file.path(root, "data/harmonized")
  dir.create(input, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  files <- c("harmonized_raw_responses.csv", "harmonized_response_ratios.csv")
  for (name in files) write.csv(data.frame(value = 1), file.path(input, name), row.names = FALSE)
  z <- list(release = "hamilton_caselle_2015_sheephead", applied = TRUE,
    source = list(reference_doi = "10.1098/rspb.2014.1817", a = 0.0144, b = 3.04,
      length_type = "TL", length_unit = "cm", weight_unit = "g", count_weighted = TRUE,
      original_lookup_md5 = "6ccfacd18fff07208f510384e914a3b2"),
    scope = list(programme = "PISCO", taxon = "Semicossyphus pulcher", response = "Bio"),
    files = setNames(lapply(files, function(name) list(rows = 1,
      after_md5 = unname(tools::md5sum(file.path(input, name))))), files))
  marker <- file.path(input, "sheephead_allometry.json")
  write_marker <- function(x) jsonlite::write_json(x, marker, auto_unbox = TRUE, digits = NA)
  expect_error <- function(expr, pattern) {
    error <- tryCatch({ force(expr); NULL }, error = identity)
    stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
  }
  expect_error(check_sheephead_allometry_inputs(root), "missing")
  write_marker(z)
  check_sheephead_allometry_inputs(root)
  for (field in c("a", "b", "length_type", "length_unit", "weight_unit", "reference_doi")) {
    wrong <- z
    wrong$source[[field]] <- "wrong"
    write_marker(wrong)
    expect_error(check_sheephead_allometry_inputs(root), "Wrong")
  }
  wrong <- z
  wrong$scope$response <- "Den"
  write_marker(wrong)
  expect_error(check_sheephead_allometry_inputs(root), "Wrong")
  wrong <- z
  wrong$applied <- FALSE
  write_marker(wrong)
  expect_error(check_sheephead_allometry_inputs(root), "Wrong")
  write_marker(z)
  # The correction must link to the preserved monthly-calibration release.
  definition <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
    identical(x[[2]], as.name("check_monthly_kelp_inputs")),
    as.list(parse("code/R/33_kelp_calibration_check.R")))
  stopifnot(length(definition) == 1L)
  eval(definition[[1]])
  calibration_dir <- file.path(root, "data/external/sbc11215")
  dir.create(calibration_dir, recursive = TRUE)
  stopifnot(file.copy("data/external/sbc11215/table5.csv", calibration_dir))
  parent <- list(package = "knb-lter-sbc.112.15", model = 1, applied = TRUE,
    mass_basis = "dry foliar", calibration_md5 = "50eda5d3664bf70780b3bead9674fac1",
    files = setNames(lapply(files, function(name) list(rows = 1, after_md5 = "parent-hash")), files))
  jsonlite::write_json(parent, file.path(input, "kelp_monthly_calibration.json"), auto_unbox = TRUE)
  for (name in files) z$files[[name]]$before_md5 <- "parent-hash"
  write_marker(z)
  check_monthly_kelp_inputs(root)
  wrong <- z
  wrong$files[[files[1]]]$before_md5 <- "wrong-parent"
  write_marker(wrong)
  expect_error(check_monthly_kelp_inputs(root), "does not descend")
  write_marker(z)
  write.csv(data.frame(value = 2), file.path(input, files[1]), row.names = FALSE)
  expect_error(check_sheephead_allometry_inputs(root), "Stale")
  cat("PASS: missing/unapplied release, coefficients, DOI, units, scope and stale fingerprints rejected.\n")
})
