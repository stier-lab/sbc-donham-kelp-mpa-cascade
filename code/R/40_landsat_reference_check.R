# Pin the corrected satellite handoff without requiring private source files.
check_landsat_reference_inputs <- function(root = here::here()) {
  path <- file.path(root, "data/harmonized/landsat_reference_correction.json")
  if (!file.exists(path)) stop("Corrected Landsat reference release marker is missing.")
  z <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  expected <- c("Skunk Pt SMR", "Santa Barbara Island SMR", "South La Jolla SMR",
                "Cabrillo SMR", "Point Dume SMCA", "Point Dume SMR")
  if (!isTRUE(z$applied) || !identical(z$release, "landsat_empty_reference_slots_2026_10_16") ||
      !setequal(unlist(z$affected_mpas), expected) || z$source_rows - z$retained_source_rows != 6L ||
      !isTRUE(z$diver_inputs_unchanged)) stop("Invalid Landsat reference correction release.")
  for (name in c("harmonized_landsat_rr.csv", "harmonized_landsat_annual_index.csv")) {
    input <- file.path(root, "data/harmonized", name)
    if (unname(tools::md5sum(input)) != z$files[[name]]$md5 ||
        nrow(read.csv(input)) != z$files[[name]]$rows) stop("Changed Landsat input: ", name)
  }
  cat("PASS: corrected Landsat reference slots and input fingerprints\n")
  check_landsat_v2_verification(root)
  invisible(z)
}

check_landsat_v2_verification <- function(root = here::here()) {
  path <- file.path(root, "data/harmonized/landsat_v2_verification.json")
  if (!file.exists(path)) stop("Landsat v2 verification manifest is missing.")
  z <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  if (!isTRUE(z$verified) || z$release != "landsat_v2_masks_2026_10_21" ||
      z$mpa_arrays != 24L || z$reference_arrays != 34L || z$empty_groups != 14L ||
      z$input_rows != 60L || z$confirmed_removed_slots != 6L || z$max_abs_difference_kg != 0)
    stop("Invalid Landsat v2 verification record.")
  required <- c("data/harmonized/harmonized_landsat_rr.csv",
    "data/harmonized/harmonized_landsat_annual_index.csv",
    "outputs/landsat_v2_verification/mask_array_checks.csv",
    "outputs/landsat_v2_verification/input_checks.csv")
  if (!setequal(vapply(z$required_files, `[[`, character(1), "path"), required))
    stop("Incomplete Landsat v2 verification fingerprints.")
  for (f in z$required_files) {
    p <- file.path(root,f$path)
    if (!file.exists(p) || unname(tools::md5sum(p)) != f$md5)
      stop("Changed Landsat v2 verified input/report: ",f$path)
  }
  # Public checkouts need no private inputs; present local sources must still match.
  for (f in z$local_source_files) {
    p <- file.path(root,f$path)
    if (file.exists(p) && unname(tools::md5sum(p)) != f$md5)
      stop("Changed local Landsat source; rerun verification: ",f$path)
  }
  cat("PASS: Landsat v2 masks verified against 24 MPA and 34 reference arrays\n")
  invisible(z)
}
