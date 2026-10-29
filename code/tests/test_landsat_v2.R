source("code/R/40_landsat_reference_check.R")
z <- check_landsat_v2_verification(normalizePath("."))
local({
  root <- tempfile("landsat-v2-test-")
  dir.create(root)
  on.exit(unlink(root, recursive=TRUE))
  files <- c("data/harmonized/landsat_v2_verification.json",
    vapply(z$required_files, `[[`, character(1), "path"))
  for (f in files) {
    dir.create(dirname(file.path(root,f)), recursive=TRUE, showWarnings=FALSE)
    stopifnot(file.copy(f,file.path(root,f)))
  }
  # A public checkout without private source files must pass the cached gate.
  check_landsat_v2_verification(root)
  report <- file.path(root,"outputs/landsat_v2_verification/mask_array_checks.csv")
  cat("changed\n",file=report,append=TRUE)
  stopifnot(inherits(try(check_landsat_v2_verification(root),silent=TRUE),"try-error"))
  file.copy("outputs/landsat_v2_verification/mask_array_checks.csv",report,overwrite=TRUE)
  input <- file.path(root,"data/harmonized/harmonized_landsat_rr.csv")
  cat("changed\n",file=input,append=TRUE)
  stopifnot(inherits(try(check_landsat_v2_verification(root),silent=TRUE),"try-error"))
  file.copy("data/harmonized/harmonized_landsat_rr.csv",input,overwrite=TRUE)
  unlink(file.path(root,"data/harmonized/landsat_v2_verification.json"))
  stopifnot(inherits(try(check_landsat_v2_verification(root),silent=TRUE),"try-error"))
})
cat("PASS: public cached verification; changed reports, changed inputs and missing verification rejected.\n")
