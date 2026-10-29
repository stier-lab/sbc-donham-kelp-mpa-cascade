# Explicit source rebuild; routine analyses consume the pinned tracked derivative.
# Rscript code/R/export_landsat_annual_index.R /path/to/data-processing-repo
local({
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) stop("Supply the data-processing repository path.")
  data_repo <- normalizePath(args[1], mustWork = TRUE)
  analysis_repo <- normalizePath(here::here(), mustWork = TRUE)
  script <- file.path(data_repo, "code/R/export_landsat_correction.R")
  if (!file.exists(script)) stop("Corrected source exporter is missing: ", script)
  old <- getwd()
  on.exit(setwd(old))
  setwd(data_repo)
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c(shQuote(script), shQuote(analysis_repo)))
  if (status != 0L) stop("Corrected Landsat source export failed.")
})
