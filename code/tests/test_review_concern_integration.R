# Focused bridge integration only; never sources run_all.R or writes repo artifacts.
local({
  run_tests <- function() {
    if (!requireNamespace("here", quietly = TRUE)) stop("Package here is required.")
    repo <- here::here()
    modules <- c("27_design_geography_sensitivity.R", "28_historical_heatwave_bacips.R",
                 "29_pair_history_atlas.R", "30_review_concern_checks.R")
    inputs <- c("data/sumstats_final.csv", "tables/table_02_meta_analysis.csv",
      file.path("data/harmonized", c("harmonized_response_ratios.csv", "harmonized_landsat_rr.csv",
        "harmonized_landsat_annual_index.csv", "harmonized_site_metadata.csv", "reference_sites.csv",
        "reference_primary_keys.csv", "mapping_provenance.csv", "checksums.md5")))
    codes <- file.path("code/R", c(modules, "00c_analysis_constants.R", "00b_color_palette.R", "review_concern_utils.R", "taxon_display.R"))
    fixture <- tempfile("review_bridge_integration_")
    for (dir in c("code/R", "data/harmonized", "data/cache", "tables", "plots", "docs"))
      dir.create(file.path(fixture, dir), recursive = TRUE)
    stopifnot(all(file.copy(file.path(repo, c(inputs, codes)), file.path(fixture, c(inputs, codes)))))
    file.create(file.path(fixture, ".here"))
    original_wd <- getwd()
    on.exit({
      setwd(original_wd)
      here::i_am("code/tests/test_review_concern_integration.R")
    }, add = TRUE)
    setwd(fixture)
    for (module in modules) {
      cat("Standalone bridge:", module, "\n")
      log <- file.path(fixture, paste0(module, ".log"))
      status <- system2(file.path(R.home("bin"), "Rscript"),
        c("--vanilla", shQuote(file.path(fixture, "code/R", module))), stdout = log, stderr = log)
      if (status != 0L) stop(module, " failed:\n", paste(readLines(log), collapse = "\n"))
    }
    files <- list.files(file.path(fixture, "tables"), pattern = "^table_s_(design|historical_heatwave|pair_history_atlas).*csv$", full.names = TRUE)
    standalone <- stats::setNames(lapply(files, utils::read.csv), basename(files))
    for (name in names(standalone)) {
      if (!isTRUE(all.equal(standalone[[name]], utils::read.csv(file.path(repo, "tables", name)),
                           check.attributes = FALSE, tolerance = 1e-9)))
        stop("Standalone bridge changed an existing table: ", name)
    }

    # Reproduce package masking and inherited pipeline objects without running the pipeline.
    suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(plyr); library(MASS) })
    shared <- new.env(parent = globalenv())
    for (name in c("input", "primary", "annual", "metadata", "pairs", "summary", "coverage", "estimates",
      "diagnostics", "reference_sites", "reference_groups", "outputs", "p", "de", "dc", "gs"))
      assign(name, paste0("pipeline sentinel: ", name), envir = shared)
    for (name in c("select", "filter", "summarise", "mutate", "group_by", "ggplot", "geom_point", "%>%"))
      assign(name, function(...) stop("Unqualified package lookup reached inherited pipeline mask"), envir = shared)
    original <- as.list(shared, all.names = TRUE)
    here::i_am("code/R/27_design_geography_sensitivity.R")
    guard <- function(paths) {
      resolved <- normalizePath(paths, mustWork = FALSE)
      if (any(!startsWith(resolved, paste0(normalizePath(fixture), "/"))))
        stop("Unexpected input/output outside the tracked fixture: ", paste(paths, collapse = ", "))
    }
    guard_name <- ".review_bridge_test_guard"
    had_guard <- exists(guard_name, envir = globalenv(), inherits = FALSE)
    previous_guard <- if (had_guard) get(guard_name, envir = globalenv()) else NULL
    assign(guard_name, guard, envir = globalenv())
    on.exit({
      if (had_guard) assign(guard_name, previous_guard, envir = globalenv()) else rm(list = guard_name, envir = globalenv())
    }, add = TRUE)
    trace("read.csv", where = asNamespace("utils"), print = FALSE,
      tracer = quote(.GlobalEnv$.review_bridge_test_guard(file)))
    trace("write.csv", where = asNamespace("utils"), print = FALSE,
      tracer = quote({
        arguments <- list(...)
        .GlobalEnv$.review_bridge_test_guard(if (is.null(arguments$file)) arguments[[2L]] else arguments$file)
      }))
    trace("md5sum", where = asNamespace("tools"), print = FALSE,
      tracer = quote(.GlobalEnv$.review_bridge_test_guard(files)))
    on.exit({
      untrace("read.csv", where = asNamespace("utils"))
      untrace("write.csv", where = asNamespace("utils"))
      untrace("md5sum", where = asNamespace("tools"))
    }, add = TRUE)
    for (module in modules) {
      cat("Shared-environment bridge:", module, "\n")
      base::source(file.path(fixture, "code/R", module), local = shared)
      stopifnot(identical(as.list(shared, all.names = TRUE), original))
    }
    for (name in names(standalone)) {
      current <- utils::read.csv(file.path(fixture, "tables", name))
      if (!isTRUE(all.equal(standalone[[name]], current, check.attributes = FALSE, tolerance = 1e-9)))
        stop("Shared-environment result changed: ", name)
    }

    expect_error <- function(expression, pattern) {
      error <- tryCatch({ force(expression); NULL }, error = identity)
      if (!inherits(error, "error") || !grepl(pattern, conditionMessage(error)))
        stop("Expected rejection matching ", pattern, "; got ", if (is.null(error)) "success" else conditionMessage(error))
    }
    corrupt_output <- function(name, column, pattern) {
      path <- file.path(fixture, "tables", name)
      original_file <- tempfile(tmpdir = fixture)
      stopifnot(file.copy(path, original_file))
      on.exit(file.copy(original_file, path, overwrite = TRUE), add = TRUE)
      data <- utils::read.csv(path)
      i <- which(is.finite(data[[column]]) & if (column == "log_lik") data$selected else TRUE)[1]
      data[[column]][i] <- data[[column]][i] + 1
      utils::write.csv(data, path, row.names = FALSE)
      expect_error(base::source(file.path(fixture, "code/R", modules[4]), local = shared), pattern)
    }
    corrupt_output("table_s_design_group_summary.csv", "se", "Group summary se")
    corrupt_output("table_s_design_leave_reference_group_out.csv", "estimate", "Influence estimate")
    corrupt_output("table_s_design_pisco_common_calendar_coverage.csv", "common_calendar_start", "Common-calendar coverage")
    corrupt_output("table_s_design_gls_fit_attempts.csv", "log_lik", "likelihood")
    checksum_path <- file.path(fixture, "data/harmonized/checksums.md5")
    checksums <- readLines(checksum_path)
    changed <- checksums
    changed[grepl(" mapping_provenance.csv$", changed)] <- paste0(strrep("0", 32), " mapping_provenance.csv")
    writeLines(changed, checksum_path)
    expect_error(base::source(file.path(fixture, "code/R", modules[4]), local = shared), "checksum manifest disagrees")
    writeLines(checksums, checksum_path)
    base::source(file.path(fixture, "code/R", modules[4]), local = shared)
    stopifnot(identical(as.list(shared, all.names = TRUE), original))
    cat("Bridge integration tests passed: standalone/shared results match, globals preserved, raw-free reads and corruption guards verified.\n")
  }
  run_tests()
})
