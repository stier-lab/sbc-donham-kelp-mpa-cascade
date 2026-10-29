# Standalone/shared-environment regression checks; all writes stay in /tmp.
local({
  run_tests <- function() {
    script <- "32_pbacips_shape_display.R"
    paths <- unlist(lapply(sys.frames(), function(x) x$ofile))
    paths <- c(paths, sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)))
    paths <- paths[basename(paths) == "test_pbacips_shape_display.R"]
    repo <- if (length(paths)) dirname(dirname(dirname(normalizePath(tail(paths, 1L))))) else getwd()
    fixture <- tempfile("pbacips_shape_tests_", tmpdir = "/tmp")
    for (directory in c("code/R", "data/harmonized", "data/cache", "tables", "plots", "docs")) {
      dir.create(file.path(fixture, directory), recursive = TRUE, showWarnings = FALSE)
    }
    code_names <- c(script, "02_pBACIPS_function.R", "08_effect_sizes.R", "00c_analysis_constants.R", "taxon_display.R")
    stopifnot(all(file.copy(file.path(repo, "code/R", code_names), file.path(fixture, "code/R", code_names))))
    years <- c(2008:2014, 2016, 2018:2023)
    ba <- ifelse(years <= 2012, "Before", "After")
    index <- c(rep(0, sum(ba == "Before")), seq_len(sum(ba == "After")) - 1L)
    supplied <- pmax(years - 2012, 0)
    noise <- rep(c(-0.01, 0.01, 0.02, -0.02), length.out = length(years))
    hill <- function(t) 0.3 + 1.7 * ((t + 0.01) / 5)^2 / (1 + ((t + 0.01) / 5)^2)
    definitions <- data.frame(
      MPA = c("Index SMR", "Step SMR", "Satellite SMR", "Campus Point SMCA", "Scorpion SMR", "CI SMCA"),
      Taxa = c("S. pulcher", "M. franciscanus", "M. pyrifera", "P. interruptus", "M. pyrifera", "S. pulcher"),
      y = c("Semicossyphus pulcher", "Mesocentrotus franciscanus", "Macrocystis pyrifera", "Panulirus interruptus", "Macrocystis pyrifera", "Semicossyphus pulcher"),
      Source = c("LTER", "KFM", "Landsat", "LTER", "KFM", "PISCO"),
      Resp = c("Bio", "Den", "Bio", "Bio", "Bio", "Den"),
      Model = c("Linear", "Step", "Linear", "Linear", "Sigmoid", "Mean"),
      AnalysisType = c(rep("pBACIPS", 5), "CI"), stringsAsFactors = FALSE)
    values <- list(0.2 + 0.1 * index + noise, ifelse(ba == "After", 1.2, 0.3) + noise,
      0.5 + 0.08 * supplied + noise, 0.3 + 0.06 * supplied + noise, hill(index), 0.4 + noise)
    annual <- do.call(rbind, lapply(seq_len(nrow(definitions)), function(i) {
      z <- definitions[i, ]
      data.frame(CA_MPA_Name_Short = z$MPA, year = years, y = z$y, source = z$Source,
        resp = z$Resp, lnDiff = values[[i]], BA = ba, time = supplied, stringsAsFactors = FALSE)
    }))
    slope_effect <- function(y, x) sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2) * 11
    effects <- c(slope_effect(values[[1]], index),
      mean(values[[2]][ba == "After"]) - mean(values[[2]][ba == "Before"]),
      slope_effect(values[[3]], supplied), slope_effect(values[[4]], supplied), hill(11) - hill(0), mean(values[[6]]))
    primary <- definitions[c("MPA", "Taxa", "Source", "Resp", "Model", "AnalysisType")]
    primary$Mean <- effects
    primary$SE <- 0.1
    primary$N <- length(years)
    primary$Primary <- "Y"
    primary$MPA_Start <- 2012L
    write_primary <- function(p = primary) {
      utils::write.csv(p, file.path(fixture, "data/sumstats_final.csv"), row.names = FALSE)
      utils::write.csv(p, file.path(fixture, "tables/table_model_selection.csv"), row.names = FALSE)
    }
    write_annual <- function(a = annual) {
      utils::write.csv(a[a$source != "Landsat", ], file.path(fixture, "data/harmonized/harmonized_response_ratios.csv"), row.names = FALSE)
      utils::write.csv(a[a$source == "Landsat", ], file.path(fixture, "data/harmonized/harmonized_landsat_rr.csv"), row.names = FALSE)
    }
    write_primary()
    write_annual()
    utils::write.csv(data.frame(CA_MPA_Name_Short = definitions$MPA, MPA_Start = 2012L),
      file.path(fixture, "data/harmonized/harmonized_site_metadata.csv"), row.names = FALSE)
    run <- function(flag = NULL, error = NULL) {
      log <- tempfile(tmpdir = fixture, fileext = ".log")
      status <- system2(file.path(R.home("bin"), "Rscript"),
        c("--vanilla", shQuote(file.path(fixture, "code/R", script)), flag), stdout = log, stderr = log)
      output <- paste(readLines(log, warn = FALSE), collapse = "\n")
      if (is.null(error) && status != 0L) stop("Unexpected script failure:\n", output)
      if (!is.null(error) && (status == 0L || !grepl(error, output))) stop("Expected ", error, "; got:\n", output)
      invisible(output)
    }
    run("--validate-only")
    stopifnot(!file.exists(file.path(fixture, "data/cache/pbacips_shape_display.rds")))
    run(error = "cache absent")
    run("--reconstruct")
    audit <- utils::read.csv(file.path(fixture, "tables/table_s_pbacips_shape_fit_audit.csv"))
    membership <- utils::read.csv(file.path(fixture, "tables/table_s_pbacips_shape_membership.csv"))
    observations <- utils::read.csv(file.path(fixture, "tables/table_s_pbacips_shape_observations.csv"))
    predictions <- utils::read.csv(file.path(fixture, "tables/table_s_pbacips_shape_predictions.csv"))
    stopifnot(nrow(audit) == 5L, all(audit$status == "verified"), all(abs(audit$difference) < 1e-7),
      nrow(membership) == 6L, sum(membership$display_scope == "not_displayed_ci") == 1L,
      membership$nonlinear_panel[membership$Model == "Sigmoid"] == "A",
      audit$clock[audit$MPA == "Satellite SMR"] == "supplied_time",
      audit$clock[audit$MPA == "Campus Point SMCA"] == "supplied_time",
      audit$clock[audit$MPA == "Index SMR"] == "observation_index",
      audit$first_after_model_time[audit$MPA == "Index SMR"] == 0,
      audit$first_after_model_time[audit$MPA == "Satellite SMR"] == 1,
      !audit$horizon_within_observed_clock[audit$MPA == "Index SMR"],
      is.na(audit$horizon[audit$Model == "Step"]),
      all(audit$missing_years == "2015;2017"),
      all(abs(predictions$time.model[predictions$MPA == "Satellite SMR"] -
        pmax(predictions$year[predictions$MPA == "Satellite SMR"] - 2012, 0)) < 1e-10),
      observations$time.model[observations$MPA == "Index SMR" & observations$year == 2016] == 2,
      observations$time.model[observations$MPA == "Satellite SMR" & observations$year == 2016] == 4)
    step_predictions <- predictions[predictions$MPA == "Step SMR", ]
    stopifnot(length(unique(step_predictions$fitted_lnRR[step_predictions$year <= 2012.5])) == 1L,
      length(unique(step_predictions$fitted_lnRR[step_predictions$year > 2012.5])) == 1L,
      abs(unique(step_predictions$fitted_lnRR[step_predictions$year > 2012.5]) -
        unique(step_predictions$fitted_lnRR[step_predictions$year <= 2012.5]) - effects[2]) < 1e-7)
    run("--check-cache")

    # Source with inherited masking/sentinels. Any fitting during cached render is fatal.
    shared <- new.env(parent = globalenv())
    for (name in c("root", "primary", "selected", "predictions", "observations", "MODEL_FIT_LOG")) {
      assign(name, paste("sentinel", name), envir = shared)
    }
    for (name in c("filter", "select", "ggplot", "lm", "nls")) {
      assign(name, function(...) stop("Inherited mask called"), envir = shared)
    }
    previous <- as.list(shared, all.names = TRUE)
    old_mode <- Sys.getenv("MPA_SHAPE_DISPLAY_MODE", unset = NA_character_)
    on.exit(if (is.na(old_mode)) Sys.unsetenv("MPA_SHAPE_DISPLAY_MODE") else
      Sys.setenv(MPA_SHAPE_DISPLAY_MODE = old_mode), add = TRUE)
    Sys.setenv(MPA_SHAPE_DISPLAY_MODE = "render")
    set.seed(32)
    old_rng <- .Random.seed
    old_options <- options()
    old_search <- search()
    Sys.setenv(MPA_SHAPE_DISPLAY_MODE = "reconstruct")
    base::source(file.path(fixture, "code/R", script), local = shared)
    stopifnot(identical(as.list(shared, all.names = TRUE), previous), identical(old_rng, .Random.seed),
      identical(old_options, options()), identical(old_search, search()))
    Sys.setenv(MPA_SHAPE_DISPLAY_MODE = "render")
    trace("lm", where = asNamespace("stats"), print = FALSE, tracer = quote(stop("Unexpected cached-render fit")))
    trace("nls", where = asNamespace("stats"), print = FALSE, tracer = quote(stop("Unexpected cached-render fit")))
    tryCatch(base::source(file.path(fixture, "code/R", script), local = shared), finally = {
      untrace("lm", where = asNamespace("stats"))
      untrace("nls", where = asNamespace("stats"))
    })
    stopifnot(identical(as.list(shared, all.names = TRUE), previous), identical(old_rng, .Random.seed),
      identical(old_options, options()), identical(old_search, search()))

    cache_path <- file.path(fixture, "data/cache/pbacips_shape_display.rds")
    cache <- readRDS(cache_path)
    bad <- cache
    bad$predictions$fitted_lnRR[1] <- NA_real_
    saveRDS(bad, cache_path)
    run("--check-cache", "Invalid cached prediction cohort")
    saveRDS(cache, cache_path)
    bad <- cache
    bad$predictions$time.model[1] <- bad$predictions$time.model[1] + 0.1
    saveRDS(bad, cache_path)
    run("--check-cache", "native clock/model labels")
    saveRDS(cache, cache_path)
    bad <- cache
    bad$audit$reconstructed_effect[1] <- bad$audit$reconstructed_effect[1] + 1
    saveRDS(bad, cache_path)
    run("--check-cache", "Cached fit audit")
    saveRDS(cache, cache_path)
    changed <- annual
    changed$lnDiff[1] <- changed$lnDiff[1] + 0.1
    write_annual(changed)
    run("--check-cache", "Stale shape cache")
    write_annual()
    write_annual(rbind(annual, annual[1, ]))
    run("--validate-only", "Duplicate source-specific annual keys")
    write_annual()
    changed <- primary
    changed$N[1] <- changed$N[1] + 1
    write_primary(changed)
    run("--validate-only", "Annual input/primary N")
    write_primary()

    # A reversed Step effect must fail, not be hidden by flipping a curve or dropping a key.
    before_cache <- unname(tools::md5sum(cache_path))
    changed <- primary
    changed$Mean[changed$Model == "Step"] <- -changed$Mean[changed$Model == "Step"]
    write_primary(changed)
    run("--reconstruct", "Shape reconstruction failed")
    failed <- utils::read.csv(file.path(fixture, "tables/table_s_pbacips_shape_fit_audit.csv"))
    stopifnot(sum(failed$status == "failed") == 1L, failed$Model[failed$status == "failed"] == "Step",
      abs(failed$reconstructed_effect[failed$status == "failed"] - effects[2]) < 1e-7,
      identical(unname(tools::md5sum(cache_path)), before_cache))
    write_primary()
    run("--check-cache")
    for (runner in c("run_all.R", "run_review_concerns.R", "run_figures_only.R")) {
      text <- paste(readLines(file.path(repo, "code/R", runner), warn = FALSE), collapse = "\n")
      stopifnot(grepl(script, text, fixed = TRUE), grepl("MPA_SHAPE_DISPLAY_MODE", text, fixed = TRUE))
      expected_mode <- if (runner == "run_figures_only.R") "render" else "reconstruct"
      stopifnot(grepl(paste0('MPA_SHAPE_DISPLAY_MODE = "', expected_mode, '"'), text, fixed = TRUE))
      expressions <- Filter(function(x) grepl(script, paste(deparse(x), collapse = "\n"), fixed = TRUE),
        as.list(parse(text = text)))
      stopifnot(length(expressions) == 1L)
      wrapper <- new.env(parent = globalenv())
      wrapper$RENDER_FIGURES <- "all"
      wrapper$source <- wrapper$source_module <- function(...) {
        wrapper$seen_mode <- Sys.getenv("MPA_SHAPE_DISPLAY_MODE")
        if (wrapper$fail) stop("Synthetic runner source failure")
        invisible(NULL)
      }
      for (setting in c("prior-setting", NA_character_)) {
        for (fail in c(FALSE, TRUE)) {
          wrapper$fail <- fail
          wrapper$seen_mode <- NULL
          if (is.na(setting)) Sys.unsetenv("MPA_SHAPE_DISPLAY_MODE") else
            Sys.setenv(MPA_SHAPE_DISPLAY_MODE = setting)
          result <- tryCatch(eval(expressions[[1]], envir = wrapper), error = identity)
          stopifnot(identical(wrapper$seen_mode, expected_mode),
            identical(Sys.getenv("MPA_SHAPE_DISPLAY_MODE", unset = NA_character_), setting),
            identical(inherits(result, "error"), fail))
        }
      }
    }
    cat("PASS: exact clocks, complete membership/panels, contrasts, cache-only rendering, shared-state preservation, stale/duplicate/N/corrupt-cache guards, reversed-Step rejection and runner modes. Fixture:", fixture, "\n")
  }
  run_tests()
})
