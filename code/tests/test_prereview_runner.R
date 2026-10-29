# Fixture-only orchestration tests: never source the full runner or fit analyses.
local({
  run_tests <- function() {
    if (!requireNamespace("here", quietly = TRUE)) stop("Package here is required.")
    here::i_am("code/tests/test_prereview_runner.R")
    root <- here::here()
    runner <- file.path(root, "code/R/run_all.R")
    expressions <- as.list(parse(runner, keep.source = FALSE))
    definition <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
      identical(x[[2]], as.name("source_module")), expressions)
    stopifnot(length(definition) == 1L)

    # Check actual top-level calls without evaluating any pipeline expressions.
    calls <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("source_module")), expressions)
    filenames <- vapply(calls, function(x) x[[2]], character(1))
    positions <- match(c("09_meta_analysis.R", "25_prereview_sensitivities.R",
                         "30_review_concern_checks.R"), filenames)
    stopifnot(!anyNA(positions), positions[1] < positions[2], positions[2] < positions[3],
      sum(filenames == "25_prereview_sensitivities.R") == 1L,
      identical(calls[[positions[2]]]$isolated, TRUE),
      identical(calls[[positions[2]]][[3]], "25"))

    protected <- c(list.files(file.path(root, "tables"), full.names = TRUE),
                   list.files(file.path(root, "plots"), full.names = TRUE, recursive = TRUE),
                   file.path(root, "data/sumstats_final.csv"))
    protected <- protected[file.info(protected)$isdir %in% FALSE]
    fingerprint <- tools::md5sum(protected)
    fixture <- tempfile("prereview runner with spaces ")
    dir.create(file.path(fixture, "code/R"), recursive = TRUE)
    stopifnot(file.copy(runner, file.path(fixture, "code/R/run_all.R")))
    file.create(file.path(fixture, ".here"))
    original_wd <- getwd()
    original_libraries <- .libPaths()
    dir.create(file.path(fixture, "fixture-library"))
    .libPaths(c(file.path(fixture, "fixture-library"), original_libraries))
    had_rng <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    original_rng <- if (had_rng) get(".Random.seed", envir = globalenv()) else NULL
    original_kind <- RNGkind()
    on.exit({
      setwd(original_wd)
      here::i_am("code/tests/test_prereview_runner.R")
      .libPaths(original_libraries)
      do.call(RNGkind, as.list(original_kind))
      if (had_rng) assign(".Random.seed", original_rng, envir = globalenv())
      else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
        rm(list = ".Random.seed", envir = globalenv())
      unlink(fixture, recursive = TRUE)
    }, add = TRUE)
    setwd(fixture)
    here::i_am("code/R/run_all.R")
    wrapper <- new.env(parent = globalenv())
    log <- character()
    wrapper$log_message <- function(...) log <<- c(log, paste(..., sep = ""))
    eval(definition[[1]], wrapper)

    path <- file.path(fixture, "code/R/25_prereview_sensitivities.R")
    success <- c(
      "stopifnot(file.exists('.here'), identical(commandArgs(TRUE), character()))",
      "stopifnot(identical(.libPaths()[1L], normalizePath('fixture-library')))",
      "ms <- 'child replacement'; out <- function(...) NULL; rr <- 'child replacement'",
      "assign('.prereview_child_only', TRUE, envir = globalenv())",
      "RNGkind('Wichmann-Hill', 'Box-Muller', 'Rejection'); set.seed(20260921); runif(1)",
      "writeLines('fresh fixture output', 'refreshed.txt')",
      "cat('fixture module output\\n')",
      "message('fixture diagnostic output')")
    writeLines(success, path)
    # Exercise a caller outside the project root, plus both RNG-presence states.
    setwd(dirname(fixture))
    for (seed_present in c(TRUE, FALSE)) {
      if (seed_present) set.seed(2718)
      else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
        rm(list = ".Random.seed", envir = globalenv())
      previous <- as.list(globalenv(), all.names = TRUE)
      kind <- RNGkind()
      wd <- getwd()
      log <- character()
      wrapper$source_module("25_prereview_sensitivities.R", "25", isolated = TRUE)
      stopifnot(identical(previous, as.list(globalenv(), all.names = TRUE)),
        identical(kind, RNGkind()), identical(wd, getwd()),
        identical(readLines(file.path(fixture, "refreshed.txt")), "fresh fixture output"),
        any(grepl("fixture module output", log, fixed = TRUE)),
        any(grepl("fixture diagnostic output", log, fixed = TRUE)),
        any(grepl("[25]", log, fixed = TRUE)),
        any(grepl("completed successfully", log, fixed = TRUE)))
    }

    for (failure in c("stop('fixture failure')", "quit(status = 7L)")) {
      writeLines(c(success[1:5], "cat('fixture before failure\\n')", failure), path)
      for (seed_present in c(TRUE, FALSE)) {
        if (seed_present) set.seed(2718)
        else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
          rm(list = ".Random.seed", envir = globalenv())
        previous <- as.list(globalenv(), all.names = TRUE)
        kind <- RNGkind()
        wd <- getwd()
        log <- character()
        error <- tryCatch({
          wrapper$source_module("25_prereview_sensitivities.R", "25", isolated = TRUE)
          NULL
        }, error = identity)
        stopifnot(inherits(error, "error"),
          grepl("Failed at module: 25_prereview_sensitivities.R", conditionMessage(error), fixed = TRUE),
          grepl("Isolated Rscript exited with status", conditionMessage(error), fixed = TRUE),
          any(grepl("fixture before failure", log, fixed = TRUE)),
          any(grepl("Pipeline halted.", log, fixed = TRUE)),
          !any(grepl("completed successfully", log, fixed = TRUE)),
          identical(previous, as.list(globalenv(), all.names = TRUE)),
          identical(kind, RNGkind()), identical(wd, getwd()))
      }
    }
    global_fixture <- file.path(fixture, "code/R/global_fixture.R")
    writeLines(".prereview_global_fixture <- 'default source preserved'", global_fixture)
    stopifnot(!exists(".prereview_global_fixture", envir = globalenv(), inherits = FALSE))
    wrapper$source_module("global_fixture.R", "fixture")
    stopifnot(identical(get(".prereview_global_fixture", envir = globalenv()), "default source preserved"))
    rm(list = ".prereview_global_fixture", envir = globalenv())
    stopifnot(identical(fingerprint, tools::md5sum(protected)))
    cat("PASS: 09 -> isolated 25 -> 30; child libraries/stdout/stderr preserved; errors halt; default global sourcing unchanged; isolated globals/RNG/working directory and production outputs preserved.\n")
  }
  run_tests()
})
