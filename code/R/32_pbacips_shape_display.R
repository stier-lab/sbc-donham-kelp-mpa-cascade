# Display retained pBACIPS shapes; never select models or alter primary estimates.
# Rscript code/R/32_pbacips_shape_display.R --validate-only: no fits or outputs.
# --reconstruct explicitly rebuilds the cache. With no flag, only render a cache.
local({
  script_name <- "32_pbacips_shape_display.R"
  source_paths <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
  command_paths <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  paths <- c(source_paths, command_paths)
  paths <- paths[basename(paths) == script_name]
  root <- if (length(paths)) {
    dirname(dirname(dirname(normalizePath(tail(paths, 1L), mustWork = TRUE))))
  } else normalizePath(getwd(), mustWork = TRUE)
  args <- if (any(basename(command_paths) == script_name)) commandArgs(TRUE) else character()
  if (any(!args %in% c("--validate-only", "--check-cache", "--reconstruct")) || length(args) > 1L) {
    stop("Use at most one of --validate-only, --check-cache or --reconstruct")
  }
  mode <- if (length(args)) sub("^--", "", args) else Sys.getenv("MPA_SHAPE_DISPLAY_MODE", "render")
  if (!mode %in% c("render", "validate-only", "check-cache", "reconstruct")) stop("Unknown display mode")

  run_display <- function() {
    source(file.path(root, "code/R/taxon_display.R"), local = TRUE)
    had_rng <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    rng <- if (had_rng) get(".Random.seed", envir = globalenv()) else NULL
    old_options <- options()
    on.exit({
      if (had_rng) assign(".Random.seed", rng, envir = globalenv())
      else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(list = ".Random.seed", envir = globalenv())
      }
      extra <- setdiff(names(options()), names(old_options))
      if (length(extra)) options(stats::setNames(rep(list(NULL), length(extra)), extra))
      options(old_options)
    }, add = TRUE)
    key_columns <- c("MPA", "Taxa", "Source", "Resp")
    key <- function(d) do.call(paste, c(d[key_columns], sep = "|"))
    need <- function(d, columns, label) {
      absent <- setdiff(columns, names(d))
      if (length(absent)) stop(label, " missing columns: ", paste(absent, collapse = ", "))
    }
    input_files <- c(primary = "data/sumstats_final.csv",
      annual = "data/harmonized/harmonized_response_ratios.csv",
      landsat = "data/harmonized/harmonized_landsat_rr.csv",
      selection = "tables/table_model_selection.csv",
      metadata = "data/harmonized/harmonized_site_metadata.csv",
      fitting_code = "code/R/02_pBACIPS_function.R",
      effect_code = "code/R/08_effect_sizes.R",
      constants = "code/R/00c_analysis_constants.R")
    inputs <- file.path(root, input_files)
    if (any(!file.exists(inputs))) stop("Shape display inputs missing: ", paste(input_files[!file.exists(inputs)], collapse = ", "))
    md5 <- stats::setNames(unname(tools::md5sum(inputs)), names(input_files))
    stable_inputs <- function() {
      if (!identical(md5, stats::setNames(unname(tools::md5sum(inputs)), names(input_files)))) {
        stop("Shape display inputs changed during execution; retry after the pipeline finishes")
      }
    }
    primary <- utils::read.csv(inputs[1], stringsAsFactors = FALSE)
    need(primary, c(key_columns, "Primary", "AnalysisType", "Model", "Mean", "SE", "N", "MPA_Start"), "Primary table")
    if (!nrow(primary) || anyNA(primary[key_columns]) || anyDuplicated(key(primary)) ||
        anyNA(primary$Primary) || any(primary$Primary != "Y")) stop("Primary keys must be unique and Primary=Y")
    primary$key <- key(primary)
    if (any(!is.finite(primary$Mean)) || any(!is.finite(primary$SE) | primary$SE <= 0) ||
        any(!is.finite(primary$N) | primary$N < 3 | primary$N != as.integer(primary$N))) stop("Invalid primary effect or N")
    selection <- utils::read.csv(inputs[4], stringsAsFactors = FALSE)
    need(selection, c(key_columns, "Model", "AnalysisType", "Mean", "N"), "Model-selection table")
    selection$key <- key(selection)
    match_selection <- match(primary$key, selection$key)
    if (anyDuplicated(selection$key) || !setequal(primary$key, selection$key) ||
        any(primary$Model != selection$Model[match_selection]) ||
        any(primary$AnalysisType != selection$AnalysisType[match_selection]) ||
        any(primary$N != selection$N[match_selection]) ||
        any(abs(primary$Mean - selection$Mean[match_selection]) > 1e-10)) {
      stop("Primary and model-selection tables disagree")
    }
    metadata <- unique(utils::read.csv(inputs[5], stringsAsFactors = FALSE)[c("CA_MPA_Name_Short", "MPA_Start")])
    if (anyDuplicated(metadata$CA_MPA_Name_Short) ||
        anyNA(match(primary$MPA, metadata$CA_MPA_Name_Short)) ||
        any(primary$MPA_Start != metadata$MPA_Start[match(primary$MPA, metadata$CA_MPA_Name_Short)])) {
      stop("Primary designation years disagree with metadata")
    }
    taxa <- c("Panulirus interruptus" = "P. interruptus", "Semicossyphus pulcher" = "S. pulcher",
      "Strongylocentrotus purpuratus" = "S. purpuratus", "Mesocentrotus franciscanus" = "M. franciscanus",
      "Macrocystis pyrifera" = "M. pyrifera")
    standardize <- function(path) {
      d <- utils::read.csv(path, stringsAsFactors = FALSE)
      need(d, c("CA_MPA_Name_Short", "y", "source", "resp", "year", "BA", "lnDiff", "time"), basename(path))
      out <- data.frame(MPA = d$CA_MPA_Name_Short, Taxa = unname(taxa[d$y]), Source = d$source,
        Resp = d$resp, year = d$year, BA = d$BA, lnRR = d$lnDiff, input_time = d$time, stringsAsFactors = FALSE)
      if (anyNA(out$Taxa)) stop("Unknown annual taxon")
      out$key <- key(out)
      out
    }
    annual <- rbind(standardize(inputs[2]), standardize(inputs[3]))
    if (anyDuplicated(annual[c("key", "year")])) stop("Duplicate source-specific annual keys")
    annual_n <- table(annual$key[is.finite(annual$lnRR)])
    primary_n <- unname(annual_n[primary$key])
    if (anyNA(primary_n) || any(primary_n != primary$N)) {
      stop("Annual input/primary N mismatch in the all-primary membership ledger")
    }
    primary <- primary[order(primary$MPA, primary$Source, primary$Resp, primary$Taxa), ]
    selected <- primary[primary$AnalysisType == "pBACIPS", ]
    if (!nrow(selected) || any(!selected$Model %in% c("Step", "Linear", "Sigmoid", "Asymptotic"))) {
      stop("Unsupported retained pBACIPS model; revise the display deliberately")
    }
    # Extract definitions as R expressions, without sourcing the effect-size pipeline.
    load_definition <- function(file, name, envir) {
      expressions <- parse(file)
      hit <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
        identical(x[[2]], as.name(name)), as.list(expressions))
      if (length(hit) != 1L) stop("Expected one definition of ", name)
      eval(hit[[1]], envir = envir)
      get(name, envir = envir, inherits = FALSE)
    }
    definitions <- new.env(parent = asNamespace("stats"))
    add_time <- load_definition(inputs[7], "add_time_columns", definitions)
    horizon <- load_definition(inputs[8], "EFFECT_SIZE_TIME_YEARS", definitions)
    if (length(horizon) != 1L || !is.finite(horizon) || horizon <= 0) stop("Invalid effect-size horizon")
    observations <- do.call(rbind, lapply(seq_len(nrow(selected)), function(i) {
      row <- selected[i, ]
      d <- annual[annual$key == row$key, ]
      d <- d[order(d$year), ]
      if (nrow(d) != row$N || any(!is.finite(d$lnRR)) || any(!is.finite(d$year)) ||
          any(d$year != as.integer(d$year)) || !all(c("Before", "After") %in% d$BA) ||
          any(d$BA != ifelse(d$year <= row$MPA_Start, "Before", "After"))) {
        stop("Annual input/primary N or designation mismatch: ", row$key)
      }
      d$CA_MPA_Name_Short <- d$MPA
      d <- add_time(d)
      d$CA_MPA_Name_Short <- NULL
      d$index_time <- d$time.model
      # Exact explicit routes in script 08; not a rule chosen from the outcome.
      uses_input_time <- row$Source == "Landsat" || row$key == "Campus Point SMCA|P. interruptus|LTER|Bio"
      d$clock <- if (uses_input_time) "supplied_time" else "observation_index"
      d$clock_units <- if (uses_input_time) "calendar_years_since_designation" else "after_observation_index"
      if (uses_input_time) d$time.model <- d$input_time
      if (any(!is.finite(d$time.model)) || any(d$time.model < 0) ||
          any(d$time.model[d$BA == "Before"] != 0)) stop("Invalid primary model clock: ", row$key)
      if (uses_input_time && any(d$time.model != pmax(d$year - row$MPA_Start, 0))) {
        stop("Supplied-time route is not the documented calendar clock: ", row$key)
      }
      d$lnDiff <- d$lnRR
      d$Model <- row$Model
      d$MPA_Start <- row$MPA_Start
      d$missing_years <- paste(setdiff(seq.int(min(d$year), max(d$year)), d$year), collapse = ";")
      d
    }))
    nonlinear <- selected[selected$Model %in% c("Sigmoid", "Asymptotic"), ]
    counts <- as.data.frame(table(selected$Model), stringsAsFactors = FALSE)
    names(counts) <- c("Model", "n_histories")
    counts$n_pbacips <- nrow(selected)
    counts$n_primary <- nrow(primary)
    counts$primary_input_md5 <- md5[["primary"]]
    stable_inputs()
    message("Retained labels: ", nrow(primary), " primary; ", nrow(selected), " pBACIPS; ",
      paste(paste(counts$Model, counts$n_histories, sep = "="), collapse = ", "))
    if (mode == "validate-only") return(invisible(NULL))
    cache_path <- file.path(root, "data/cache/pbacips_shape_display.rds")
    table_dir <- file.path(root, "tables")
    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    tolerance <- 1e-7
    if (mode == "reconstruct") {
      packages <- c("minpack.lm", "nls2")
      absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
      if (length(absent)) stop("Required reconstruction packages missing: ", paste(absent, collapse = ", "))
      # Load only function assignments. No libraries, global log or model-selection call.
      for (expression in as.list(parse(inputs[6]))) {
        if (is.call(expression) && identical(expression[[1]], as.name("<-")) &&
            is.call(expression[[3]]) && identical(expression[[3]][[1]], as.name("function"))) {
          eval(expression, envir = definitions)
        }
      }
      definitions$MODEL_FIT_LOG <- list()
      definitions$nls.lm <- minpack.lm::nls.lm
      definitions$nls.lm.control <- minpack.lm::nls.lm.control
      definitions$nls2 <- nls2::nls2
      prediction_rows <- audit_rows <- vector("list", nrow(selected))
      for (i in seq_len(nrow(selected))) {
        row <- selected[i, ]
        d <- observations[observations$key == row$key, ]
        d$BA <- factor(d$BA, levels = c("Before", "After"))
        warnings <- character()
        effect <- NA_real_
        outcome <- tryCatch(withCallingHandlers({
          linear <- if (row$Model != "Step") {
            stats::lm(stats::as.formula("lnDiff ~ time.model", env = asNamespace("stats")), data = d)
          } else NULL
          fit <- switch(row$Model,
            Linear = linear,
            Step = stats::lm(stats::as.formula("lnDiff ~ BA", env = asNamespace("stats")), data = d),
            Sigmoid = definitions$mySIGfun_standalone(d$lnRR, d$time.model, max(which(d$time.model == 0)), d$time.true),
            Asymptotic = definitions$myASYfun_standalone(d$lnRR, d$time.model, max(which(d$time.model == 0)), d$time.true))
          if (is.null(fit) || any(!is.finite(stats::coef(fit)))) stop("Retained fit failed")
          fallback <- attr(fit, "model_fallback")
          if (row$Model %in% c("Sigmoid", "Asymptotic") &&
              (!inherits(fit, "nls") || (!is.null(fallback) && !fallback %in% c("4pl_reparam", "selfstart", "stable_reparam")) ||
               (!is.null(fit$convInfo) && !isTRUE(fit$convInfo$isConv)))) {
            stop("Nonconverged/non-NLS substitute rejected")
          }
          newdata <- function(t, ba = "After") data.frame(time.model = t, time_offset = t + 0.01,
            time = t, time.true = max(which(d$time.model == 0)) + t,
            BA = factor(rep(ba, length(t)), levels = c("Before", "After")))
          effect <- if (row$Model == "Step") {
            as.numeric(stats::predict(fit, newdata(0, "After")) - stats::predict(fit, newdata(0, "Before")))
          } else diff(as.numeric(stats::predict(fit, newdata(c(0, horizon)))))
          if (!is.finite(effect) || abs(effect - row$Mean) > tolerance) {
            stop("Reconstructed contrast does not match the primary effect")
          }
          # Calendar display coordinates preserve each primary route's native clock.
          post <- d[d$BA == "After", ]
          years <- sort(unique(c(d$year, seq(min(d$year), max(d$year), length.out = 401L), row$MPA_Start + c(0.49, 0.51))))
          times <- if (d$clock[1] == "supplied_time") pmax(years - row$MPA_Start, 0) else
            stats::approx(post$year, post$time.model, xout = years, rule = 2)$y
          nd <- newdata(times)
          nd$BA <- factor(ifelse(years <= row$MPA_Start + 0.5, "Before", "After"), levels = c("Before", "After"))
          pred <- row[rep(1L, length(years)), c(key_columns, "key", "Model")]
          pred$year <- years
          pred$time.model <- times
          pred$clock <- d$clock[1]
          pred$clock_units <- d$clock_units[1]
          pred$BA <- as.character(nd$BA)
          pred$curve <- "Retained fit"
          pred$fitted_lnRR <- as.numeric(stats::predict(fit, newdata = nd))
          if (row$Model %in% c("Sigmoid", "Asymptotic")) {
            comparator <- pred
            comparator$curve <- "Linear comparator"
            comparator$fitted_lnRR <- as.numeric(stats::predict(linear, newdata = nd))
            pred <- rbind(pred, comparator)
          }
          if (any(!is.finite(pred$fitted_lnRR))) stop("Nonfinite curve predictions")
          audit <- data.frame(row[key_columns], key = row$key, Model = row$Model,
            N = nrow(d), primary_N = row$N, primary_effect = row$Mean, reconstructed_effect = effect,
            difference = effect - row$Mean, horizon = if (row$Model == "Step") NA_real_ else horizon,
            horizon_within_observed_clock = if (row$Model == "Step") NA else max(d$time.model) >= horizon,
            contrast_definition = if (row$Model == "Step") "After minus Before" else paste0("f(model time ", horizon, ") minus f(model time 0)"),
            calendar_year_at_horizon = if (row$Model == "Step" || !any(d$time.model[d$BA == "After"] == horizon)) NA_real_ else
              d$year[d$BA == "After" & d$time.model == horizon][1L],
            year_min = min(d$year), year_max = max(d$year), first_after_year = min(post$year),
            model_time_max = max(d$time.model), first_after_model_time = min(post$time.model),
            n_before = sum(d$BA == "Before"), n_after = sum(d$BA == "After"),
            missing_years = d$missing_years[1], clock = d$clock[1], clock_units = d$clock_units[1],
            parameter_names = paste(names(stats::coef(fit)), collapse = ";"),
            parameter_values = paste(format(stats::coef(fit), digits = 17), collapse = ";"),
            fit_class = class(fit)[1], fallback = if (is.null(fallback)) "none" else fallback,
            status = "verified", error = "", stringsAsFactors = FALSE)
          list(predictions = pred, audit = audit)
        }, warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }), error = identity)
        if (inherits(outcome, "error")) {
          audit_rows[[i]] <- data.frame(row[key_columns], key = row$key, Model = row$Model,
            N = nrow(d), primary_N = row$N, primary_effect = row$Mean, reconstructed_effect = effect,
            difference = effect - row$Mean, horizon = if (row$Model == "Step") NA_real_ else horizon,
            horizon_within_observed_clock = if (row$Model == "Step") NA else max(d$time.model) >= horizon,
            contrast_definition = if (row$Model == "Step") "After minus Before" else paste0("f(model time ", horizon, ") minus f(model time 0)"),
            calendar_year_at_horizon = if (row$Model == "Step" || !any(d$time.model[d$BA == "After"] == horizon)) NA_real_ else
              d$year[d$BA == "After" & d$time.model == horizon][1L],
            year_min = min(d$year), year_max = max(d$year), first_after_year = min(d$year[d$BA == "After"]),
            model_time_max = max(d$time.model), first_after_model_time = min(d$time.model[d$BA == "After"]),
            n_before = sum(d$BA == "Before"), n_after = sum(d$BA == "After"),
            missing_years = d$missing_years[1], clock = d$clock[1], clock_units = d$clock_units[1],
            parameter_names = "", parameter_values = "", fit_class = "",
            fallback = "", status = "failed", error = conditionMessage(outcome), stringsAsFactors = FALSE)
        } else {
          prediction_rows[[i]] <- outcome$predictions
          audit_rows[[i]] <- outcome$audit
        }
        audit_rows[[i]]$warnings <- paste(unique(warnings), collapse = " | ")
      }
      audit <- do.call(rbind, audit_rows)
      stable_inputs()
      utils::write.csv(audit, file.path(table_dir, "table_s_pbacips_shape_fit_audit.csv"), row.names = FALSE)
      if (any(audit$status != "verified")) stop("Shape reconstruction failed; see table_s_pbacips_shape_fit_audit.csv. No cache or figures replaced.")
      cache <- list(version = 3L, input_md5 = md5, observations = observations,
        predictions = do.call(rbind, prediction_rows), audit = audit)
      dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
      saveRDS(cache, cache_path)
    } else {
      if (!file.exists(cache_path)) stop("Verified shape cache absent. After the full pipeline finishes, explicitly run --reconstruct; ordinary execution never fits models.")
      cache <- readRDS(cache_path)
    }
    if (!identical(cache$version, 3L) || !identical(cache$input_md5, md5)) stop("Stale shape cache; input/code MD5 changed. Explicit reconstruction required after the pipeline finishes.")
    if (!isTRUE(all.equal(cache$observations, observations, check.attributes = FALSE, tolerance = 1e-12))) stop("Cached observations differ from harmonized inputs")
    audit <- cache$audit
    predictions <- cache$predictions
    need(audit, c("key", "Model", "N", "primary_effect", "reconstructed_effect", "difference", "status", "clock", "clock_units"), "Fit audit")
    need(predictions, c(key_columns, "key", "Model", "year", "time.model", "BA", "curve", "fitted_lnRR", "clock", "clock_units"), "Predictions")
    a <- match(selected$key, audit$key)
    if (anyDuplicated(audit$key) || !setequal(audit$key, selected$key) || any(audit$status != "verified") ||
        anyNA(a) || any(audit$Model[a] != selected$Model) || any(audit$N[a] != selected$N) ||
        any(!is.finite(audit$reconstructed_effect[a])) || any(abs(audit$primary_effect[a] - selected$Mean) > tolerance) ||
        any(abs(audit$reconstructed_effect[a] - selected$Mean) > tolerance)) stop("Cached fit audit does not match primary effects")
    retained <- predictions[predictions$curve == "Retained fit", ]
    comparator <- predictions[predictions$curve == "Linear comparator", ]
    if (!setequal(retained$key, selected$key) || !setequal(comparator$key, nonlinear$key) ||
        any(!predictions$curve %in% c("Retained fit", "Linear comparator")) ||
        any(!is.finite(predictions$fitted_lnRR)) || any(!is.finite(predictions$year)) ||
        any(!is.finite(predictions$time.model)) || anyDuplicated(predictions[c("key", "curve", "year")]) ||
        any(predictions$key != key(predictions))) stop("Invalid cached prediction cohort")
    for (i in seq_len(nrow(selected))) {
      d <- observations[observations$key == selected$key[i], ]
      p <- predictions[predictions$key == selected$key[i], ]
      after <- d[d$BA == "After", ]
      expected_time <- if (d$clock[1] == "supplied_time") pmax(p$year - selected$MPA_Start[i], 0) else
        stats::approx(after$year, after$time.model, xout = p$year, rule = 2)$y
      if (any(p$Model != selected$Model[i]) || any(p$year < min(d$year) | p$year > max(d$year)) ||
          any(p$time.model < 0 | p$time.model > max(d$time.model)) ||
          any(abs(p$time.model - expected_time) > 1e-10) ||
          any(p$BA != ifelse(p$year <= selected$MPA_Start[i] + 0.5, "Before", "After")) ||
          any(p$clock != d$clock[1] | p$clock_units != d$clock_units[1]) ||
          audit$clock[a[i]] != d$clock[1] || audit$clock_units[a[i]] != d$clock_units[1]) {
        stop("Cached curve exceeds observed support or native clock/model labels")
      }
    }
    stable_inputs()
    if (mode == "check-cache") { message("Shape cache verified; no fits or outputs"); return(invisible(NULL)) }
    packages <- c("ggplot2", "patchwork", "ragg", "scales", "ggrepel")
    absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(absent)) stop("Required rendering packages missing: ", paste(absent, collapse = ", "))
    membership <- primary[c(key_columns, "key", "Model", "AnalysisType", "N", "MPA_Start")]
    membership$display_scope <- ifelse(membership$AnalysisType == "pBACIPS", "all_pbacips", "not_displayed_ci")
    membership$nonlinear_comparison <- membership$key %in% nonlinear$key
    membership$nonlinear_panel <- LETTERS[match(membership$key, nonlinear$key)]
    membership$n_finite_annual_inputs <- unname(annual_n[membership$key])
    membership$annual_N_matches_primary <- membership$n_finite_annual_inputs == membership$N
    member_index <- match(membership$key, selected$key)
    membership$page <- ifelse(is.na(member_index), NA_integer_, (member_index - 1L) %/% 4L + 1L)
    membership$panel <- ifelse(is.na(member_index), NA_character_, LETTERS[(member_index - 1L) %% 4L + 1L])
    membership$primary_input_md5 <- md5[["primary"]]
    membership$clock <- observations$clock[match(membership$key, observations$key)]
    membership$clock_units <- observations$clock_units[match(membership$key, observations$key)]
    make_panel <- function(row, letter, show_comparator = FALSE) {
      d <- observations[observations$key == row$key, ]
      p <- predictions[predictions$key == row$key & (show_comparator | predictions$curve == "Retained fit"), ]
      mpa_label <- paste(strwrap(row$MPA, width = 22L), collapse = "\n")
      y_label <- paste(mpa_label, paste(taxon_display(row$Taxa), ifelse(row$Resp == "Bio", "biomass", "density")),
        paste(row$Source, "lnRR (unitless)"), sep = "\n")
      last <- p[!duplicated(p$curve, fromLast = TRUE), ]
      last$label <- ifelse(last$curve == "Retained fit", row$Model, "Linear comparator")
      plot <- ggplot2::ggplot(d, ggplot2::aes(year, lnRR)) +
        ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.25, linetype = "dashed") +
        ggplot2::geom_vline(xintercept = row$MPA_Start + 0.5, color = "grey55", linewidth = 0.3, linetype = "dotted") +
        ggplot2::geom_line(data = p, ggplot2::aes(y = fitted_lnRR, color = curve, linetype = curve), linewidth = 0.65) +
        ggplot2::geom_point(ggplot2::aes(color = BA, shape = BA), size = 1.6, alpha = 0.85) +
        ggrepel::geom_text_repel(data = last, ggplot2::aes(x = max(d$year) + 0.5, y = fitted_lnRR, label = label, color = curve),
          inherit.aes = FALSE, hjust = 0, size = 2.6, direction = "y", seed = 32L,
          box.padding = 0.2, point.padding = 0.2, max.overlaps = Inf, max.iter = 10000,
          max.time = Inf, min.segment.length = 0, segment.size = 0.25) +
        ggplot2::scale_color_manual(values = c(Before = "#0072B2", After = "#D55E00", "Retained fit" = "#202020", "Linear comparator" = "#009E73")) +
        ggplot2::scale_shape_manual(values = c(Before = 16, After = 17)) +
        ggplot2::scale_linetype_manual(values = c("Retained fit" = "solid", "Linear comparator" = "dashed")) +
        ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 5), expand = ggplot2::expansion(mult = c(0.02, 0.23))) +
        ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(n = 4), expand = ggplot2::expansion(mult = c(0.1, 0.1))) +
        ggplot2::labs(x = "Year", y = y_label, tag = letter) +
        ggplot2::theme_classic(base_size = 9, base_family = "sans") +
        ggplot2::theme(legend.position = "none", plot.tag = ggplot2::element_text(face = "bold", size = 11),
          plot.tag.position = "topleft", plot.margin = ggplot2::margin(5, 4, 4, 4, "mm"),
          axis.title.y = ggplot2::element_text(size = 8, lineheight = 0.95),
          axis.title.x = ggplot2::element_text(size = 9), axis.text = ggplot2::element_text(size = 8))
      plot
    }
    plot_dir <- file.path(root, "plots/pbacips_shape_histories")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
    save_plot <- function(plot, stem, height) {
      ggplot2::ggsave(paste0(stem, ".pdf"), plot, device = pdf_device, width = 170, height = height, units = "mm", bg = "white")
      ggplot2::ggsave(paste0(stem, ".png"), plot, device = ragg::agg_png, width = 170, height = height, units = "mm", dpi = 300, bg = "white")
    }
    pages <- split(seq_len(nrow(selected)), (seq_len(nrow(selected)) - 1L) %/% 4L + 1L)
    page_plots <- lapply(seq_along(pages), function(page) {
      rows <- pages[[page]]
      panels <- lapply(seq_along(rows), function(j) make_panel(selected[rows[j], ], LETTERS[j]))
      p <- patchwork::wrap_plots(panels, ncol = 1)
      save_plot(p, file.path(plot_dir, sprintf("page_%03d", page)), 210)
      p
    })
    combined_pdf <- file.path(root, "plots/fig_s_pbacips_selected_histories.pdf")
    pdf_device(combined_pdf, width = 170 / 25.4, height = 210 / 25.4, onefile = TRUE)
    tryCatch(for (p in page_plots) print(p), finally = grDevices::dev.off())
    if (nrow(nonlinear)) {
      panels <- lapply(seq_len(nrow(nonlinear)), function(i) make_panel(nonlinear[i, ], LETTERS[i], TRUE))
      save_plot(patchwork::wrap_plots(panels, ncol = 1), file.path(root, "plots/fig_s_pbacips_nonlinear_comparison"), 65 * nrow(nonlinear))
    }
    outputs <- list(membership = membership, observations = observations, predictions = predictions, fit_audit = audit, counts = counts)
    stable_inputs()
    for (name in names(outputs)) {
      utils::write.csv(outputs[[name]], file.path(table_dir, paste0("table_s_pbacips_shape_", name, ".csv")), row.names = FALSE)
    }
    provenance <- data.frame(input = input_files, md5 = unname(md5), stringsAsFactors = FALSE)
    utils::write.csv(provenance, file.path(table_dir, "table_s_pbacips_shape_provenance.csv"), row.names = FALSE)
    message("Shape display saved: all ", nrow(selected), " retained pBACIPS histories and all ", nrow(nonlinear), " nonlinear comparisons. No model selection or primary estimates changed.")
  }
  run_display()
})
