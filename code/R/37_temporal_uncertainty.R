# Calendar-year AR1 on the exact Figure 4 rows, including crossed source/MPA terms.
local({
  root <- here::here()
  source(file.path(root, "code/R/temporal_uncertainty_utils.R"), local = TRUE)
  cache <- readRDS(file.path(root, "data/cache/temporal_primary_models.rds"))
  input <- read.csv(file.path(root, "data/harmonized/harmonized_response_ratios.csv"))
  primary <- read.csv(file.path(root, "data/sumstats_final.csv"))
  expected <- temporal_primary_data(input, primary)
  stopifnot(identical(as.character(expected$history), as.character(cache$annual$history)),
    identical(expected$year, cache$annual$year), identical(expected$lnDiff, cache$annual$lnDiff))
  d <- cache$annual[cache$annual$resp == "Bio", ]
  d$Species <- factor(d$Species, levels = levels(cache$bio@frame$Species))
  d$history <- factor(d$history)
  d$year_f <- factor(d$year, levels = seq.int(min(d$year), max(d$year)))
  stopifnot(nrow(d) == nobs(cache$bio), !anyDuplicated(d[c("history", "year")]),
    !lme4::isSingular(cache$bio))
  candidates <- c(correlated_slopes_source = "(1 + time | CA_MPA_Name_Short) + (1 | source)",
    independent_slopes_source = "(1 + time || CA_MPA_Name_Short) + (1 | source)",
    intercepts_source = "(1 | CA_MPA_Name_Short) + (1 | source)",
    intercepts = "(1 | CA_MPA_Name_Short)")
  fits <- audits <- vector("list", length(candidates))
  for (i in seq_along(candidates)) {
    warnings <- character()
    f <- as.formula(paste("lnDiff ~ Species * time +", candidates[i], "+ ar1(year_f + 0 | history)"))
    m <- tryCatch(withCallingHandlers(glmmTMB::glmmTMB(f, data = d, family = gaussian(), REML = TRUE,
      control = glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 1000, eval.max = 1500))),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) e)
    valid <- !inherits(m, "error") && m$fit$convergence == 0 && isTRUE(m$sdr$pdHess) &&
      all(is.finite(glmmTMB::fixef(m)$cond)) && all(is.finite(vcov(m)$cond))
    variance <- if (valid) glmmTMB::VarCorr(m)$cond else NULL
    process_sd <- if (valid) attr(variance$history, "stddev")[1] else NA_real_
    phi <- if (valid) attr(variance$history, "correlation")[1, 2] else NA_real_
    error_sd <- if (valid) sigma(m) else NA_real_
    audits[[i]] <- data.frame(candidate = names(candidates)[i], valid,
      convergence = if (inherits(m, "error")) NA_integer_ else m$fit$convergence,
      pdHess = if (inherits(m, "error")) FALSE else isTRUE(m$sdr$pdHess),
      n_obs = nrow(d), n_history = length(unique(d$history)), n_mpa = length(unique(d$CA_MPA_Name_Short)),
      logLik = if (valid) as.numeric(logLik(m)) else NA_real_,
      max_abs_gradient = if (valid) max(abs(m$sdr$gradient.fixed)) else NA_real_,
      source_sd = if (valid && "source" %in% names(variance)) attr(variance$source, "stddev")[1] else NA_real_,
      ar1_process_sd = process_sd, annual_error_sd = error_sd, ar1_process_phi = phi,
      conditional_lag1_correlation = phi * process_sd^2 / (process_sd^2 + error_sd^2),
      warnings = if (inherits(m, "error")) conditionMessage(m) else paste(unique(warnings), collapse = "; "),
      selected = FALSE)
    if (valid) fits[[i]] <- m
    message("Temporal AR1 candidate ", names(candidates)[i], ": valid=", valid)
  }
  audit <- do.call(rbind, audits)
  write.csv(audit, file.path(root, "tables/table_s_temporal_ar1_attempts.csv"), row.names = FALSE)
  selected <- which(audit$valid)[1]
  if (is.na(selected)) stop("No valid calendar-AR1 temporal fit; inspect recorded attempts.")
  audit$selected[selected] <- TRUE
  m <- fits[[selected]]
  beta <- glmmTMB::fixef(m)$cond; v <- vcov(m)$cond
  vc <- glmmTMB::VarCorr(m)$cond$history
  rho <- attr(vc, "correlation")[1, 2]
  grid <- do.call(rbind, lapply(levels(d$Species), function(sp) {
    nd <- data.frame(Species = factor(sp, levels = levels(d$Species)), time = seq(0, 15, by = .1))
    x <- model.matrix(~ Species * time, nd)[, names(beta), drop = FALSE]
    est <- drop(x %*% beta); se <- sqrt(rowSums((x %*% v) * x))
    data.frame(Species = sp, time = nd$time, estimate = est, se,
      ci_low = est - qnorm(.975) * se, ci_high = est + qnorm(.975) * se,
      inference = "conditional model-based Wald; not a small-MPA calibration")
  }))
  primary_grid <- read.csv(file.path(root, "tables/table_s_temporal_predictions.csv"))
  primary_grid <- primary_grid[primary_grid$resp == "Bio", ]
  comparison <- merge(primary_grid, grid, by = c("Species", "time"), suffixes = c("_primary", "_ar1"))
  comparison$ar1_phi <- rho
  broad <- input[input$BA == "After" & input$time >= 0 & input$resp == "Bio", ]
  broad$Species <- factor(broad$y, levels = levels(d$Species))
  old <- lme4::lmer(lnDiff ~ Species * time + (1 | CA_MPA_Name_Short) + (1 | source), broad)
  cohort <- do.call(rbind, lapply(levels(d$Species), function(sp) {
    rbind(cbind(cohort = "retained diver histories", n_obs = nrow(d), Species = sp,
      temporal_contrast(cache$bio, sp, slope = TRUE)),
      cbind(cohort = "former broader histories/intercept fallback", n_obs = nrow(broad), Species = sp,
        temporal_contrast(old, sp, slope = TRUE)))
  }))
  slopes <- do.call(rbind, lapply(levels(d$Species), function(sp) {
    nd <- data.frame(Species = factor(sp, levels = levels(d$Species)), time = c(0, 1))
    x <- model.matrix(~ Species * time, nd)[, names(beta), drop = FALSE]
    g <- x[2, ] - x[1, ]; estimate <- sum(g * beta); se <- sqrt(drop(t(g) %*% v %*% g))
    base <- temporal_contrast(cache$bio, sp, slope = TRUE)
    data.frame(Species = sp, estimate, se, ci_low = estimate - qnorm(.975) * se,
      ci_high = estimate + qnorm(.975) * se, primary_slope = base$estimate, primary_se = base$se,
      primary_ci_low = base$ci_low, primary_ci_high = base$ci_high,
      sign_changed = sign(estimate) != sign(base$estimate), ar1_phi = rho,
      primary_df = base$df, uncertainty = "AR1 model-based Wald; primary exact-contrast Satterthwaite")
  }))
  output <- function(x, name) write.csv(x, file.path(root, "tables", paste0(name, ".csv")), row.names = FALSE)
  output(audit, "table_s_temporal_ar1_attempts")
  output(data.frame(statistic = NA_real_,
    note = "Retired pooled DW/nested reconstruction; use table_s_temporal_ar1_slopes and its exact retained-diver cohort."),
    "table_s_ar1_sensitivity")
  output(slopes, "table_s_temporal_ar1_slopes")
  output(comparison, "table_s_temporal_ar1_predictions")
  output(cohort, "table_s_temporal_cohort_comparison")
  output(d[c("history", "CA_MPA_Name_Short", "source", "resp", "y", "year", "time", "lnDiff")],
    "table_s_temporal_primary_annual")
  saveRDS(list(ar1 = m, primary = cache$bio, annual = d, attempts = audit),
    file.path(root, "data/cache/temporal_ar1_check.rds"))
  cat("Calendar AR1 sensitivity uses the exact Figure 4 annual rows; process rho=", rho, "\n")
  print(slopes)
})
