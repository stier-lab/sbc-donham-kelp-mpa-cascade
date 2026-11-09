# Review-only period refits on corrected annual ratios. Calendar-year gaps
# enter the AR1 factor explicitly; source is a fixed effect. These tests do
# not replace the manuscript's original heatwave models or identify cause.
suppressPackageStartupMessages(library(glmmTMB))
root <- normalizePath(getwd(), mustWork = TRUE)
folder <- file.path(root, "outputs/whole_dataset_biomass_update_20261105")
d <- read.csv(file.path(root, "outputs/integrated_review_20261103/annual_values.csv"))
d <- d[d$response_supported == "True" & is.finite(d$ln_ratio), ]
response <- c("M. pyrifera" = "Bio", "M. franciscanus" = "Den",
              "S. purpuratus" = "Den", "P. interruptus" = "Den", "S. pulcher" = "Den")
rows <- list()
for (taxon in names(response)) {
  z <- d[d$Taxa == taxon & d$Resp == response[[taxon]], ]
  z$period <- factor(ifelse(z$year <= 2013, "before", ifelse(z$year <= 2016, "during", "after")),
                     levels = c("before", "during", "after"))
  z$year_factor <- factor(z$year, levels = seq(min(z$year), max(z$year)))
  z$MPA <- factor(z$MPA)
  z$Source <- factor(z$Source)
  z$series_id <- factor(z$series_id)
  status <- if (length(unique(z$period)) < 3L) "missing_period" else "fitted"
  fit <- NULL
  warning_text <- character()
  if (status == "fitted") {
    fit <- tryCatch(withCallingHandlers(glmmTMB(
      ln_ratio ~ t + period + Source + (1 | MPA) + ar1(year_factor + 0 | series_id),
      data = z, REML = TRUE,
      control = glmmTMBControl(optCtrl = list(iter.max = 10000, eval.max = 10000))),
      warning = function(w) { warning_text <<- c(warning_text, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) e)
    if (inherits(fit, "error")) { status <- paste0("fit_failed: ", conditionMessage(fit)); fit <- NULL }
    else if (!isTRUE(fit$sdr$pdHess) || fit$fit$convergence != 0L) status <- "unstable_fit"
  }
  co <- if (is.null(fit)) NULL else summary(fit)$coefficients$cond
  for (term in c("periodduring", "periodafter", "t")) {
    ok <- !is.null(co) && term %in% rownames(co)
    estimate <- if (ok) as.numeric(co[term, "Estimate"]) else NA_real_
    se <- if (ok) as.numeric(co[term, "Std. Error"]) else NA_real_
    rows[[length(rows) + 1L]] <- data.frame(Taxa = taxon, Resp = response[[taxon]],
      term = term, status = status, n_annual = nrow(z), n_histories = length(unique(z$series_id)),
      n_MPAs = length(unique(z$MPA)), sources = paste(sort(unique(as.character(z$Source))), collapse = ";"),
      estimate = estimate, se = se,
      low = if (ok) estimate - 1.96 * se else NA_real_,
      high = if (ok) estimate + 1.96 * se else NA_real_,
      p = if (ok) as.numeric(co[term, ncol(co)]) else NA_real_,
      warning = paste(unique(warning_text), collapse = "; "),
      scope = "annual protected/reference log ratio; linear recovery plus heatwave period; source fixed effect; MPA intercept; calendar-gap AR1 within history")
  }
}
out <- do.call(rbind, rows)
out$q <- NA_real_
idx <- which(out$term %in% c("periodduring", "periodafter") & is.finite(out$p))
if (length(idx)) out$q[idx] <- p.adjust(out$p[idx], method = "BH", n = 10L)
write.csv(out, file.path(folder, "heatwave_period_refits.csv"), row.names = FALSE)
cat("Heatwave period fits:", sum(out$status == "fitted") / 3, "of", length(response), "taxa\n")
