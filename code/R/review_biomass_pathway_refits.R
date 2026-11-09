# Conditional pathway refit for the corrected integrated review.
# Same physical stations, zones and years are selected in the Python input
# builder. This never replaces the manuscript's Table 3 or frozen figures.
suppressPackageStartupMessages({library(metafor); library(sandwich)})
root <- normalizePath(getwd(), mustWork = TRUE)
source(file.path(root, "code/R/primary_uncertainty_utils.R"), local = TRUE)
folder <- file.path(root, "outputs/whole_dataset_biomass_update_20261105")
annual <- read.csv(file.path(folder, "matched_year_pathway_annual.csv"))
stopifnot(all(c("panel", "MPA", "Source", "year", "x", "y") %in% names(annual)))
stopifnot(!anyDuplicated(annual[c("panel", "MPA", "Source", "year")]))
rows <- list()
for (key in split(annual, interaction(annual$panel, annual$MPA, annual$Source, drop = TRUE))) {
  key <- key[order(key$year), ]
  if (nrow(key) < 3L) next
  centered <- key$year - mean(key$year)
  xfit <- lm(key$x ~ centered)
  yfit <- lm(key$y ~ centered)
  xcov <- tryCatch(uncertainty_hac(xfit, key$year, 3L), error = function(e) NULL)
  ycov <- tryCatch(uncertainty_hac(yfit, key$year, 3L), error = function(e) NULL)
  rows[[length(rows) + 1L]] <- data.frame(panel = key$panel[1], MPA = key$MPA[1],
    Source = key$Source[1], n_years = nrow(key), first_year = min(key$year),
    last_year = max(key$year), x_slope = unname(coef(xfit)[2]),
    y_slope = unname(coef(yfit)[2]),
    x_vi = if (is.null(xcov)) NA_real_ else xcov[2, 2],
    y_vi = if (is.null(ycov)) NA_real_ else ycov[2, 2],
    uncertainty = if (is.null(xcov) || is.null(ycov)) "HAC unavailable" else "calendar-gap HAC lag 3")
}
effects <- do.call(rbind, rows)
write.csv(effects, file.path(folder, "matched_year_pathway_effects.csv"), row.names = FALSE)
fits <- list()
for (panel in letters[1:6]) for (src in sort(unique(effects$Source))) {
  d <- effects[effects$panel == panel & effects$Source == src, ]
  valid <- d[is.finite(d$x_slope) & is.finite(d$y_slope) &
    is.finite(d$y_vi) & d$y_vi > 0, ]
  status <- if (nrow(valid) < 4L) "fewer_than_four_MPAs" else if (
    sd(valid$x_slope) < 1e-12) "constant_predictor" else "fitted"
  model <- if (status == "fitted") tryCatch(
    suppressWarnings(metafor::rma(yi = y_slope, vi = y_vi, mods = ~ x_slope,
      data = valid, method = "REML", test = "knha")), error = function(e) e) else NULL
  if (inherits(model, "error")) { status <- paste0("fit_failed: ", conditionMessage(model)); model <- NULL }
  fits[[length(fits) + 1L]] <- data.frame(panel = panel, Source = src,
    candidate_MPAs = nrow(d), fitted_MPAs = nrow(valid),
    matched_MPA_years = sum(valid$n_years), status = status,
    slope = if (is.null(model)) NA_real_ else unname(coef(model)[2]),
    se = if (is.null(model)) NA_real_ else model$se[2],
    ci_low = if (is.null(model)) NA_real_ else model$ci.lb[2],
    ci_high = if (is.null(model)) NA_real_ else model$ci.ub[2],
    p = if (is.null(model)) NA_real_ else model$pval[2],
    tau2 = if (is.null(model)) NA_real_ else model$tau2,
    scope = "Within-program association of matched-year post-protection slopes; outcome HAC uncertainty; predictor and cross-taxon error not fitted")
}
result <- do.call(rbind, fits)
result$q <- NA_real_
for (src in unique(result$Source)) {
  idx <- which(result$Source == src & is.finite(result$p))
  if (length(idx)) result$q[idx] <- p.adjust(result$p[idx], method = "BH", n = 6L)
}
write.csv(result, file.path(folder, "matched_year_pathway_fits.csv"), row.names = FALSE)
stopifnot(nrow(result) == 6L * length(unique(effects$Source)))
cat("Matched-year pathway fits:", sum(result$status == "fitted"), "of", nrow(result), "\n")
