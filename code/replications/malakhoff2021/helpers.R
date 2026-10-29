urchin_mass <- function(size, species) {
  stopifnot(all(is.finite(size) & size > 0), all(species %in% c("red", "purple")))
  0.00059 * size^ifelse(species == "red", 2.917, 2.870)
}

size_class <- function(size, definition) {
  if (definition == "code") {
    ifelse(size < 36, "small", ifelse(size < 51, "medium", ifelse(size < 71, "large", "very_large")))
  } else {
    ifelse(size < 35, "small", ifelse(size <= 50, "medium", ifelse(size <= 70, "large", "very_large")))
  }
}

fit_baci <- function(input, endpoint, end_year, first_post, formulation, offset = 1) {
  d <- as.data.frame(input[input$year >= 1991 & input$year <= end_year, ])
  d <- d[is.finite(d$value), ]
  d <- d[order(d$site, d$year), ]
  d$site <- factor(d$site)
  d$period <- factor(ifelse(d$year >= first_post, "After", "Before"),
                     levels = c("Before", "After"))
  d$status <- factor(ifelse(d$site %in% sprintf("a-k-%02d", c(2, 6, 9, 14)),
                            "Reserve", "Control"), levels = c("Control", "Reserve"))
  d$response <- log(d$value + offset)
  formula <- response ~ status * period
  warnings <- character()
  result <- tryCatch(withCallingHandlers({
    if (formulation == "original_year_site") {
      nlme::lme(formula, random = ~1 | year/site, correlation = nlme::corAR1(),
                data = d, control = nlme::lmeControl(msMaxIter = 200))
    } else if (formulation == "site_independent") {
      nlme::lme(formula, random = ~1 | site, data = d,
                control = nlme::lmeControl(msMaxIter = 200))
    } else {
      nlme::lme(formula, random = ~1 | site,
                correlation = nlme::corAR1(form = ~year | site), data = d,
                control = nlme::lmeControl(msMaxIter = 200))
    }
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  row <- data.frame(endpoint = endpoint, end_year = end_year,
                    first_post = first_post, formulation = formulation,
                    n = nrow(d), n_sites = length(unique(d$site)),
                    n_reserve_sites = length(unique(d$site[d$status == "Reserve"])),
                    estimate = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_,
                    p = NA_real_, phi = NA_real_, min_random_sd = NA_real_,
                    residual_lag1 = NA_real_, status = "failed", diagnostic = "",
                    warnings = paste(unique(warnings), collapse = "; "))
  if (inherits(result, "error")) {
    row$diagnostic <- conditionMessage(result)
    return(list(row = row, fit = NULL, data = d))
  }
  # Original single-species year/site groups have one annual observation each:
  # their AR parameter does not represent within-site temporal dependence.
  tt <- summary(result)$tTable
  term <- "statusReserve:periodAfter"
  stopifnot(term %in% rownames(tt))
  df <- tt[term, "DF"]
  row$estimate <- tt[term, "Value"]
  row$se <- tt[term, "Std.Error"]
  row$lower <- row$estimate - stats::qt(.975, df) * row$se
  row$upper <- row$estimate + stats::qt(.975, df) * row$se
  row$p <- tt[term, "p-value"]
  if (!is.null(result$modelStruct$corStruct)) {
    row$phi <- as.numeric(stats::coef(result$modelStruct$corStruct, unconstrained = FALSE))[1]
  }
  sd <- suppressWarnings(as.numeric(nlme::VarCorr(result)[, "StdDev"]))
  row$min_random_sd <- min(sd[is.finite(sd)], na.rm = TRUE)
  d$residual <- residuals(result, type = "normalized")
  residual_pairs <- lapply(split(d, d$site), function(x) {
    if (nrow(x) < 2) return(NULL)
    ok <- diff(x$year) == 1
    data.frame(previous = head(x$residual, -1)[ok], next_year = tail(x$residual, -1)[ok])
  })
  rp <- do.call(rbind, residual_pairs)
  if (nrow(rp) > 2) row$residual_lag1 <- stats::cor(rp$previous, rp$next_year)
  row$status <- if (formulation == "original_year_site") "not_temporal_AR1" else "fitted"
  row$diagnostic <- if (formulation == "original_year_site")
    "AR1 grouped within singleton year/site cells, not a within-site time series" else
    if (row$min_random_sd < 1e-4) "near-zero variance component; inspect" else "see residual diagnostics"
  list(row = row, fit = result, data = d)
}
