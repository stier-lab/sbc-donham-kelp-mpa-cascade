# Shared estimators for the designation and historical-disturbance checks.

review_summary <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (!n) return(data.frame(n_mpa = 0L, estimate = NA_real_, se = NA_real_,
                            ci_low = NA_real_, ci_high = NA_real_, p = NA_real_))
  estimate <- mean(x)
  if (n < 2L) return(data.frame(n_mpa = n, estimate = estimate, se = NA_real_,
                               ci_low = NA_real_, ci_high = NA_real_, p = NA_real_))
  se <- stats::sd(x) / sqrt(n)
  critical <- stats::qt(0.975, n - 1L)
  p <- if (se > 0) 2 * stats::pt(-abs(estimate / se), n - 1L) else NA_real_
  data.frame(n_mpa = n, estimate = estimate, se = se,
             ci_low = estimate - critical * se, ci_high = estimate + critical * se, p = p)
}

review_signflip <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L || length(x) > 20L) return(NA_real_)
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), length(x))))
  mean(abs(drop(signs %*% x)) >= abs(sum(x)) - 1e-12)
}

review_reference_analysis_files <- function() {
  c(harmonized_response_ratios = "data/harmonized/harmonized_response_ratios.csv",
    harmonized_landsat_rr = "data/harmonized/harmonized_landsat_rr.csv",
    primary_results = "data/sumstats_final.csv")
}

review_design_input <- function(input, primary, taxa_short) {
  input$Taxa <- unname(taxa_short[as.character(input$y)])
  primary_key <- dplyr::distinct(primary, MPA, Taxa, Source, Resp, AnalysisType, Model)
  input <- input %>% dplyr::filter(is.finite(lnDiff), resp %in% c("Bio", "Den"), !is.na(Taxa)) %>%
    dplyr::inner_join(primary_key,
      by = c("CA_MPA_Name_Short" = "MPA", "Taxa", "source" = "Source", "resp" = "Resp")) %>%
    dplyr::mutate(Region = ifelse(Location == "C", "Islands", "Mainland"),
      series_id = paste(CA_MPA_Name_Short, Taxa, source, resp, sep = "|"))
  observed <- input %>% dplyr::distinct(MPA = CA_MPA_Name_Short, Taxa, Source = source,
                                        Resp = resp, AnalysisType, Model)
  if (nrow(dplyr::anti_join(primary_key, observed, by = names(primary_key))))
    stop("A current primary key has no finite harmonized annual history.")
  if (anyDuplicated(input[c("series_id", "year")])) stop("Duplicate annual primary-series keys.")
  if (anyNA(input$year) || any(!input$BA %in% c("Before", "After")))
    stop("Missing year or invalid before/after classification in primary histories.")
  input
}

review_reference_key_manifest <- function(input) {
  input %>% dplyr::group_by(series_id, MPA = CA_MPA_Name_Short, Taxa, Source = source,
                           Resp = resp, Region, AnalysisType, Model) %>%
    dplyr::summarise(n_years = dplyr::n(),
      observed_years = paste(sort(year), collapse = ","),
      before_years = paste(sort(year[BA == "Before"]), collapse = ","),
      after_years = paste(sort(year[BA == "After"]), collapse = ","), .groups = "drop") %>%
    dplyr::arrange(series_id)
}

review_read_reference_registry <- function(input, root = here::here()) {
  rebuild <- " Rebuild explicitly with code/R/export_reference_registry.R and the source data workspace."
  fail <- function(message) stop(paste0("Stale or invalid reference registry: ", message, rebuild), call. = FALSE)
  registry_files <- c(reference_sites = "data/harmonized/reference_sites.csv",
                      reference_primary_keys = "data/harmonized/reference_primary_keys.csv")
  provenance_file <- file.path(root, "data/harmonized/mapping_provenance.csv")
  if (!file.exists(provenance_file)) fail("mapping_provenance.csv is missing.")
  provenance <- utils::read.csv(provenance_file, stringsAsFactors = FALSE)
  if (!all(c("role", "kind", "path", "md5") %in% names(provenance)) || anyDuplicated(provenance$role) ||
      anyNA(provenance$md5) || any(!grepl("^[0-9a-f]{32}$", provenance$md5)))
    fail("provenance schema or checksums are invalid.")
  checked <- c(review_reference_analysis_files(), registry_files)
  for (role in names(checked)) {
    row <- provenance[provenance$role == role, ]
    kind <- if (role %in% names(registry_files)) "derived_registry" else "analysis_input"
    if (nrow(row) != 1L || row$kind != kind || row$path != checked[[role]])
      fail(paste0("missing or incorrect provenance for ", role, "."))
    path <- file.path(root, checked[[role]])
    if (!file.exists(path) || !identical(unname(tools::md5sum(path)), row$md5))
      fail(paste0("input/footprint MD5 mismatch for ", checked[[role]], "."))
  }
  keys <- utils::read.csv(file.path(root, registry_files[["reference_primary_keys"]]), stringsAsFactors = FALSE)
  current <- as.data.frame(review_reference_key_manifest(input))
  current <- current[order(current$series_id), ]
  if (!identical(names(keys), names(current)) || anyDuplicated(keys$series_id) ||
      !isTRUE(all.equal(keys[order(keys$series_id), ], current, check.attributes = FALSE)))
    fail("primary key, model, before/after classification or calendar coverage changed.")
  sites <- utils::read.csv(file.path(root, registry_files[["reference_sites"]]), stringsAsFactors = FALSE)
  required <- c("series_id", "MPA", "Taxa", "Source", "Resp", "site_id", "reference_id",
    "reference_years", "n_reference_years", "included", "excluded", "identity_verified",
    "exclusion_reason", "contribution_audit", "mapping_file", "mapping_md5")
  if (!all(required %in% names(sites)) || anyNA(sites$series_id) ||
      anyNA(sites$reference_id) || any(!nzchar(sites$reference_id)) ||
      anyDuplicated(sites[c("series_id", "reference_id")]) ||
      !setequal(unique(sites$series_id), keys$series_id)) fail("site registry does not cover every primary series exactly.")
  meta <- sites %>% dplyr::distinct(series_id, MPA, Taxa, Source, Resp)
  if (nrow(dplyr::anti_join(meta, keys, by = names(meta)))) fail("site registry metadata disagrees with primary keys.")
  if (!all(vapply(sites[c("included", "excluded", "identity_verified")], is.logical, logical(1))) ||
      anyNA(sites[c("included", "excluded", "identity_verified", "reference_years", "n_reference_years")]))
    fail("invalid inclusion flags or reference-year coverage.")
  years_by_series <- split(input$year, input$series_id)
  valid_years <- vapply(seq_len(nrow(sites)), function(i) {
    raw <- sites$reference_years[i]
    if (nzchar(raw) && !grepl("^[0-9]{4}(,[0-9]{4})*$", raw)) return(FALSE)
    years <- if (nzchar(raw)) as.integer(strsplit(raw, ",", fixed = TRUE)[[1]]) else integer()
    length(years) == sites$n_reference_years[i] && !anyDuplicated(years) &&
      all(years %in% years_by_series[[sites$series_id[i]]]) &&
      sites$included[i] == (!sites$excluded[i] && length(years) > 0L)
  }, logical(1))
  if (!all(valid_years)) fail("reference survey footprint is inconsistent with current calendar coverage.")
  list(reference_sites = sites, mapping_provenance = provenance, primary_keys = keys)
}

review_gls_fit_audit <- function(fit) {
  phi <- if (is.null(fit$modelStruct$corStruct)) 0 else
    as.numeric(stats::coef(fit$modelStruct$corStruct, unconstrained = FALSE))
  apvar <- fit$apVar
  apvar_not_applicable <- is.null(apvar) && is.null(fit$modelStruct$corStruct) && is.null(fit$modelStruct$varStruct)
  apvar_ok <- if (apvar_not_applicable) NA else is.matrix(apvar) && all(is.finite(apvar)) &&
    min(eigen(apvar, symmetric = TRUE, only.values = TRUE)$values) > 0
  data.frame(ar1_phi = phi, correlation_boundary = abs(phi) > 0.98,
    apvar_ok = apvar_ok,
    apvar_status = if (apvar_not_applicable) "not applicable: no correlation/variance-function parameters" else if (apvar_ok) "positive definite" else if (is.character(apvar))
      paste(apvar, collapse = "; ") else if (is.null(apvar)) "unavailable" else "not positive definite")
}

review_fit_gls <- function(d, formula, correlation = "AR1") {
  stopifnot(correlation %in% c("AR1", "IID"))
  d <- d[order(d$year), , drop = FALSE]
  if (anyDuplicated(d$year)) stop("Duplicate annual observations in GLS series.")
  starts <- if (correlation == "AR1") c(-0.75, -0.25, 0, 0.5, 0.9) else NA_real_
  configurations <- expand.grid(start_phi = starts,
    optimizer = c("nlminb", "optim_BFGS"), stringsAsFactors = FALSE)
  fits <- vector("list", nrow(configurations)); attempts <- vector("list", nrow(configurations))
  for (i in seq_len(nrow(configurations))) {
    cfg <- configurations[i, ]; warnings <- character()
    optimizer <- if (cfg$optimizer == "nlminb") "nlminb" else "optim"
    method <- sub("^optim_", "", cfg$optimizer)
    fit <- tryCatch(withCallingHandlers(
      nlme::gls(formula, data = d, method = "REML",
        correlation = if (correlation == "AR1") nlme::corAR1(cfg$start_phi, form = ~ year) else NULL,
        control = nlme::glsControl(msMaxIter = 200, opt = optimizer,
                                  optimMethod = if (optimizer == "optim") method else "BFGS")),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) e)
    error <- if (inherits(fit, "error")) conditionMessage(fit) else ""
    if (inherits(fit, "error") && !nzchar(error)) error <- "Optimizer failed with an empty error message"
    valid <- !inherits(fit, "error") && all(is.finite(stats::coef(fit))) &&
      is.finite(as.numeric(stats::logLik(fit))) && is.finite(fit$sigma) && fit$sigma > 0 &&
      all(is.finite(stats::vcov(fit))) && min(eigen(stats::vcov(fit), symmetric = TRUE, only.values = TRUE)$values) > 0
    if (!inherits(fit, "error") && !valid) error <- "Nonfinite likelihood/coefficients or invalid coefficient covariance"
    audit <- if (valid) review_gls_fit_audit(fit) else data.frame(
      ar1_phi = NA_real_, correlation_boundary = NA, apvar_ok = NA, apvar_status = "fit failed")
    attempts[[i]] <- cbind(attempt_id = i, cfg, success = valid,
      log_lik = if (valid) as.numeric(stats::logLik(fit)) else NA_real_, audit,
      error = error, warnings = paste(unique(warnings), collapse = "; "), selected = FALSE)
    if (valid) fits[[i]] <- fit
  }
  attempts <- do.call(rbind, attempts)
  good <- which(attempts$success)
  if (!length(good)) return(list(model = NULL,
    error = paste(unique(attempts$error), collapse = "; "), warnings = unique(attempts$warnings[nzchar(attempts$warnings)]),
    attempts = attempts, selection = data.frame(n_attempts = nrow(attempts), n_successful = 0L,
      selected_optimizer = NA_character_, selected_start_phi = NA_real_, best_log_lik = NA_real_,
      optimizer_loglik_gap = NA_real_, near_best_coef_range = NA_real_, optimizer_instability = NA)))
  best <- good[which.max(attempts$log_lik[good])]
  attempts$selected[best] <- TRUE
  near <- good[attempts$log_lik[good] >= attempts$log_lik[best] - 1e-4]
  coefficients <- do.call(rbind, lapply(fits[near], coef))
  coef_range <- max(apply(coefficients, 2, function(x) diff(range(x))))
  optimizer_best <- tapply(attempts$log_lik[good], attempts$optimizer[good], max)
  optimizer_gap <- diff(range(optimizer_best))
  selection <- data.frame(n_attempts = nrow(attempts), n_successful = length(good),
    selected_optimizer = attempts$optimizer[best], selected_start_phi = attempts$start_phi[best],
    best_log_lik = attempts$log_lik[best], optimizer_loglik_gap = optimizer_gap,
    near_best_coef_range = coef_range,
    optimizer_instability = length(optimizer_best) < 2L || optimizer_gap > 0.01 || coef_range > 0.05)
  list(model = fits[[best]], error = "", warnings = attempts$warnings[best],
       attempts = attempts, selection = selection)
}

review_gls_contrast <- function(fit, contrast) {
  beta <- stats::coef(fit)
  contrast <- contrast[names(beta)]
  estimate <- sum(contrast * beta)
  se <- sqrt(drop(t(contrast) %*% stats::vcov(fit) %*% contrast))
  df <- fit$dims$N - length(beta)
  data.frame(estimate = estimate, se = se,
             ci_low = estimate - stats::qt(0.975, df) * se,
             ci_high = estimate + stats::qt(0.975, df) * se,
             p = if (se > 0) 2 * stats::pt(-abs(estimate / se), df) else NA_real_)
}

review_gls_diagnostics <- function(fit, d, label) {
  r <- as.numeric(stats::residuals(fit, type = "normalized"))
  adjacent <- which(diff(d$year) == 1L)
  lag1 <- if (length(adjacent) >= 3L && stats::sd(r[adjacent]) > 0 && stats::sd(r[adjacent + 1L]) > 0)
    stats::cor(r[adjacent], r[adjacent + 1L]) else NA_real_
  lag_at <- function(lag) {
    idx <- which(outer(d$year, d$year, "-") == lag, arr.ind = TRUE)
    if (nrow(idx) < 4L || stats::sd(r[idx[, 1]]) == 0 || stats::sd(r[idx[, 2]]) == 0) return(NA_real_)
    stats::cor(r[idx[, 1]], r[idx[, 2]])
  }
  rp <- as.numeric(stats::residuals(fit, type = "pearson"))
  before <- if ("BA" %in% names(d)) rp[d$BA == "Before"] else numeric()
  after <- if ("BA" %in% names(d)) rp[d$BA == "After"] else numeric()
  variance_ratio <- if (length(before) >= 3L && length(after) >= 3L && stats::var(before) > 0)
    stats::var(after) / stats::var(before) else NA_real_
  cbind(data.frame(model = label, n_years = nrow(d), year_min = min(d$year), year_max = max(d$year),
             residual_lag1 = lag1,
             residual_lag1_n_pairs = length(adjacent), residual_lag2 = lag_at(2L), residual_lag3 = lag_at(3L),
             n_abs_normalized_residual_gt3 = sum(abs(r) > 3),
             after_before_residual_variance_ratio = variance_ratio,
             shapiro_p = if (length(r) >= 3 && length(r) <= 5000 && stats::sd(r) > 0)
               stats::shapiro.test(r)$p.value else NA_real_,
             residual_fitted_undefined = stats::sd(stats::fitted(fit)) == 0,
             abs_residual_fitted_r = if (stats::sd(abs(r)) > 0 && stats::sd(stats::fitted(fit)) > 0)
               stats::cor(abs(r), stats::fitted(fit)) else NA_real_), review_gls_fit_audit(fit))
}

review_reference_components <- function(mpa, reference_ids) {
  stopifnot(length(mpa) == length(reference_ids), !anyDuplicated(mpa))
  parent <- stats::setNames(mpa, mpa)
  root <- function(x) { while (parent[[x]] != x) x <- parent[[x]]; x }
  ids <- strsplit(reference_ids, ";", fixed = TRUE)
  references <- unique(unlist(ids))
  for (id in references[nzchar(references)]) {
    members <- mpa[vapply(ids, function(x) id %in% x, logical(1))]
    roots <- vapply(members, root, character(1))
    parent[roots] <- sort(roots)[1]
  }
  paste0("RG|", vapply(mpa, root, character(1)))
}

review_group_summary <- function(x, groups, mapping_complete = TRUE) {
  keep <- is.finite(x); x <- x[keep]; groups <- groups[keep]
  naive <- review_summary(x)
  n <- length(x); sizes <- table(groups); g <- length(sizes)
  out <- data.frame(n_mpa = n, estimate = naive$estimate,
    se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_,
    iid_se = naive$se, iid_ci_low = naive$ci_low, iid_ci_high = naive$ci_high, iid_p = naive$p,
    n_reference_groups = g, reference_group_df = if (g) g - 1L else NA_integer_,
    largest_reference_group_n = if (g) max(sizes) else 0L,
    reference_group_size_ratio = if (g) max(sizes) / min(sizes) else NA_real_,
    few_reference_groups = g < 10L, reference_mapping_complete = all(mapping_complete),
    uncertainty_method = "delete-reference-group jackknife with t(G-1) interval",
    uncertainty_status = if (g < 3L) "fewer than 3 groups: interval unavailable" else if (g < 10L)
      "few groups: approximate descriptive uncertainty" else "approximate group-aware uncertainty")
  if (!out$reference_mapping_complete) out$uncertainty_status <-
    paste0(out$uncertainty_status, "; reference identity/overlap audit incomplete")
  if (g < 3L) return(out)
  deleted <- vapply(names(sizes), function(id) mean(x[groups != id]), numeric(1))
  se <- sqrt((g - 1) / g * sum((deleted - mean(deleted))^2))
  out$se <- se
  out$ci_low <- out$estimate - stats::qt(0.975, g - 1L) * se
  out$ci_high <- out$estimate + stats::qt(0.975, g - 1L) * se
  out$p <- if (se > 0) 2 * stats::pt(-abs(out$estimate / se), g - 1L) else NA_real_
  out
}

review_export_figure <- function(p, stem, width = 170, height = 125) {
  ggplot2::ggsave(here::here("plots", paste0(stem, ".pdf")), p,
                  width = width, height = height, units = "mm",
                  device = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf")
  ggplot2::ggsave(here::here("plots", paste0(stem, ".png")), p,
                  width = width, height = height, units = "mm", dpi = 600, bg = "white")
}
