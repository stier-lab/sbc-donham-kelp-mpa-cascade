# Arithmetic, coverage and integration checks for the coauthor-requested analyses.
local({
  if (!requireNamespace("here", quietly = TRUE)) stop("Package here is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package dplyr is required.")
  `%>%` <- dplyr::`%>%`
  source(here::here("code", "R", "00c_analysis_constants.R"), local = TRUE)
  source(here::here("code", "R", "review_concern_utils.R"), local = TRUE)
  csv <- function(stem) read.csv(here::here("tables", paste0(stem, ".csv")), stringsAsFactors = FALSE)
  near <- function(a, b, label, tolerance = 1e-9) {
    a <- as.numeric(a); b <- as.numeric(b)
    if (length(a) != length(b) || !identical(is.na(a), is.na(b)) ||
        !identical(is.finite(a), is.finite(b))) stop(label, " disagree.")
    finite <- is.finite(a)
    infinite <- !is.na(a) & !finite
    if (any(abs(a[finite] - b[finite]) > tolerance * pmax(1, abs(a[finite]), abs(b[finite]))) ||
        any(a[infinite] != b[infinite])) stop(label, " disagree.")
  }
  require_columns <- function(x, columns, label) {
    missing <- setdiff(columns, names(x))
    if (length(missing)) stop(label, " missing columns: ", paste(missing, collapse = ", "))
  }
  key <- function(x, columns) do.call(paste, c(x[columns], sep = "|"))
  same_keys <- function(a, b, columns, label) {
    require_columns(a, columns, label); require_columns(b, columns, label)
    if (anyNA(a[columns]) || anyNA(b[columns]) || anyDuplicated(a[columns]) ||
        anyDuplicated(b[columns]) || !setequal(key(a, columns), key(b, columns)))
      stop(label, " have missing, duplicate or unequal keys.")
    b[match(key(a, columns), key(b, columns)), , drop = FALSE]
  }
  select_rows <- function(x, row, columns) {
    keep <- rep(TRUE, nrow(x))
    for (column in columns) keep <- keep & x[[column]] == row[[column]]
    x[keep, , drop = FALSE]
  }
  same_table <- function(a, b, columns, label) {
    b <- same_keys(a, b, columns, label)
    if (!identical(names(a), names(b)) || !isTRUE(all.equal(a, b, check.attributes = FALSE, tolerance = 1e-9)))
      stop(label, " disagree.")
  }
  # Independent graph traversal checks the production union-find implementation.
  connected_groups <- function(mpas, ids) {
    references <- strsplit(ids, ";", fixed = TRUE)
    group <- rep(NA_character_, length(mpas))
    for (i in seq_along(mpas)) {
      if (!is.na(group[i])) next
      members <- i
      repeat {
        pooled <- unique(unlist(references[members]))
        expanded <- which(vapply(references, function(x) length(intersect(x, pooled)) > 0L, logical(1)))
        if (setequal(members, expanded)) break
        members <- expanded
      }
      group[members] <- paste0("RG|", sort(mpas[members])[1])
    }
    group
  }
  group_stats <- function(x, groups) {
    n <- length(x); g <- length(unique(groups)); estimate <- if (n) mean(x) else NA_real_
    se <- low <- high <- p <- NA_real_
    if (g >= 3L) {
      deleted <- vapply(unique(groups), function(id) mean(x[groups != id]), numeric(1))
      se <- sqrt((g - 1L) / g * sum((deleted - mean(deleted))^2))
      low <- estimate - stats::qt(.975, g - 1L) * se
      high <- estimate + stats::qt(.975, g - 1L) * se
      if (se > 0) p <- 2 * stats::pt(-abs(estimate / se), g - 1L)
    }
    c(n_mpa = n, n_reference_groups = g, estimate = estimate, se = se, ci_low = low, ci_high = high, p = p)
  }
  primary <- read.csv(here::here("data", "sumstats_final.csv"))
  meta <- csv("table_02_meta_analysis")
  near(meta$RR, exp(meta$Estimate), "Full-precision response ratio")
  near(meta$RR_lower, exp(meta$CI_lower), "Full-precision lower interval")
  near(meta$RR_upper, exp(meta$CI_upper), "Full-precision upper interval")
  old <- read.csv(here::here("data", "harmonized", "harmonized_landsat_rr.csv"))
  annual <- read.csv(here::here("data", "harmonized", "harmonized_landsat_annual_index.csv"))
  stopifnot(!anyDuplicated(annual[c("CA_MPA_Name_Short", "year", "status")]))
  for (status in c("mpa", "reference")) {
    z <- annual[annual$status == status, ]
    z <- merge(old[c("CA_MPA_Name_Short", "year", status)], z, by = c("CA_MPA_Name_Short", "year"))
    stopifnot(nrow(z) == nrow(old))
    near(z[[status]], z$index + z$log_offset, "Uncorrected Landsat export")
  }
  hp <- csv("table_s_historical_heatwave_pairs")
  hc <- csv("table_s_historical_heatwave_coverage")
  stopifnot(all(hp$MPA %in% primary$MPA[primary$Source == "Landsat"]),
            !anyDuplicated(hp[c("specification", "MPA", "metric")]),
            sum(hc$eligible) == nrow(hp))
  near(hp$historical_log_advantage, log(hp$historical_mpa / hp$historical_reference), "Historical advantage")
  near(hp$modern_log_advantage, log(hp$modern_mpa / hp$modern_reference), "Modern advantage")
  near(hp$change_log_advantage, hp$modern_log_advantage - hp$historical_log_advantage, "Within-pair event contrast")
  hs <- csv("table_s_historical_heatwave_summary")
  for (i in seq_len(nrow(hs))) {
    x <- hp[hp$specification == hs$specification[i] & hp$metric == hs$metric[i], ]
    if (hs$Region[i] != "All") x <- x[x$Region == hs$Region[i], ]
    stopifnot(hs$n_mpa[i] == nrow(x))
    near(hs$estimate[i], mean(x$change_log_advantage), "Historical summary")
  }
  for (family in unique(hs$test_family)) {
    ix <- hs$test_family == family
    near(hs$p_fdr[ix], p.adjust(hs$p[ix], "BH"), "Historical multiplicity family")
  }
  common <- csv("table_s_historical_heatwave_common_cohort")
  stopifnot(!anyDuplicated(common[c("sensitivity", "specification", "metric")]))
  for (i in seq_len(nrow(common))) {
    ids <- Reduce(intersect, lapply(c("primary", common$sensitivity[i]), function(spec)
      hp$MPA[hp$specification == spec & hp$metric == common$metric[i]]))
    x <- hp[hp$specification == common$specification[i] & hp$metric == common$metric[i] & hp$MPA %in% ids, ]
    stopifnot(common$n_mpa[i] == length(ids))
    near(common$estimate[i], mean(x$change_log_advantage), "Fixed-cohort window comparison")
  }
  de <- csv("table_s_design_matched_estimates")
  dp <- csv("table_s_design_matched_pairs")
  require_columns(de, c("series_id", "MPA", "Taxa", "Source", "Resp", "correlation", "design", "analysis_window",
    "estimate", "se", "p", "fit_status", "fit_error", "reference_ids", "reference_mapping_complete",
    "n_fit_years", "fit_year_min", "fit_year_max", "n_before", "n_after", "n_attempts", "n_successful",
    "best_log_lik", "selected_optimizer", "selected_start_phi", "optimizer_loglik_gap", "near_best_coef_range",
    "optimizer_instability", "ar1_phi", "correlation_boundary", "apvar_ok", "apvar_status"), "Source-estimate audit")
  keys <- c("MPA", "Taxa", "Source", "Resp")
  stopifnot(nrow(merge(unique(de[keys]), unique(primary[keys]), by = keys)) == nrow(unique(de[keys])),
            all(de$n_after >= 3L),
            all(de$n_before[de$design == "BACIPS_step"] >= 3L),
            !anyDuplicated(dp[c("MPA", "Taxa", "Resp", "correlation", "source_subset", "design")]))
  rr <- review_design_input(dplyr::bind_rows(
    read.csv(here::here("data", "harmonized", "harmonized_response_ratios.csv")), old), primary, RESILIENCE_TAXA_SHORT)
  registry <- review_read_reference_registry(rr)
  same_table(csv("table_s_design_reference_sites"), registry$reference_sites,
             c("series_id", "reference_id"), "Exported reference registry")
  same_table(csv("table_s_design_reference_provenance"), registry$mapping_provenance,
             "role", "Exported reference provenance")
  manifest <- utils::read.table(here::here("data", "harmonized", "checksums.md5"),
    col.names = c("md5", "file"), stringsAsFactors = FALSE)
  for (file in c("reference_sites.csv", "reference_primary_keys.csv", "mapping_provenance.csv")) {
    entry <- manifest[manifest$file == file, ]
    if (nrow(entry) != 1L || entry$md5 != unname(tools::md5sum(here::here("data", "harmonized", file))))
      stop("Harmonized checksum manifest disagrees with registry file: ", file)
  }

  dc <- csv("table_s_design_matched_coverage")
  same_keys(dc, registry$primary_keys, "series_id", "Design coverage")
  expected_annual <- common_expected <- list()
  for (sid in registry$primary_keys$series_id) {
    d <- rr[rr$series_id == sid, ]; row <- dc[dc$series_id == sid, ]
    before <- sum(d$BA == "Before"); after <- sum(d$BA == "After")
    stopifnot(row$n_before == before, row$n_after == after,
              row$matched_eligible == (before >= 3L && after >= 3L))
    common_years <- sort(d$year[d$BA == "After" & d$year >= 2013L & d$year <= 2023L])
    if (d$source[1] == "PISCO" && after >= 3L) common_expected[[sid]] <- data.frame(
      series_id = sid, n_common_years = length(common_years), expected_years = 11L,
      common_calendar_start = 2013L, common_calendar_end = 2023L,
      observed_years = paste(common_years, collapse = ","), eligible = length(common_years) >= 3L,
      complete_calendar = identical(common_years, 2013:2023))
    if (after < 3L) next
    designs <- c("CI_after", if (before >= 3L) c("BACIPS_step", "Before_level"))
    expected_annual[[sid]] <- data.frame(series_id = sid, design = designs, analysis_window = "Full history")
    if (d$source[1] == "PISCO" && length(common_years) >= 3L) expected_annual[[paste0(sid, " common")]] <-
      data.frame(series_id = sid, design = "CI_after", analysis_window = "PISCO 2013-2023")
  }
  ac <- de[de$correlation == "Annual means", ]
  same_keys(ac, dplyr::bind_rows(expected_annual), c("series_id", "design", "analysis_window"), "Arithmetic primary cohorts")
  stopifnot(all(ac$fit_status == "ok"), all(is.na(ac$se)), all(is.na(ac$p)))
  arithmetic <- de[(de$correlation == "Annual means" | de$correlation == "IID") & de$fit_status == "ok", ]
  for (i in seq_len(nrow(arithmetic))) {
    row <- arithmetic[i, ]; d <- rr[rr$series_id == row$series_id, ]
    if (row$analysis_window == "PISCO 2013-2023") d <- d[d$year >= 2013L & d$year <= 2023L, ]
    expected <- switch(row$design,
      CI_after = mean(d$lnDiff[d$BA == "After"]),
      BACIPS_step = mean(d$lnDiff[d$BA == "After"]) - mean(d$lnDiff[d$BA == "Before"]),
      Before_level = mean(d$lnDiff[d$BA == "Before"]), Before_trend = NULL)
    if (!is.null(expected)) near(row$estimate, expected, "Annual-mean/IID arithmetic")
  }
  cc <- csv("table_s_design_pisco_common_calendar_coverage")
  expected <- same_keys(cc, dplyr::bind_rows(common_expected), "series_id", "PISCO common-calendar coverage")
  for (column in setdiff(names(expected), "series_id"))
    if (!identical(cc[[column]], expected[[column]])) stop("Common-calendar coverage disagrees: ", column)

  # Verify every source's reference set against the vendored footprint, including failed GLS rows.
  for (i in seq_len(nrow(de))) {
    row <- de[i, ]; d <- rr[rr$series_id == row$series_id, ]
    if (row$analysis_window == "PISCO 2013-2023") d <- d[d$year >= 2013L & d$year <= 2023L, ]
    if (row$design == "CI_after") d <- d[d$BA == "After", ]
    if (grepl("^Before", row$design)) d <- d[d$BA == "Before", ]
    sites <- registry$reference_sites[registry$reference_sites$series_id == row$series_id & registry$reference_sites$included, ]
    if (nrow(sites)) sites <- sites[vapply(strsplit(sites$reference_years, ",", fixed = TRUE), function(years)
      any(as.integer(years) %in% d$year), logical(1)), ]
    ids <- if (nrow(sites)) paste(sort(unique(sites$reference_id)), collapse = ";") else paste0("unmapped|", row$series_id)
    stopifnot(row$reference_ids == ids, row$reference_mapping_complete == (nrow(sites) > 0L && all(sites$identity_verified)),
      row$n_fit_years == nrow(d), row$fit_year_min == min(d$year), row$fit_year_max == max(d$year))
  }

  check_pairs <- function(actual, window, matched) {
    source <- de[de$analysis_window == window & de$fit_status == "ok" & de$design %in% c("CI_after", "BACIPS_step"), ]
    if (matched) source <- source %>% dplyr::group_by(series_id, correlation) %>%
      dplyr::filter(dplyr::n_distinct(design) == 2L) %>% dplyr::ungroup() else source <- source[source$design == "CI_after", ]
    subsets <- if (window == "PISCO 2013-2023") "PISCO only" else c("All programmes", "PISCO only", "Divers only", "Landsat only")
    expected <- dplyr::bind_rows(lapply(subsets, function(subset) {
      d <- switch(subset, "All programmes" = source, "PISCO only" = source[source$Source == "PISCO", ],
        "Divers only" = source[source$Source != "Landsat", ], "Landsat only" = source[source$Source == "Landsat", ])
      d %>% dplyr::group_by(MPA, Taxa, Resp, Region, correlation, design) %>%
        dplyr::summarise(estimate = mean(estimate), n_sources = dplyr::n_distinct(Source),
          Sources = paste(sort(unique(Source)), collapse = ";"),
          reference_ids = paste(sort(unique(unlist(strsplit(reference_ids, ";", fixed = TRUE)))), collapse = ";"),
          reference_mapping_complete = all(reference_mapping_complete),
          n_boundary_sources = sum(correlation_boundary, na.rm = TRUE),
          n_apvar_problem_sources = sum(!apvar_ok, na.rm = TRUE),
          n_unstable_sources = sum(optimizer_instability, na.rm = TRUE), .groups = "drop") %>%
        dplyr::mutate(source_subset = subset)
    }))
    expected <- same_keys(actual, expected, c("MPA", "Taxa", "Resp", "Region", "correlation", "design", "source_subset"), "Equal-MPA source averages")
    near(actual$estimate, expected$estimate, "Equal-MPA source averages")
    for (column in c("n_boundary_sources", "n_apvar_problem_sources", "n_unstable_sources"))
      near(actual[[column]], expected[[column]], paste("Source audit propagation", column))
    stopifnot(identical(actual$n_sources, expected$n_sources), identical(actual$Sources, expected$Sources),
      identical(actual$reference_ids, expected$reference_ids), identical(actual$reference_mapping_complete, expected$reference_mapping_complete))
  }
  all_ci <- csv("table_s_design_all_ci_pairs"); common_ci <- csv("table_s_design_pisco_common_calendar_pairs")
  check_pairs(dp, "Full history", TRUE); check_pairs(all_ci, "Full history", FALSE)
  check_pairs(common_ci, "PISCO 2013-2023", FALSE)
  stopifnot(all(common_ci$cohort == "PISCO 2013-2023 after-only"), all(common_ci$source_subset == "PISCO only"))
  analysis_pairs <- dplyr::bind_rows(dplyr::mutate(dp, cohort = "Matched before/after histories"), all_ci, common_ci)
  view_keys <- c("Taxa", "Resp", "correlation", "source_subset", "design", "cohort", "Region_scope")
  views <- dplyr::bind_rows(dplyr::mutate(analysis_pairs, Region_scope = "All"),
                           dplyr::mutate(analysis_pairs, Region_scope = Region)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(view_keys))) %>%
    dplyr::group_modify(~ dplyr::mutate(.x, reference_group_id = connected_groups(MPA, reference_ids))) %>% dplyr::ungroup()
  rg <- csv("table_s_design_reference_groups")
  expected <- same_keys(rg, views, c(view_keys, "MPA"), "Reference-group membership")
  stopifnot(identical(rg$reference_group_id, expected$reference_group_id), identical(rg$reference_ids, expected$reference_ids))
  sizes <- table(key(rg, c(view_keys, "reference_group_id")))
  near(rg$n_mpa_in_reference_group, as.numeric(sizes[key(rg, c(view_keys, "reference_group_id"))]), "Reference-group sizes")
  gs <- csv("table_s_design_group_summary")
  summary_keys <- c(setdiff(view_keys, "Region_scope"), "Region")
  require_columns(gs, c(summary_keys, "n_mpa", "estimate", "se", "ci_low", "ci_high", "p", "p_group_approx", "p_fdr",
    "iid_se", "iid_ci_low", "iid_ci_high", "iid_p", "n_reference_groups", "reference_group_df", "few_reference_groups",
    "reference_mapping_complete", "uncertainty_method", "uncertainty_status", "fdr_family_n",
    "n_boundary_sources", "n_apvar_problem_sources", "n_unstable_sources"), "Group summary")
  expected_keys <- dplyr::distinct(dplyr::rename(views[view_keys], Region = Region_scope))
  same_keys(gs, expected_keys, summary_keys, "Group-summary coverage")
  for (i in seq_len(nrow(gs))) {
    row <- gs[i, ]; row$Region_scope <- row$Region; d <- select_rows(views, row, view_keys)
    values <- group_stats(d$estimate, d$reference_group_id)
    for (column in setdiff(names(values), "p")) near(row[[column]], values[[column]], paste("Group summary", column))
    near(row$p_group_approx, values[["p"]], "Descriptive group p audit")
    near(row$p, if (row$correlation == "Annual means") values[["p"]] else NA_real_, "Primary/descriptive p separation")
    iid <- review_summary(d$estimate)
    for (column in c("se", "ci_low", "ci_high", "p")) near(row[[paste0("iid_", column)]], iid[[column]], paste("IID audit", column))
    stopifnot(row$few_reference_groups == (values[["n_reference_groups"]] < 10L),
      row$reference_group_df == values[["n_reference_groups"]] - 1L,
      row$reference_mapping_complete == all(d$reference_mapping_complete))
    for (column in c("n_boundary_sources", "n_apvar_problem_sources", "n_unstable_sources"))
      near(row[[column]], sum(d[[column]]), paste("Group audit propagation", column))
    if (values[["n_reference_groups"]] < 3L) stopifnot(is.na(row$se), is.na(row$ci_low), is.na(row$p))
  }
  for (family in split(gs, key(gs, c("correlation", "source_subset", "design", "cohort", "Region")))) {
    near(family$p_fdr, stats::p.adjust(family$p, "BH"), "Design multiplicity family")
    stopifnot(all(family$fdr_family_n == sum(is.finite(family$p))))
  }
  for (spec in list(c("table_s_design_matched_summary", "Matched before/after histories"),
    c("table_s_design_all_ci_summary", "All primary histories with >=3 after years"),
    c("table_s_design_pisco_common_calendar_summary", "PISCO 2013-2023 after-only")))
    same_table(csv(spec[1]), gs[gs$cohort == spec[2], ], summary_keys, "Cohort-specific summary view")

  annual_views <- views[views$correlation == "Annual means", ]
  for (kind in c("mpa", "reference_group")) {
    influence <- csv(paste0("table_s_design_leave_", kind, "_out"))
    require_columns(influence, c(view_keys, "omitted_unit", "omitted_mpas", "n_removed_mpa", "full_n_mpa",
      "full_n_reference_groups", "full_estimate", "estimate_change", "sign_changed", "n_mpa", "n_reference_groups",
      "estimate", "se", "ci_low", "ci_high", "reference_group_df", "few_reference_groups", "reference_mapping_complete"), "Influence audit")
    expected_keys <- dplyr::distinct(annual_views[c(view_keys, if (kind == "mpa") "MPA" else "reference_group_id")])
    names(expected_keys)[ncol(expected_keys)] <- "omitted_unit"
    same_keys(influence, expected_keys, c(view_keys, "omitted_unit"), "Influence coverage")
    stopifnot(!any(c("p", "iid_p") %in% names(influence)))
    for (i in seq_len(nrow(influence))) {
      row <- influence[i, ]; d <- select_rows(annual_views, row, view_keys)
      omitted <- if (kind == "mpa") d$MPA == row$omitted_unit else d$reference_group_id == row$omitted_unit
      remaining <- d[!omitted, ]; values <- group_stats(remaining$estimate, remaining$reference_group_id)
      stopifnot(row$omitted_mpas == paste(sort(d$MPA[omitted]), collapse = ";"), row$n_removed_mpa == sum(omitted),
        row$full_n_mpa == nrow(d), row$full_n_reference_groups == length(unique(d$reference_group_id)))
      near(row$reference_group_df, if (values[["n_reference_groups"]]) values[["n_reference_groups"]] - 1L else NA_real_, "Influence group df")
      stopifnot(row$few_reference_groups == (values[["n_reference_groups"]] < 10L),
        row$reference_mapping_complete == all(remaining$reference_mapping_complete))
      for (column in setdiff(names(values), "p")) near(row[[column]], values[[column]], paste("Influence", column))
      near(row$full_estimate, mean(d$estimate), "Influence full mean")
      near(row$estimate_change, values[["estimate"]] - mean(d$estimate), "Influence estimate change")
      expected_sign <- if (is.finite(values[["estimate"]])) sign(values[["estimate"]]) != sign(mean(d$estimate)) else NA
      stopifnot(identical(row$sign_changed, expected_sign))
    }
  }

  geography <- csv("table_s_design_geography")
  geo_keys <- c("Taxa", "Resp", "correlation", "source_subset", "design", "cohort")
  same_keys(geography, unique(analysis_pairs[geo_keys]), geo_keys, "Geographic cohort coverage")
  for (i in seq_len(nrow(geography))) {
    row <- geography[i, ]; d <- select_rows(analysis_pairs, row, geo_keys)
    groups <- connected_groups(d$MPA, d$reference_ids)
    island <- d$Region == "Islands"; mainland <- d$Region == "Mainland"
    a <- group_stats(d$estimate[island], groups[island]); b <- group_stats(d$estimate[mainland], groups[mainland])
    cross <- length(intersect(groups[island], groups[mainland]))
    stopifnot(row$n_islands == sum(island), row$n_mainland == sum(mainland),
      row$n_islands_reference_groups == a[["n_reference_groups"]],
      row$n_mainland_reference_groups == b[["n_reference_groups"]], row$cross_region_reference_groups == cross,
      row$few_reference_groups == (a[["n_reference_groups"]] < 10L || b[["n_reference_groups"]] < 10L),
      row$reference_mapping_complete == all(d$reference_mapping_complete))
    difference <- a[["estimate"]] - b[["estimate"]]
    near(row$islands_minus_mainland, difference, "Equal-MPA geographic difference")
    se <- df <- low <- high <- p <- NA_real_
    variance <- a[["se"]]^2 + b[["se"]]^2
    if (is.finite(variance) && variance > 0 && cross == 0L) {
      se <- sqrt(variance)
      df <- variance^2 / (a[["se"]]^4 / (a[["n_reference_groups"]] - 1L) +
        b[["se"]]^4 / (b[["n_reference_groups"]] - 1L))
      low <- difference - stats::qt(.975, df) * se; high <- difference + stats::qt(.975, df) * se
      p <- 2 * stats::pt(-abs(difference / se), df)
    }
    near(row$se, se, "Geographic group SE"); near(row$reference_group_df, df, "Geographic group df")
    near(row$ci_low, low, "Geographic lower interval"); near(row$ci_high, high, "Geographic upper interval")
    near(row$p_group_approx, p, "Geographic descriptive p audit")
    near(row$p, if (row$correlation == "Annual means") p else NA_real_, "Geographic p separation")
  }
  for (family in split(geography, key(geography, c("correlation", "source_subset", "design", "cohort")))) {
    near(family$p_fdr, stats::p.adjust(family$p, "BH"), "Geographic multiplicity family")
    stopifnot(all(family$fdr_family_n == sum(is.finite(family$p))))
  }

  attempts <- csv("table_s_design_gls_fit_attempts")
  fit_keys <- c("series_id", "correlation", "design", "analysis_window")
  require_columns(attempts, c(fit_keys, "attempt_id", "start_phi", "optimizer", "success", "log_lik", "selected",
    "ar1_phi", "correlation_boundary", "apvar_ok", "apvar_status", "error", "warnings"), "GLS attempts")
  gls <- de[de$correlation %in% c("AR1", "IID"), ]
  same_keys(gls, unique(attempts[fit_keys]), fit_keys, "GLS-attempt coverage")
  for (i in seq_len(nrow(gls))) {
    row <- gls[i, ]; a <- select_rows(attempts, row, fit_keys); selected <- a[a$selected, ]
    expected_n <- if (row$correlation == "AR1") 10L else 2L
    stopifnot(nrow(a) == expected_n, row$n_attempts == expected_n, row$n_successful == sum(a$success),
      !anyDuplicated(a$attempt_id), sum(a$selected) == as.integer(any(a$success)),
      row$fit_status == if (any(a$success)) "ok" else "failed",
      all(nzchar(a$error[!a$success])), all(is.finite(a$log_lik[a$success])))
    for (optimizer in c("nlminb", "optim_BFGS")) {
      starts <- a$start_phi[a$optimizer == optimizer]
      if (row$correlation == "AR1") stopifnot(setequal(starts, c(-.75, -.25, 0, .5, .9))) else
        stopifnot(length(starts) == 1L, is.na(starts))
    }
    stopifnot(all(a$correlation_boundary[a$success] == (abs(a$ar1_phi[a$success]) > .98)))
    if (nrow(selected)) {
      near(row$best_log_lik, max(a$log_lik[a$success]), "Maximum-likelihood selection")
      near(selected$log_lik, row$best_log_lik, "Selected likelihood audit")
      near(selected$ar1_phi, row$ar1_phi, "Selected correlation audit")
      near(selected$start_phi, row$selected_start_phi, "Selected starting-correlation audit")
      optimizer_best <- tapply(a$log_lik[a$success], a$optimizer[a$success], max)
      near(row$optimizer_loglik_gap, diff(range(optimizer_best)), "Optimizer likelihood disagreement")
      stopifnot(row$selected_optimizer == selected$optimizer, identical(row$apvar_ok, selected$apvar_ok),
        row$apvar_status == selected$apvar_status, row$correlation_boundary == (abs(row$ar1_phi) > .98),
        row$optimizer_instability == (length(optimizer_best) < 2L || row$optimizer_loglik_gap > .01 || row$near_best_coef_range > .05))
      if (row$correlation == "IID") stopifnot(is.na(row$apvar_ok), row$ar1_phi == 0)
    } else stopifnot(is.na(row$estimate), nzchar(row$fit_error))
  }
  atlas <- csv("table_s_pair_history_atlas_index")
  same_keys(atlas, primary, keys, "Atlas primary-series coverage")
  stopifnot(setequal(unique(atlas$MPA), unique(primary$MPA)),
            nrow(atlas) == nrow(primary), all(file.exists(here::here(atlas$png))))
  required <- c("plots/fig_s_design_matched.pdf", "plots/fig_s_design_matched.png",
                "plots/fig_s_historical_heatwave_bacips.pdf", "plots/fig_s_historical_heatwave_bacips.png",
                "plots/fig_s_pair_history_atlas.pdf")
  stopifnot(all(file.exists(here::here(required))), all(file.info(here::here(required))$size > 0))
  cat("Review-concern checks passed: registry integrity, matched estimands, group uncertainty/influence, common calendar, GLS audits, event arithmetic, multiplicity and MPA coverage.\n")
})
