# Compare MPA-designation estimands on identical annual pair histories.
# Run after 08_effect_sizes.R; source-specific series remain distinct until
# equal-weight MPA summaries, so multiple programmes do not replicate an MPA.

local({
  suppressPackageStartupMessages({ library(here); library(dplyr); library(tidyr); library(ggplot2) })
  `%>%` <- dplyr::`%>%`
  source(here::here("code", "R", "00c_analysis_constants.R"), local = TRUE)
  source(here::here("code", "R", "00b_color_palette.R"), local = TRUE)
  source(here::here("code", "R", "review_concern_utils.R"), local = TRUE)
  cat("Running designation-design and geography sensitivity (script 27)...\n")
  input <- dplyr::bind_rows(
    utils::read.csv(here::here("data", "harmonized", "harmonized_response_ratios.csv")),
    utils::read.csv(here::here("data", "harmonized", "harmonized_landsat_rr.csv")))
  primary <- utils::read.csv(here::here("data", "sumstats_final.csv"))
  input <- review_design_input(input, primary, RESILIENCE_TAXA_SHORT)
  # Routine analyses use only checksum-bound, vendored identities and footprints.
  registry <- review_read_reference_registry(input)
  reference_sites <- registry$reference_sites
  mapping_provenance <- registry$mapping_provenance
  reference_details <- function(sid, years) {
    z <- reference_sites %>% dplyr::filter(series_id == sid, included)
    if (nrow(z)) z <- z[vapply(strsplit(z$reference_years, ",", fixed = TRUE), function(y)
      any(as.integer(y) %in% years), logical(1)), ]
    if (!nrow(z)) return(data.frame(reference_ids = paste0("unmapped|", sid),
      reference_mapping_complete = FALSE, reference_mapping_note = "No mapped reference in analysis years"))
    data.frame(reference_ids = paste(sort(unique(z$reference_id)), collapse = ";"),
      reference_mapping_complete = all(z$identity_verified),
      reference_mapping_note = paste(sort(unique(z$contribution_audit)), collapse = "; "))
  }
  estimates <- diagnostics <- coverage <- fit_attempts <- common_coverage <- list()
  fitted_series <- function(d, meta, correlation, design, n_before, n_after, window = "Full history") {
    formula <- if (design == "BACIPS_step") lnDiff ~ BA else if (design == "Before_trend")
      lnDiff ~ I(year - min(year)) else lnDiff ~ 1
    result <- review_fit_gls(d, formula, correlation)
    failed <- is.null(result$model)
    values <- if (failed) data.frame(estimate = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_) else {
      contrast <- stats::setNames(rep(0, length(stats::coef(result$model))), names(stats::coef(result$model)))
      contrast[if (design == "BACIPS_step") "BAAfter" else if (design == "Before_trend")
        names(contrast)[2] else "(Intercept)"] <- 1
      review_gls_contrast(result$model, contrast)
    }
    key <- paste(meta$series_id, window, correlation, design)
    audit <- if (failed) data.frame(ar1_phi = NA_real_, correlation_boundary = NA,
      apvar_ok = NA, apvar_status = "fit failed") else review_gls_fit_audit(result$model)
    estimates[[key]] <<- cbind(meta, correlation = correlation, design = design, analysis_window = window,
      n_before = n_before, n_after = n_after, n_fit_years = nrow(d), fit_year_min = min(d$year), fit_year_max = max(d$year),
      values, fit_status = if (failed) "failed" else "ok", fit_error = result$error,
      fit_warnings = paste(result$warnings, collapse = "; "), result$selection, audit,
      inference_status = "descriptive GLS sensitivity; conditional model SE/CI/p only",
      reference_details(meta$series_id, d$year))
    fit_attempts[[key]] <<- cbind(meta, correlation = correlation, design = design,
      analysis_window = window, result$attempts)
    if (!failed) diagnostics[[key]] <<- cbind(meta, correlation = correlation, analysis_window = window,
      review_gls_diagnostics(result$model, d, design), result$selection)
  }
  annual_estimate <- function(d, meta, design, value, n_before, n_after, window = "Full history") {
    key <- paste(meta$series_id, window, "Annual means", design)
    estimates[[key]] <<- cbind(meta, correlation = "Annual means", design = design, analysis_window = window,
      n_before = n_before, n_after = n_after, n_fit_years = nrow(d), fit_year_min = min(d$year), fit_year_max = max(d$year),
      estimate = unname(value), se = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_,
      fit_status = "ok", fit_error = "", fit_warnings = "Uncertainty estimated among reference groups; no annual iid inference",
      n_attempts = 0L, n_successful = 0L, selected_optimizer = NA_character_, selected_start_phi = NA_real_,
      best_log_lik = NA_real_, optimizer_loglik_gap = NA_real_, near_best_coef_range = NA_real_,
      optimizer_instability = FALSE, ar1_phi = NA_real_, correlation_boundary = FALSE,
      apvar_ok = NA, apvar_status = "not applicable: arithmetic mean",
      inference_status = "primary equal-MPA annual mean; approximate group-aware uncertainty",
      reference_details(meta$series_id, d$year))
  }
  for (sid in sort(unique(input$series_id))) {
    d <- input[input$series_id == sid, , drop = FALSE]
    d <- d[order(d$year), , drop = FALSE]
    d$BA <- factor(d$BA, levels = c("Before", "After"))
    n_before <- sum(d$BA == "Before"); n_after <- sum(d$BA == "After")
    meta <- data.frame(series_id = sid, MPA = d$CA_MPA_Name_Short[1], Taxa = d$Taxa[1],
                       Source = d$source[1], Resp = d$resp[1], Region = d$Region[1],
                       primary_design = d$AnalysisType[1], primary_model = d$Model[1])
    eligible <- n_before >= 3L && n_after >= 3L
    coverage[[sid]] <- cbind(meta, n_before = n_before, n_after = n_after,
                              year_min = min(d$year), year_max = max(d$year),
                              matched_eligible = eligible,
                              reason = if (eligible) "eligible" else "fewer than 3 paired years in before or after period")
    if (n_after < 3L) next
    # Keep arithmetic period means and equal-MPA weighting as the primary bridge.
    after_mean <- mean(d$lnDiff[d$BA == "After"])
    mean_values <- c(CI_after = after_mean)
    if (eligible) mean_values <- c(mean_values,
      BACIPS_step = after_mean - mean(d$lnDiff[d$BA == "Before"]),
      Before_level = mean(d$lnDiff[d$BA == "Before"]))
    for (design in names(mean_values)) {
      dd <- if (design == "CI_after") d[d$BA == "After", ] else if (design == "Before_level") d[d$BA == "Before", ] else d
      annual_estimate(dd, meta, design, mean_values[design], n_before, n_after)
    }
    for (correlation in c("AR1", "IID")) {
      designs <- c("CI_after", if (eligible) c("BACIPS_step", "Before_level", "Before_trend"))
      for (design in designs) {
        dd <- if (design == "CI_after") d[d$BA == "After", ] else if (grepl("^Before", design)) d[d$BA == "Before", ] else d
        fitted_series(dd, meta, correlation, design, n_before, n_after)
      }
    }
    if (meta$Source == "PISCO") {
      dd <- d[d$BA == "After" & d$year >= 2013L & d$year <= 2023L, ]
      common_coverage[[sid]] <- cbind(meta, common_calendar_start = 2013L, common_calendar_end = 2023L,
        n_common_years = nrow(dd), expected_years = 11L, observed_years = paste(dd$year, collapse = ","),
        eligible = nrow(dd) >= 3L, complete_calendar = nrow(dd) == 11L)
      if (nrow(dd) >= 3L) {
        annual_estimate(dd, meta, "CI_after", mean(dd$lnDiff), n_before, nrow(dd), "PISCO 2013-2023")
        for (correlation in c("AR1", "IID"))
          fitted_series(dd, meta, correlation, "CI_after", n_before, nrow(dd), "PISCO 2013-2023")
      }
    }
  }
  coverage <- dplyr::bind_rows(coverage); estimates <- dplyr::bind_rows(estimates); diagnostics <- dplyr::bind_rows(diagnostics)
  fit_attempts <- dplyr::bind_rows(fit_attempts); common_coverage <- dplyr::bind_rows(common_coverage)
  # Both design estimates must fit for a source to enter a matched comparison.
  matched <- estimates %>% dplyr::filter(analysis_window == "Full history", design %in% c("BACIPS_step", "CI_after"), fit_status == "ok") %>%
    dplyr::group_by(series_id, correlation) %>% dplyr::filter(dplyr::n_distinct(design) == 2L) %>% dplyr::ungroup()
  make_pairs <- function(x, source_subset) {
    x %>% dplyr::group_by(MPA, Taxa, Resp, Region, correlation, design) %>%
      dplyr::summarise(estimate = mean(estimate), n_sources = dplyr::n_distinct(Source),
        Sources = paste(sort(unique(Source)), collapse = ";"),
        reference_ids = paste(sort(unique(unlist(strsplit(reference_ids, ";", fixed = TRUE)))), collapse = ";"),
        reference_mapping_complete = all(reference_mapping_complete),
        reference_mapping_note = paste(sort(unique(reference_mapping_note)), collapse = "; "),
        n_boundary_sources = sum(correlation_boundary, na.rm = TRUE),
        n_apvar_problem_sources = sum(!apvar_ok, na.rm = TRUE),
        n_unstable_sources = sum(optimizer_instability, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(source_subset = source_subset)
  }
  pairs <- dplyr::bind_rows(make_pairs(matched, "All programmes"),
                     make_pairs(dplyr::filter(matched, Source == "PISCO"), "PISCO only"),
                     make_pairs(dplyr::filter(matched, Source != "Landsat"), "Divers only"),
                     make_pairs(dplyr::filter(matched, Source == "Landsat"), "Landsat only"))
  all_ci <- estimates %>% dplyr::filter(analysis_window == "Full history", design == "CI_after", fit_status == "ok")
  all_ci_pairs <- dplyr::bind_rows(make_pairs(all_ci, "All programmes"),
                            make_pairs(dplyr::filter(all_ci, Source == "PISCO"), "PISCO only"),
                            make_pairs(dplyr::filter(all_ci, Source != "Landsat"), "Divers only"),
                            make_pairs(dplyr::filter(all_ci, Source == "Landsat"), "Landsat only")) %>%
    dplyr::mutate(cohort = "All primary histories with >=3 after years")
  common_pairs <- make_pairs(dplyr::filter(estimates, analysis_window == "PISCO 2013-2023", fit_status == "ok"),
                             "PISCO only") %>% dplyr::mutate(cohort = "PISCO 2013-2023 after-only")
  design_difference <- pairs %>% dplyr::group_by(MPA, Taxa, Resp, Region, correlation, source_subset) %>%
    dplyr::summarise(design_gap = estimate[design == "CI_after"] - estimate[design == "BACIPS_step"],
      reference_ids = paste(sort(unique(unlist(strsplit(reference_ids, ";", fixed = TRUE)))), collapse = ";"),
      reference_mapping_complete = all(reference_mapping_complete), .groups = "drop")
  summarize_pairs <- function(x, value, grouping) {
    x %>% dplyr::group_by(dplyr::across(dplyr::all_of(grouping))) %>%
      dplyr::group_modify(~ review_group_summary(.x[[value]],
        review_reference_components(.x$MPA, .x$reference_ids), .x$reference_mapping_complete)) %>% dplyr::ungroup()
  }
  analysis_pairs <- dplyr::bind_rows(dplyr::mutate(pairs, cohort = "Matched before/after histories"), all_ci_pairs, common_pairs)
  grouping <- c("Taxa", "Resp", "correlation", "source_subset", "design", "cohort", "Region_scope")
  views <- dplyr::bind_rows(dplyr::mutate(analysis_pairs, Region_scope = "All"),
                    dplyr::mutate(analysis_pairs, Region_scope = Region)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping))) %>% dplyr::group_modify(~ dplyr::mutate(.x,
      reference_group_id = review_reference_components(MPA, reference_ids))) %>% dplyr::ungroup()
  reference_groups <- views %>% dplyr::group_by(dplyr::across(dplyr::all_of(c(grouping, "reference_group_id")))) %>%
    dplyr::mutate(n_mpa_in_reference_group = dplyr::n()) %>% dplyr::ungroup() %>%
    dplyr::select(dplyr::all_of(grouping), MPA, Region, Sources, n_sources, reference_group_id,
           n_mpa_in_reference_group, reference_ids, reference_mapping_complete, reference_mapping_note)
  summaries <- views %>% dplyr::group_by(dplyr::across(dplyr::all_of(grouping))) %>% dplyr::group_modify(~ {
    out <- review_group_summary(.x$estimate, .x$reference_group_id, .x$reference_mapping_complete)
    out$p_group_approx <- out$p
    if (.y$correlation != "Annual means") out$p <- NA_real_
    out$inference_status <- if (.y$correlation != "Annual means") "descriptive GLS sensitivity" else
      "exploratory annual-mean contrast; group interval is approximate"
    out$n_boundary_sources <- sum(.x$n_boundary_sources)
    out$n_apvar_problem_sources <- sum(.x$n_apvar_problem_sources)
    out$n_unstable_sources <- sum(.x$n_unstable_sources)
    out
  }) %>% dplyr::ungroup() %>% dplyr::rename(Region = Region_scope) %>%
    dplyr::group_by(correlation, source_subset, design, cohort, Region) %>%
    dplyr::mutate(p_fdr = stats::p.adjust(p, method = "BH"), fdr_family_n = sum(is.finite(p))) %>% dplyr::ungroup()
  summary <- dplyr::filter(summaries, cohort == "Matched before/after histories")
  all_ci_summary <- dplyr::filter(summaries, cohort == "All primary histories with >=3 after years")
  common_summary <- dplyr::filter(summaries, cohort == "PISCO 2013-2023 after-only")
  gap_summary <- dplyr::bind_rows(
    summarize_pairs(design_difference, "design_gap", c("Taxa", "Resp", "correlation", "source_subset")) %>% dplyr::mutate(Region = "All"),
    summarize_pairs(design_difference, "design_gap", c("Taxa", "Resp", "correlation", "source_subset", "Region")))
  gap_summary <- dplyr::select(gap_summary, -p, -iid_p)  # A gap is not a test of method validity.
  geography_summary <- function(d) {
    d$reference_group_id <- review_reference_components(d$MPA, d$reference_ids)
    a <- d[d$Region == "Islands", ]; b <- d[d$Region == "Mainland", ]
    sa <- review_group_summary(a$estimate, a$reference_group_id, a$reference_mapping_complete)
    sb <- review_group_summary(b$estimate, b$reference_group_id, b$reference_mapping_complete)
    cross_region <- length(base::intersect(unique(a$reference_group_id), unique(b$reference_group_id)))
    naive_ok <- nrow(a) >= 2L && nrow(b) >= 2L && stats::sd(a$estimate) > 0 && stats::sd(b$estimate) > 0
    naive <- if (naive_ok) stats::t.test(a$estimate, b$estimate) else NULL
    difference <- sa$estimate - sb$estimate
    se <- df <- p <- low <- high <- NA_real_
    ok <- is.finite(sa$se) && is.finite(sb$se) && cross_region == 0L && sa$se^2 + sb$se^2 > 0
    if (ok) {
      se <- sqrt(sa$se^2 + sb$se^2)
      df <- (sa$se^2 + sb$se^2)^2 /
        (sa$se^4 / sa$reference_group_df + sb$se^4 / sb$reference_group_df)
      low <- difference - stats::qt(0.975, df) * se; high <- difference + stats::qt(0.975, df) * se
      p <- 2 * stats::pt(-abs(difference / se), df)
    }
    data.frame(n_islands = nrow(a), n_mainland = nrow(b),
      n_islands_reference_groups = sa$n_reference_groups, n_mainland_reference_groups = sb$n_reference_groups,
      cross_region_reference_groups = cross_region, islands_minus_mainland = difference,
      se = se, reference_group_df = df, ci_low = low, ci_high = high, p = p,
      iid_ci_low = if (naive_ok) naive$conf.int[1] else NA_real_,
      iid_ci_high = if (naive_ok) naive$conf.int[2] else NA_real_, iid_p = if (naive_ok) naive$p.value else NA_real_,
      few_reference_groups = sa$few_reference_groups || sb$few_reference_groups,
      reference_mapping_complete = all(d$reference_mapping_complete),
      uncertainty_method = "regional delete-reference-group jackknife variances; Welch t degrees of freedom",
      uncertainty_status = if (cross_region > 0L) "cross-region reference dependence: interval unavailable" else if (!ok)
        "fewer than 3 groups in a region or zero variance: interval unavailable" else if (sa$few_reference_groups || sb$few_reference_groups)
        "few groups: approximate descriptive regional uncertainty" else "approximate group-aware regional uncertainty")
  }
  geographic <- analysis_pairs %>% dplyr::group_by(Taxa, Resp, correlation, source_subset, design, cohort) %>%
    dplyr::group_modify(~ geography_summary(.x)) %>% dplyr::ungroup() %>%
    dplyr::mutate(p_group_approx = p, p = ifelse(correlation == "Annual means", p, NA_real_),
           inference_status = ifelse(correlation == "Annual means", "exploratory geographic comparison", "descriptive GLS sensitivity")) %>%
    dplyr::group_by(correlation, source_subset, design, cohort) %>%
    dplyr::mutate(p_fdr = stats::p.adjust(p, method = "BH"), fdr_family_n = sum(is.finite(p))) %>% dplyr::ungroup()
  influence <- function(d, kind) {
    full <- review_group_summary(d$estimate, d$reference_group_id, d$reference_mapping_complete)
    units <- if (kind == "MPA") d$MPA else unique(d$reference_group_id)
    dplyr::bind_rows(lapply(units, function(id) {
      omitted <- if (kind == "MPA") d$MPA == id else d$reference_group_id == id
      remaining <- d[!omitted, ]
      out <- review_group_summary(remaining$estimate, remaining$reference_group_id, remaining$reference_mapping_complete)
      out$p <- NULL; out$iid_p <- NULL
      cbind(omitted_unit = id, omitted_mpas = paste(sort(d$MPA[omitted]), collapse = ";"),
        n_removed_mpa = sum(omitted), full_n_mpa = full$n_mpa, full_n_reference_groups = full$n_reference_groups,
        full_estimate = full$estimate, estimate_change = out$estimate - full$estimate,
        sign_changed = if (is.finite(out$estimate)) sign(out$estimate) != sign(full$estimate) else NA,
        inference_status = "descriptive deletion influence; no additional hypothesis tests", out)
    }))
  }
  annual_views <- dplyr::filter(views, correlation == "Annual means")
  leave_mpa <- annual_views %>% dplyr::group_by(dplyr::across(dplyr::all_of(grouping))) %>% dplyr::group_modify(~ influence(.x, "MPA")) %>% dplyr::ungroup()
  leave_reference_group <- annual_views %>% dplyr::group_by(dplyr::across(dplyr::all_of(grouping))) %>%
    dplyr::group_modify(~ influence(.x, "reference group")) %>% dplyr::ungroup()
  primary_influence_diagnostics <- annual_views %>% dplyr::group_by(dplyr::across(dplyr::all_of(grouping))) %>% dplyr::group_modify(~ {
    x <- .x$estimate
    data.frame(n_mpa = length(x), n_reference_groups = dplyr::n_distinct(.x$reference_group_id),
      shapiro_p = if (length(x) >= 3L && stats::sd(x) > 0) stats::shapiro.test(x)$p.value else NA_real_,
      max_abs_standardized_mpa_contrast = if (length(x) >= 2L && stats::sd(x) > 0) max(abs((x - mean(x)) / stats::sd(x))) else NA_real_,
      diagnostic_status = "distribution screening; small n limits power; p is not a hypothesis family")
  }) %>% dplyr::ungroup()
  primary_summary <- primary %>% dplyr::mutate(Region = ifelse(Location == "C", "Islands", "Mainland")) %>%
    dplyr::group_by(MPA, Taxa, Resp, Region, AnalysisType) %>% dplyr::summarise(estimate = mean(Mean), .groups = "drop") %>%
    dplyr::group_by(Taxa, Resp, Region, AnalysisType) %>% dplyr::group_modify(~ review_summary(.x$estimate)) %>% dplyr::ungroup() %>%
    dplyr::rename(iid_se = se, iid_ci_low = ci_low, iid_ci_high = ci_high, iid_p = p) %>%
    dplyr::mutate(inference_status = "descriptive original-model subgroup; different estimands/cohorts; iid uncertainty audit only")
  source_support <- input %>% dplyr::distinct(CA_MPA_Name_Short, Taxa, resp, source, Region) %>%
    dplyr::count(Taxa, resp, source, Region, name = "n_mpa")
  outputs <- list(table_s_design_matched_coverage = coverage, table_s_design_matched_estimates = estimates,
                  table_s_design_matched_diagnostics = diagnostics, table_s_design_matched_pairs = pairs,
                  table_s_design_matched_summary = summary, table_s_design_matched_gap = gap_summary,
                  table_s_design_geography = geographic, table_s_design_primary_subsets = primary_summary,
                  table_s_design_source_support = source_support, table_s_design_all_ci_pairs = all_ci_pairs,
                  table_s_design_all_ci_summary = all_ci_summary,
                  table_s_design_group_summary = summaries,
                  table_s_design_reference_sites = reference_sites,
                  table_s_design_reference_groups = reference_groups,
                  table_s_design_reference_provenance = mapping_provenance,
                  table_s_design_leave_mpa_out = leave_mpa,
                  table_s_design_leave_reference_group_out = leave_reference_group,
                  table_s_design_annual_mean_diagnostics = primary_influence_diagnostics,
                  table_s_design_gls_fit_attempts = fit_attempts,
                  table_s_design_pisco_common_calendar_coverage = common_coverage,
                  table_s_design_pisco_common_calendar_pairs = common_pairs,
                  table_s_design_pisco_common_calendar_summary = common_summary)
  for (name in names(outputs)) utils::write.csv(outputs[[name]], here::here("tables", paste0(name, ".csv")), row.names = FALSE)
  saveRDS(list(input = input, matched = matched, pairs = pairs, diagnostics = diagnostics,
               reference_sites = reference_sites, reference_groups = reference_groups,
               fit_attempts = fit_attempts, group_summary = summaries, common_pairs = common_pairs),
          here::here("data", "cache", "designation_design_checks.rds"))

  # Each point is one MPA, with programmes averaged within the MPA.
  pd <- pairs %>% dplyr::filter(correlation == "Annual means", source_subset == "All programmes") %>%
    dplyr::mutate(label = paste(taxon_display(Taxa), ifelse(Resp == "Bio", "biomass", "density"), sep = "\n"),
           design = factor(design, c("CI_after", "BACIPS_step"), c("After-only CI", "Before-after BACIPS")))
  p <- ggplot2::ggplot(pd, ggplot2::aes(design, estimate, group = MPA)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3) +
    ggplot2::geom_line(color = "grey75", linewidth = 0.35) +
    ggplot2::geom_point(ggplot2::aes(color = Region, shape = Region), size = 1.4) +
    ggplot2::facet_wrap(~ label, ncol = 3, scales = "free_y") +
    ggplot2::scale_color_manual(values = c(Islands = "#0072B2", Mainland = "#D55E00")) +
    ggplot2::labs(x = NULL, y = "Log response ratio / change in log response ratio") + theme_mpa(base_size = 8) +
    ggplot2::theme(legend.position = "bottom", axis.text.x = ggplot2::element_text(size = 6.5), plot.margin = ggplot2::margin(6, 6, 6, 6))
  review_export_figure(p, "fig_s_design_matched", height = 185)
  writeLines(c("# Matched designation-design comparison", "",
    "Each connected pair of points is one MPA. After-only CI estimates the post-designation level of ln(MPA/reference); stationary BACIPS estimates its before-to-after change. These are different estimands: their discrepancy can reflect a pre-existing reef difference and is not a ranking of statistical methods.",
    "Primary estimates are arithmetic means of annual source-specific lnRR: after-only CI is the after-period mean; stationary BACIPS subtracts the before-period mean. Source-specific estimates are averaged within MPA before plotting or network summaries. Blue circles mark islands; orange triangles mark mainland MPAs. Zero denotes parity (CI) or no change (BACIPS). Facets distinguish taxa and biomass/density responses; y-axis ranges differ.",
    "The figure shows points and connecting lines, without confidence intervals. In table_s_design_matched_summary.csv and the cohort-indexed table_s_design_group_summary.csv, se/ci_low/ci_high are approximate reference-group-aware intervals. Connected groups are derived from declared shared reference-site IDs among the sources entering each comparison. The equal-MPA point estimate is unchanged. The uncertainty is a delete-connected-reference-group jackknife: SE squared = (G-1)/G times the sum of squared deviations of leave-group-out estimates from their mean; intervals use t critical values with G-1 degrees of freedom. The previous MPA t results are retained only as iid_se/iid_ci_low/iid_ci_high/iid_p audit columns. Annual years and overlapping programmes do not increase MPA replication.",
    "Fewer than three reference groups yields no group interval or p-value. Fewer than ten groups is flagged as few_reference_groups; the reported intervals and p-values are approximate descriptive sensitivities, not well-calibrated evidence with few groups. Unequal group sizes further limit the jackknife approximation. The mapped groups account for documented sharing, not all spatial dependence. PISCO reference assignments are checked against survey-site/year footprints without using response values. Biomass site availability and individual MBON reference contributions are not reconstructed. Landsat records lack globally identifiable reference polygons: MPA labels remain distinct and polygon sharing/overlap is explicitly unverified. reference_mapping_complete flags this limitation. Site assignments, inclusion/exclusion decisions, comparison-specific groups, source files and MD5s are exported.",
    "AR(1) GLS indexed by calendar year and IID GLS remain descriptive model sensitivities. Each AR(1) fit tries starting correlations -0.75, -0.25, 0, 0.5 and 0.9 with nlminb and optim/BFGS, selecting the largest finite REML log likelihood with valid coefficient covariance. IID fits try both optimizers without changing the covariance model. Every attempt, error and warning is exported; no AR(1)-to-IID fallback occurs. Boundary correlation (absolute phi > 0.98), approximate variance-covariance failure (apvar_ok), and optimizer instability are separate flags. Optimizer instability denotes fewer than two successful optimizers, best-by-optimizer log-likelihood disagreement >0.01, or coefficient range >0.05 among attempts within 0.0001 of the best log likelihood. Conditional GLS SE/CI/p values remain in the source-estimate audit. GLS network summaries have p/p_fdr unavailable and p_group_approx retained only as a descriptive audit. Failed fits require both designs to succeed in their matched cohort; arithmetic primary estimates never drop these series because of GLS failure.",
    "The pre-existing pBACIPS primary models retain their progressive, year-11 or step estimands. Their separate subgroup table is descriptive and uses different cohorts; it is not the matched-design experiment.",
    "Regional comparisons are descriptive because programmes and designation dates covary with geography. Matched PISCO before/after histories provide no island comparison. A separate PISCO 2013-2023 after-only cohort uses the same calendar window in both regions, retaining at least three observed paired years per source series and auditing missing years. It does not require all eleven years or equate reserve age. Regional intervals combine the two regional group-jackknife variances with Welch degrees of freedom; no interval is emitted with fewer than three groups in either region or a documented cross-region shared group.",
    "Leave-one-MPA and leave-one-reference-group outputs retain the full cohort's connected groups for each deletion and recompute the equal-MPA mean among retained MPAs. Deleting a group removes all its MPA programme averages, rather than deleting only its PISCO contribution. These are influence checks, without additional hypothesis-test p-values. Annual-mean distribution diagnostics supplement the GLS normalized-residual normality, calendar-lag correlation and period-variance screens. Diagnostic p-values are screens, not inferential endpoints.",
    "", "## Output schema", "",
    "Summary keys: Taxa, Resp, correlation, source_subset, design, cohort, Region. Filter correlation = Annual means for the primary bridge and cohort = Matched before/after histories for the matched comparison. uncertainty_status, few_reference_groups, reference_mapping_complete and reference_group_size_ratio must accompany interval interpretation. BH families contain the finite annual-mean endpoint p-values within correlation/source_subset/design/cohort/Region; source and window subsets are exploratory and no cross-subset multiplicity claim is made.",
    "", "- table_s_design_group_summary.csv: preferred cohort-indexed group-aware summaries; matched_summary/all_ci_summary/pisco_common_calendar_summary are cohort-specific views.",
    "- table_s_design_reference_sites.csv: series_id, MPA, Source, Taxa, Resp, site_id, reference_id, reference_years, included, excluded, exclusion_reason, identity_verified, contribution_audit, mapping_file, mapping_md5. PISCO reference_years are surveyed years overlapping the harmonized series; MBON/Landsat years are series coverage only.",
    "- table_s_design_reference_groups.csv: cohort, Region_scope, MPA, Sources, reference_group_id, n_mpa_in_reference_group, reference_ids and mapping limitations, with the other summary keys.",
    "Routine script 27 reads only tracked data/harmonized/reference_sites.csv, reference_primary_keys.csv and mapping_provenance.csv for this audit. It validates current input MD5s, every primary key/model/year classification and the vendored footprint checksums before fitting. The separate export_reference_registry.R is an explicit one-time build from source mappings, survey footprints and exclusion definitions; it is not a routine-pipeline dependency. Source_audit paths are provenance only and are never opened by script 27. See data/harmonized/REFERENCE_REGISTRY.md for rebuild and audit limitations.",
    "- table_s_design_reference_provenance.csv: role, kind, path and md5 for the harmonized inputs, primary-results snapshot, vendored registry/key files, source-audit mappings/footprints/exclusion definitions and exporter.",
    "- table_s_design_leave_mpa_out.csv / table_s_design_leave_reference_group_out.csv: annual-mean summary keys with Region_scope, omitted_unit, omitted_mpas, n_removed_mpa, full_n_mpa, full_n_reference_groups, full_estimate, estimate, estimate_change, sign_changed, retained group counts and approximate intervals.",
    "- table_s_design_gls_fit_attempts.csv: series_id, correlation, design, analysis_window, attempt_id, optimizer, start_phi, success, log_lik, selected, ar1_phi, correlation_boundary, apvar_ok, apvar_status, error and warnings. matched_estimates adds selection diagnostics and source-level reference identities.",
    "- table_s_design_pisco_common_calendar_coverage.csv / pairs.csv / summary.csv: 2013-2023 post-only eligibility, equal-MPA estimates and group-aware uncertainty. The geography table includes this cohort.",
    "- table_s_design_annual_mean_diagnostics.csv: MPA-level distribution screens for every arithmetic cohort/source/region/design view."),
    here::here("docs", "design_geography_legend.md"))
  cat("  Eligible histories:", sum(coverage$matched_eligible), "of", nrow(coverage),
      "; AR1 fit failures:", sum(estimates$correlation == "AR1" & estimates$fit_status == "failed"), "\n")
})
