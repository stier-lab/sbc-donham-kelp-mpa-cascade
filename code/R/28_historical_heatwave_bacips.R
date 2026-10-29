# Compare kelp disturbance responses at the same Landsat pairs before and after
# MPA designation. Event differences are observational, not randomized MPA effects.

local({
  suppressPackageStartupMessages({ library(here); library(dplyr); library(tidyr); library(ggplot2) })
  `%>%` <- dplyr::`%>%`
  source(here::here("code", "R", "00c_analysis_constants.R"), local = TRUE)
  source(here::here("code", "R", "00b_color_palette.R"), local = TRUE)
  source(here::here("code", "R", "review_concern_utils.R"), local = TRUE)
  cat("Running historical-event before/after protection check (script 28)...\n")
  annual <- utils::read.csv(here::here("data", "harmonized", "harmonized_landsat_annual_index.csv"))
  metadata <- utils::read.csv(here::here("data", "harmonized", "harmonized_site_metadata.csv"))
  primary <- utils::read.csv(here::here("data", "sumstats_final.csv"))
  mpas <- unique(primary$MPA[primary$Source == "Landsat"])
  annual <- annual %>% dplyr::filter(CA_MPA_Name_Short %in% mpas) %>%
    dplyr::left_join(metadata %>% dplyr::select(CA_MPA_Name_Short, MPA_Start, Location), by = "CA_MPA_Name_Short") %>%
    dplyr::mutate(Region = ifelse(Location == "C", "Islands", "Mainland"))
  stopifnot(!anyNA(annual$MPA_Start), !anyDuplicated(annual[c("CA_MPA_Name_Short", "year", "status")]),
            all(annual$index >= 0), all(annual$index <= 1 + 1e-12))
  # The historical recovery window must end before designation. The modern
  # disturbance must begin after designation; baseline overlap is documented.
  specifications <- list(
    primary = list(historical_base = 1993:1996, historical_event = 1997:1998, historical_recovery = 1999:2001,
                   modern_base = 2010:2013, modern_event = 2014:2016, modern_recovery = 2017:2019, minimum = 2L),
    equal_event_duration = list(historical_base = 1993:1996, historical_event = 1997:1998, historical_recovery = 1999:2001,
                   modern_base = 2010:2013, modern_event = 2014:2015, modern_recovery = 2017:2019, minimum = 2L),
    strict_complete_windows = list(historical_base = 1993:1996, historical_event = 1997:1998, historical_recovery = 1999:2001,
                   modern_base = 2010:2013, modern_event = 2014:2016, modern_recovery = 2017:2019, minimum = Inf),
    post_designation_baseline = list(historical_base = 1993:1996, historical_event = 1997:1998, historical_recovery = 1999:2001,
                   modern_base = 2013L, modern_event = 2014:2016, modern_recovery = 2017:2019, minimum = 2L),
    earlier_historical_baseline = list(historical_base = 1991:1996, historical_event = 1997:1998, historical_recovery = 1999:2001,
                   modern_base = 2010:2013, modern_event = 2014:2016, modern_recovery = 2017:2019, minimum = 2L))
  coverage <- pairs <- annual_windows <- list()
  for (spec_name in names(specifications)) {
    spec <- specifications[[spec_name]]
    for (mpa in sort(mpas)) {
      d <- annual[annual$CA_MPA_Name_Short == mpa, ]
      start <- unique(d$MPA_Start)
      if (length(start) != 1L) stop("Missing or ambiguous MPA start: ", mpa)
      timing_ok <- max(spec$historical_recovery) <= start && min(spec$modern_event) > start
      for (metric in c("resistance", "recovery")) {
        window_names <- c("historical_base", if (metric == "resistance") "historical_event" else "historical_recovery",
                          "modern_base", if (metric == "resistance") "modern_event" else "modern_recovery")
        w <- dplyr::bind_rows(lapply(window_names, function(wn) {
          d %>% dplyr::filter(year %in% spec[[wn]]) %>% dplyr::group_by(status) %>%
            dplyr::summarise(mean_index = mean(index), n_years = dplyr::n_distinct(year),
                      years = paste(sort(year), collapse = ","), .groups = "drop") %>%
            dplyr::mutate(window = wn, expected_years = length(spec[[wn]]),
                   required_years = min(spec$minimum, length(spec[[wn]])))
        }))
        counts_ok <- nrow(w) == 8L && all(w$n_years >= w$required_years)
        positive_ok <- counts_ok && all(is.finite(w$mean_index) & w$mean_index > 0)
        eligible <- timing_ok && positive_ok
        reason <- if (!timing_ok) "event/recovery does not bracket designation" else if (!counts_ok)
          "incomplete window coverage" else if (!positive_ok) "zero or nonfinite window mean; log contrast undefined" else "eligible"
        coverage[[paste(spec_name, mpa, metric)]] <- data.frame(specification = spec_name, MPA = mpa,
          Region = d$Region[1], MPA_Start = start, metric = metric, eligible = eligible, reason = reason,
          modern_baseline_contains_designation = start %in% spec$modern_base,
          minimum_years_required = if (is.infinite(spec$minimum)) "all" else as.character(spec$minimum),
          minimum_years_observed = if (nrow(w)) min(w$n_years) else 0L)
        annual_windows[[paste(spec_name, mpa, metric)]] <- cbind(specification = spec_name, MPA = mpa, metric = metric, w)
        if (!eligible) next
        val <- function(status, window) w$mean_index[w$status == status & w$window == window]
        htarget <- window_names[2]; mtarget <- window_names[4]
        hmpa <- val("mpa", htarget) / val("mpa", "historical_base")
        href <- val("reference", htarget) / val("reference", "historical_base")
        mmpa <- val("mpa", mtarget) / val("mpa", "modern_base")
        mref <- val("reference", mtarget) / val("reference", "modern_base")
        h <- log(hmpa / href); m <- log(mmpa / mref)
        pairs[[paste(spec_name, mpa, metric)]] <- data.frame(specification = spec_name, MPA = mpa,
          Region = d$Region[1], MPA_Start = start, metric = metric,
          historical_mpa = hmpa, historical_reference = href,
          modern_mpa = mmpa, modern_reference = mref,
          historical_log_advantage = h, modern_log_advantage = m,
          change_log_advantage = m - h)
      }
    }
  }
  coverage <- dplyr::bind_rows(coverage); pairs <- dplyr::bind_rows(pairs); annual_windows <- dplyr::bind_rows(annual_windows)
  if (!nrow(pairs)) stop("No eligible historical/modern Landsat pairs.")
  summarize <- function(d) {
    out <- review_summary(d$change_log_advantage)
    out$historical_mean_log_advantage <- mean(d$historical_log_advantage)
    out$modern_mean_log_advantage <- mean(d$modern_log_advantage)
    out$ratio_of_advantages <- exp(out$estimate)
    out$ratio_ci_low <- exp(out$ci_low); out$ratio_ci_high <- exp(out$ci_high)
    out$p_signflip <- review_signflip(d$change_log_advantage)
    out$p_wilcoxon <- if (nrow(d) >= 2L) suppressWarnings(stats::wilcox.test(d$change_log_advantage, exact = FALSE)$p.value) else NA_real_
    out$n_positive_change <- sum(d$change_log_advantage > 0)
    out$shapiro_p <- if (nrow(d) >= 3L && stats::sd(d$change_log_advantage) > 0)
      stats::shapiro.test(d$change_log_advantage)$p.value else NA_real_
    out
  }
  summary <- dplyr::bind_rows(
    pairs %>% dplyr::group_by(specification, metric) %>% dplyr::group_modify(~ summarize(.x)) %>% dplyr::mutate(Region = "All"),
    pairs %>% dplyr::group_by(specification, metric, Region) %>% dplyr::group_modify(~ summarize(.x))) %>% dplyr::ungroup() %>%
    dplyr::mutate(test_family = dplyr::case_when(
      Region == "All" & specification == "primary" ~ "Primary historical comparison (2 endpoints)",
      Region == "All" ~ "Exploratory window sensitivities (8 endpoints)",
      TRUE ~ "Exploratory geographic subsets (20 endpoints)")) %>%
    dplyr::group_by(test_family) %>% dplyr::mutate(p_fdr = stats::p.adjust(p, method = "BH"),
                                    fdr_family_n = sum(is.finite(p))) %>% dplyr::ungroup()
  loo <- pairs %>% dplyr::group_by(specification, metric) %>% dplyr::group_modify(~ {
    if (nrow(.x) < 3L) return(data.frame())
    dplyr::bind_rows(lapply(seq_len(nrow(.x)), function(i) cbind(omitted_mpa = .x$MPA[i], summarize(.x[-i, ]))))
  }) %>% dplyr::ungroup()
  # Protection/reference changes within each event, using the same eligible pairs.
  event_comparison <- pairs %>% dplyr::group_by(specification, metric, Region) %>% dplyr::group_modify(~ {
    dplyr::bind_rows(lapply(c("historical", "modern"), function(event) {
      x <- .x[[paste0(event, "_log_advantage")]]
      cbind(event = event, review_summary(x), p_signflip = review_signflip(x))
    }))
  }) %>% dplyr::ungroup()
  # Hold eligible MPAs fixed to distinguish baseline choice from changing cohort.
  common_cohort <- dplyr::bind_rows(lapply(base::setdiff(names(specifications), "primary"), function(spec_name) {
    dplyr::bind_rows(lapply(c("resistance", "recovery"), function(mt) {
      p <- pairs %>% dplyr::filter(specification == "primary", metric == mt)
      s <- pairs %>% dplyr::filter(specification == spec_name, metric == mt)
      common <- base::intersect(p$MPA, s$MPA)
      dplyr::bind_rows(lapply(c("primary", spec_name), function(spec) {
        d <- pairs %>% dplyr::filter(specification == spec, metric == mt, MPA %in% common)
        cbind(sensitivity = spec_name, specification = spec, metric = mt,
              cohort = "Intersection of primary and sensitivity", summarize(d))
      }))
    }))
  }))
  outputs <- list(table_s_historical_heatwave_coverage = coverage, table_s_historical_heatwave_windows = annual_windows,
                  table_s_historical_heatwave_pairs = pairs, table_s_historical_heatwave_summary = summary,
                  table_s_historical_heatwave_leave_one_out = loo, table_s_historical_heatwave_event_comparison = event_comparison,
                  table_s_historical_heatwave_common_cohort = common_cohort)
  for (name in names(outputs)) utils::write.csv(outputs[[name]], here::here("tables", paste0(name, ".csv")), row.names = FALSE)
  saveRDS(list(annual = annual, specifications = specifications, pairs = pairs),
          here::here("data", "cache", "historical_heatwave_checks.rds"))

  pp <- pairs %>% dplyr::filter(specification == "primary") %>%
    dplyr::select(MPA, Region, metric, historical_log_advantage, modern_log_advantage) %>%
    tidyr::pivot_longer(dplyr::ends_with("log_advantage"), names_to = "event", values_to = "value") %>%
    dplyr::mutate(event = factor(event, c("historical_log_advantage", "modern_log_advantage"),
                          c("1997-1998\nBefore designation", "2014-2016\nAfter designation")),
           metric = factor(metric, c("resistance", "recovery"), c("Resistance", "Recovery")))
  p <- ggplot2::ggplot(pp, ggplot2::aes(event, value, group = MPA)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3) +
    ggplot2::geom_line(color = "grey70", linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(color = Region, shape = Region), size = 1.7) +
    ggplot2::facet_wrap(~ metric, nrow = 1) +
    ggplot2::scale_color_manual(values = c(Islands = "#0072B2", Mainland = "#D55E00")) +
    ggplot2::labs(x = NULL, y = "Log(MPA / reference) baseline-relative response") +
    theme_mpa(base_size = 9) + ggplot2::theme(legend.position = "bottom", plot.margin = ggplot2::margin(8, 8, 6, 6))
  review_export_figure(p, "fig_s_historical_heatwave_bacips", height = 95)
  writeLines(c("# Historical versus modern disturbance response", "",
    "Each line joins the same MPA-reference polygon pair's response to the 1997-1998 El Nino before designation and the 2014-2016 marine heatwave after designation. Blue circles denote islands and orange triangles denote mainland MPAs. Each point is log[(MPA event-window mean / MPA baseline mean)/(reference event-window mean / reference baseline mean)]. Zero means equal baseline-relative response. Panels show resistance and early recovery.",
    "Historical baseline: 1993-1996; disturbance: 1997-1998; recovery: 1999-2001. Modern baseline: 2010-2013; disturbance: 2014-2016; recovery: 2017-2019. At least two observed paired years are required in every window for the primary analysis. The complete historical recovery must precede designation, and the modern disturbance must follow it. Only the 16 Landsat MPAs already in the primary paper cohort are candidates.",
    "The MPA unit enters once in each metric, and the within-pair event difference is the unit of inference. Report equal-weight pair mean differences, t confidence intervals, exact sign-flip sensitivity, Wilcoxon sensitivity and leave-one-pair-out checks. The two primary endpoint p-values form one BH family. The eight alternative-window endpoints form a separate exploratory family; all 20 geographic-subset endpoints form a third family. These analyses were added during coauthor revision. Sign-flip inference assumes symmetry of pair-level differences; it is not randomization-based causal inference.",
    "True zero annual observations remain zero; no offset is added. A zero window mean is audited and excluded only because its log contrast is undefined. This is an observational comparison of different historical events. Event intensity, duration, regional oceanography and contemporaneous disease differ; fishing effort is unavailable. The 2010-2013 baseline also overlaps 2012 designation for mainland MPAs. A 2013-only post-designation baseline is shown as a sensitivity with reduced baseline precision; other windows still require two years. The common-cohort table holds MPA membership fixed for each window sensitivity. Cabrillo drops from the 2013-baseline comparison because its baseline is zero; it is retained in the primary analysis. Leave-one-out results and Shapiro tests disclose skew and influence, not grounds for deleting observations. Persistence of an advantage does not establish a causal protection effect."),
    here::here("docs", "historical_heatwave_legend.md"))
  cat("  Primary comparisons:\n")
  print(summary %>% dplyr::filter(specification == "primary", Region == "All") %>%
          dplyr::select(metric, n_mpa, estimate, ci_low, ci_high, p, p_signflip, p_fdr))
})
