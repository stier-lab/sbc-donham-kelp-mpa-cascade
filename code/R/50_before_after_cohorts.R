# Coverage-defined cohorts, not a comparison of selected CI/pBACIPS models.
run_before_after_cohorts <- function(root = here::here()) {
  suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})
  source(file.path(root, "code/R/review_concern_utils.R"), local = TRUE)
  source(file.path(root, "code/R/00c_analysis_constants.R"), local = TRUE)
  paths <- c("data/harmonized/harmonized_response_ratios.csv", "data/harmonized/harmonized_landsat_rr.csv",
    "data/sumstats_final.csv", "data/harmonized/reference_sites.csv", "data/harmonized/reference_primary_keys.csv",
    "data/harmonized/mapping_provenance.csv", "code/R/50_before_after_cohorts.R", "code/R/review_concern_utils.R")
  hashes <- tools::md5sum(file.path(root, paths))
  read <- function(p) read.csv(file.path(root, p))
  out <- function(x, suffix) write.csv(x, file.path(root, "tables", paste0("table_s_coverage_cohort_", suffix, ".csv")), row.names = FALSE)
  input <- review_design_input(bind_rows(read(paths[1]), read(paths[2])), read(paths[3]), RESILIENCE_TAXA_SHORT)
  registry <- review_read_reference_registry(input, root)
  ref <- registry$reference_sites %>% filter(included) %>% group_by(series_id) %>%
    summarise(reference_ids = paste(sort(unique(reference_id)), collapse = ";"),
      reference_mapping_complete = all(identity_verified), .groups = "drop")
  series <- input %>% group_by(series_id, MPA = CA_MPA_Name_Short, Taxa, Resp = resp, Source = source, Region) %>%
    summarise(n_before = sum(BA == "Before"), n_after = sum(BA == "After"),
      before_mean = if (any(BA == "Before")) mean(lnDiff[BA == "Before"]) else NA_real_,
      after_mean = mean(lnDiff[BA == "After"]),
      before_years = paste(year[BA == "Before"], collapse = ";"),
      after_years = paste(year[BA == "After"], collapse = ";"), .groups = "drop") %>%
    left_join(ref, by = "series_id") %>% mutate(change = after_mean - before_mean,
      eligible_after = n_after >= 3L, before_coverage = case_when(n_before == 0 ~ "No before years",
        n_before < 3 ~ "One or two before years", TRUE ~ "At least three before years"))
  stopifnot(nrow(series) == n_distinct(input$series_id), !anyNA(series$reference_ids))
  out(series, "series")
  subsets <- list("All programmes" = series, "Divers only" = filter(series, Source != "Landsat"),
    "Landsat only" = filter(series, Source == "Landsat"))
  for (s in unique(series$Source)) subsets[[paste(s, "only")]] <- filter(series, Source == s)
  rows <- rosters <- list()
  for (view in names(subsets)) for (minimum in c(3L, 1L)) {
    d <- subsets[[view]] %>% filter(eligible_after) %>% mutate(eligible_before = n_before >= minimum)
    roster <- d %>% group_by(MPA, Taxa, Resp, Region) %>% summarise(
      n_sources = n(), n_eligible_sources = sum(eligible_before),
      max_before_years = max(n_before), min_before_years = min(n_before),
      Sources = paste(sort(Source), collapse = ";"),
      reference_ids = paste(sort(unique(unlist(strsplit(reference_ids, ";", fixed = TRUE)))), collapse = ";"),
      reference_mapping_complete = all(reference_mapping_complete), .groups = "drop") %>%
      mutate(cohort = ifelse(n_eligible_sources > 0, "With before-after coverage", "Without sufficient before coverage"),
        mixed_coverage = n_eligible_sources > 0 & n_eligible_sources < n_sources,
        source_subset = view, min_before = minimum)
    rosters[[length(rosters) + 1L]] <- roster
    d <- left_join(d, select(roster, MPA, Taxa, Resp, cohort), by = c("MPA", "Taxa", "Resp"))
    for (metric in c("Post-designation difference", "Post difference, eligible sources", "Before-after change")) {
      z <- if (metric == "Post-designation difference") d else filter(d, eligible_before)
      z$value <- if (metric == "Before-after change") z$change else z$after_mean
      rows[[length(rows) + 1L]] <- z %>% group_by(MPA, Taxa, Resp, Region, cohort) %>%
        summarise(estimate = mean(value), n_sources = n(), Sources = paste(sort(Source), collapse = ";"),
          reference_ids = paste(sort(unique(unlist(strsplit(reference_ids, ";", fixed = TRUE)))), collapse = ";"),
          reference_mapping_complete = all(reference_mapping_complete), .groups = "drop") %>%
        mutate(metric = metric, source_subset = view, min_before = minimum)
    }
  }
  roster <- bind_rows(rosters); pairs <- bind_rows(rows)
  stopifnot(!anyDuplicated(roster[c("MPA", "Taxa", "Resp", "source_subset", "min_before")]))
  views <- bind_rows(mutate(pairs, region_scope = "All"), mutate(pairs, region_scope = Region))
  summary <- views %>% group_by(Taxa, Resp, source_subset, min_before, region_scope, cohort, metric) %>%
    group_modify(~ review_group_summary(.x$estimate, review_reference_components(.x$MPA, .x$reference_ids),
      .x$reference_mapping_complete)) %>% ungroup() %>% select(-p, -iid_p) %>%
    mutate(inference = "Exploratory descriptive cohort comparison; no significance tests")
  # Delete shared-reference components jointly across both coverage cohorts.
  difference <- views %>% filter(metric == "Post-designation difference") %>%
    group_by(Taxa, Resp, source_subset, min_before, region_scope) %>% group_modify(~ {
      d <- .x; g <- review_reference_components(d$MPA, d$reference_ids)
      yes <- d$cohort == "With before-after coverage"
      contrast <- function(keep) mean(d$estimate[keep & yes]) - mean(d$estimate[keep & !yes])
      groups <- unique(g); ng1 <- n_distinct(g[yes]); ng0 <- n_distinct(g[!yes])
      est <- if (any(yes) && any(!yes)) contrast(rep(TRUE, nrow(d))) else NA_real_
      se <- low <- high <- NA_real_
      if (ng1 >= 3 && ng0 >= 3) {
        jack <- vapply(groups, function(id) contrast(g != id), numeric(1))
        if (all(is.finite(jack))) {
          se <- sqrt((length(groups) - 1) / length(groups) * sum((jack - mean(jack))^2))
          low <- est - qt(.975, length(groups) - 1) * se
          high <- est + qt(.975, length(groups) - 1) * se
        }
      }
      data.frame(n_with = sum(yes), n_without = sum(!yes), groups_with = ng1, groups_without = ng0,
        cross_cohort_groups = length(intersect(g[yes], g[!yes])), estimate = est, se, ci_low = low, ci_high = high,
        reference_mapping_complete = all(d$reference_mapping_complete),
        uncertainty_status = if (is.finite(se)) "Approximate joint reference-group jackknife; descriptive" else
          "Interval unavailable: fewer than three reference groups in a cohort")
    }) %>% ungroup()
  out(roster, "roster"); out(pairs, "pairs"); out(summary, "summary"); out(difference, "difference")
  coverage <- roster %>% group_by(source_subset, min_before, Taxa, Resp, cohort, Region) %>%
    summarise(n_mpa = n(), n_mixed = sum(mixed_coverage), n_zero_before = sum(max_before_years == 0),
      n_short_before = sum(max_before_years > 0 & max_before_years < min_before), .groups = "drop")
  out(coverage, "coverage")
  labels <- c("M. pyrifera" = "Giant kelp", "M. franciscanus" = "Red urchin", "S. purpuratus" = "Purple urchin",
    "S. pulcher" = "Sheephead", "P. interruptus" = "Lobster")
  plot_rows <- pairs %>% filter(min_before == 3, source_subset == "All programmes", metric == "Post-designation difference") %>%
    mutate(display = ifelse(cohort == "With before-after coverage", "With before\ndata", "Limited before\ndata"),
      display = factor(display, c("With before\ndata", "Limited before\ndata")),
      label = paste(unname(labels[Taxa]), ifelse(Resp == "Bio", "biomass", "density")))
  stopifnot(!anyNA(labels[pairs$Taxa]))
  means <- plot_rows %>% group_by(label, display) %>% summarise(estimate = mean(estimate), n = n(), .groups = "drop") %>%
    complete(label, display, fill = list(n = 0L))
  p <- ggplot(plot_rows, aes(display, estimate)) + geom_hline(yintercept = 0, color = "grey65", linetype = 2) +
    geom_point(aes(color = Region, shape = Region), position = position_jitter(width = .13, height = 0, seed = 50), size = 2, alpha = .8) +
    geom_point(data = means, shape = 95, size = 9, color = "black", na.rm = TRUE) +
    geom_text(data = means, aes(y = Inf, label = paste0("n=", n)), vjust = 1.4, size = 3) +
    facet_wrap(~label, ncol = 3, scales = "free_y") +
    scale_y_continuous(expand = expansion(mult = c(.06, .22))) +
    scale_color_manual(values = c("Islands" = "#007F73", "Mainland" = "#BA4A58")) +
    scale_shape_manual(values = c("Islands" = 16, "Mainland" = 17)) +
    labs(x = NULL, y = "Post-designation log protected/reference ratio", color = NULL, shape = NULL) +
    theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = "bottom",
      axis.text.x = element_text(size = 9), strip.background = element_rect(fill = "grey96"))
  ggsave(file.path(root, "plots/fig_s_coverage_cohorts.png"), p, width = 7, height = 6.3, dpi = 300, bg = "white")
  ggsave(file.path(root, "plots/fig_s_coverage_cohorts.pdf"), p, width = 7, height = 6.3)
  stopifnot(identical(unname(hashes), unname(tools::md5sum(file.path(root, paths)))))
  out(data.frame(path = paths, md5 = unname(hashes)), "inputs")
  out(data.frame(status = "complete", n_primary_histories = nrow(series), n_before_eligible = sum(series$n_before >= 3 & series$n_after >= 3),
    n_after_ineligible = sum(!series$eligible_after)), "status")
  print(summary %>% filter(min_before == 3, source_subset == "All programmes", region_scope == "All") %>%
    select(Taxa, Resp, cohort, metric, n_mpa, estimate, ci_low, ci_high))
  invisible(list(series = series, roster = roster, pairs = pairs, summary = summary, difference = difference))
}
if (sys.nframe() == 0L) run_before_after_cohorts()
