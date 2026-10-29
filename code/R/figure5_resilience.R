# Display and coverage exports for the descriptive and paired resilience panels.
build_figure5_resilience <- function(kelp, pair_tab, trajectory) {
  annual <- data.table::dcast(
    kelp[year >= 2008 & year <= 2023],
    CA_MPA_Name_Short + year ~ status, value.var = "value")
  grid <- data.table::CJ(CA_MPA_Name_Short = sort(unique(annual$CA_MPA_Name_Short)),
                        year = 2008:2023)
  coverage <- merge(grid, annual, all.x = TRUE,
                    by = c("CA_MPA_Name_Short", "year"))
  coverage[, `:=`(mpa_present = is.finite(mpa),
                  reference_present = is.finite(reference))]
  coverage[, paired_present := mpa_present & reference_present]
  counts <- coverage[, .(n_mpa = sum(mpa_present),
                         n_reference = sum(reference_present),
                         n_pairs = sum(paired_present),
                         contributing_pairs = paste(CA_MPA_Name_Short[paired_present],
                                                    collapse = "; ")), by = year]
  stopifnot(all(coverage$mpa_present == coverage$reference_present))
  write.csv(coverage, here::here("tables", "table_figure5_annual_membership.csv"),
            row.names = FALSE)
  write.csv(counts, here::here("tables", "table_figure5_annual_coverage.csv"),
            row.names = FALSE)

  pairs <- data.table::copy(pair_tab)
  periods <- c("2014-2016", "2017-2019", "2020-2023")
  pairs[, ratio := inside / outside]
  pairs[, period_position := match(period, rev(periods))]
  # Deterministic offsets separate individual observations, not summary estimates.
  pairs[, offset := seq(-0.25, 0.05, length.out = .N),
        by = period]
  pairs[, direction := factor(ifelse(ratio > 1, "MPA", "Reference"),
                              levels = c("MPA", "Reference"))]
  summary <- pairs[, {
    interval <- stats::t.test(log(ratio))$conf.int
    list(n_pairs = .N, ratio = exp(mean(log(ratio))),
         lower = exp(interval[1]), upper = exp(interval[2]))
  }, by = .(metric, period, period_position)]
  # These are the same geometric means and log-ratio t intervals as the
  # baseline sensitivity table, not confidence intervals for the signed-rank test.
  write.csv(pairs, here::here("tables", "table_figure5_period_pairs.csv"), row.names = FALSE)
  write.csv(summary, here::here("tables", "table_figure5_period_summary.csv"), row.names = FALSE)

  paired_plot <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey55",
                       linewidth = 0.4) +
    ggplot2::geom_point(data = pairs,
      ggplot2::aes(ratio, period_position + offset, color = direction), size = 1.8,
      alpha = 0.85) +
    ggplot2::geom_errorbar(data = summary,
      ggplot2::aes(xmin = lower, xmax = upper, y = period_position + 0.20),
      orientation = "y", width = 0.13, linewidth = 0.65) +
    ggplot2::geom_point(data = summary,
      ggplot2::aes(ratio, period_position + 0.20), shape = 23, size = 2.8,
      fill = "white", color = "black", stroke = 0.7) +
    ggplot2::scale_color_manual(values = c(MPA = "#0072B2", Reference = "#D55E00"),
                               guide = "none") +
    ggplot2::scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 50, 100, 300),
                         labels = c("0.5", "1", "2", "5", "10", "50", "100", "300"),
                         expand = ggplot2::expansion(mult = c(0.06, 0.06))) +
    ggplot2::scale_y_continuous(breaks = 1:3, labels = rev(periods),
                               limits = c(0.5, 3.5), expand = c(0, 0)) +
    ggplot2::labs(x = "Baseline-relative biomass ratio\n(MPA / reference)", y = NULL) +
    theme_mpa(base_size = 9) +
    ggplot2::theme(plot.margin = ggplot2::margin(6, 4, 4, 6))
  plot <- patchwork::wrap_plots(trajectory, paired_plot, nrow = 1) +
    patchwork::plot_annotation(tag_levels = "A") &
    ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 11))
  list(plot = plot, coverage = coverage, counts = counts, pairs = pairs, summary = summary)
}
