size_refit_display <- function(root) {
  read <- function(stem) read.csv(file.path(root, "tables", paste0(stem, ".csv")))
  label <- function(d) d %>% mutate(Taxa = factor(Taxa,
    levels = c("S. pulcher", "P. interruptus"), labels = c("Sheephead", "Spiny lobster")),
    Resp = factor(Resp, levels = c("Den", "Bio"), labels = c("Density", "Biomass")),
    size = factor(size, levels = c("all", "large"), labels = c("All sizes", "Large individuals")))
  palette <- c("All sizes" = "#697D85", "Large individuals" = "#0072B2")
  common <- list(ggplot2::scale_colour_manual(values = palette), ggplot2::scale_fill_manual(values = palette),
    ggplot2::scale_linetype_manual(values = c("All sizes" = "dashed", "Large individuals" = "solid")),
    ggplot2::theme_classic(base_size = 10, base_family = "Times New Roman"),
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold"),
      panel.spacing = grid::unit(5, "mm")),
    ggplot2::guides(fill = "none", linetype = "none"))
  save <- function(plot, stem) {
    ggplot2::ggsave(file.path(root, "plots", paste0(stem, ".png")), plot,
      width = 180, height = 150, units = "mm", dpi = 300, bg = "white")
    ggplot2::ggsave(file.path(root, "plots", paste0(stem, ".pdf")), plot,
      width = 180, height = 150, units = "mm", device = grDevices::cairo_pdf, bg = "white")
  }
  annual <- read("table_s_size_refit_temporal_annual") %>%
    filter(scope == "expanded_size_cohort", Taxa %in% c("P. interruptus", "S. pulcher"), time <= 15) %>% label()
  prediction <- read("table_s_size_refit_temporal_predictions") %>%
    filter(scope == "expanded_size_cohort", Species %in% c("P. interruptus", "S. pulcher"), model == "calendar_AR1") %>%
    rename(Taxa = Species, Resp = resp) %>% label()
  support <- annual %>% group_by(Taxa, Resp) %>% summarise(first = min(time), last = max(time), .groups = "drop")
  prediction <- inner_join(prediction, support, by = c("Taxa", "Resp")) %>% filter(time >= first, time <= last)
  p <- ggplot2::ggplot(annual, ggplot2::aes(time, lnDiff, colour = size)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey75", linewidth = .3) +
    ggplot2::geom_point(alpha = .20, size = .9) +
    ggplot2::geom_ribbon(data = prediction, ggplot2::aes(y = NULL, ymin = ci_low, ymax = ci_high, fill = size),
      colour = NA, alpha = .13, inherit.aes = TRUE) +
    ggplot2::geom_line(data = prediction, ggplot2::aes(y = estimate, linetype = size), linewidth = .7) +
    ggplot2::facet_grid(Taxa ~ Resp) + common +
    ggplot2::labs(x = "Years since MPA designation", y = "Log MPA/reference ratio")
  save(p, "fig_s_size_refit_temporal")

  effects <- read("table_s_size_refit_pathway_effects") %>%
    filter(scope == "expanded_size_cohort", estimand == "post_only_mean")
  summary <- read("table_s_size_refit_pathways")
  grids <- effects %>% group_by(Taxa, Resp, size) %>% group_modify(~ {
    m <- metafor::rma(y, vi = vi, mods = ~x, data = .x, method = "REML", test = "knha")
    expected <- summary %>% filter(scope == "expanded_size_cohort", estimand == "post_only_mean",
      Taxa == .y$Taxa, Resp == .y$Resp, size == .y$size)
    stopifnot(nrow(expected) == 1, abs(as.numeric(coef(m)[2]) - expected$slope) < 1e-10)
    x <- seq(min(.x$x), max(.x$x), length.out = 100)
    pred <- predict(m, newmods = x)
    data.frame(x, fit = pred$pred, ci_low = pred$ci.lb, ci_high = pred$ci.ub)
  }) %>% ungroup() %>% label()
  p <- ggplot2::ggplot(label(effects), ggplot2::aes(x, y, colour = size)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey80", linewidth = .3) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey80", linewidth = .3) +
    ggplot2::geom_ribbon(data = grids, ggplot2::aes(y = NULL, ymin = ci_low, ymax = ci_high, fill = size),
      colour = NA, alpha = .13) +
    ggplot2::geom_line(data = grids, ggplot2::aes(y = fit, linetype = size), linewidth = .7) +
    ggplot2::geom_point(size = 1.7, alpha = .85) + ggplot2::facet_grid(Taxa ~ Resp, scales = "free") + common +
    ggplot2::labs(x = "Predator mean log MPA/reference ratio", y = "Purple urchin mean log MPA/reference ratio")
  save(p, "fig_s_size_refit_pathways")
  writeLines(c("# Size-specific refit figures", "",
    "Temporal display: points are observed annual PISCO log ratios, with all/large values matched by MPA and year. Lines and nominal 95% confidence bands are from the calendar-AR1 fits, not the independent-error fits. Density and biomass use different valid measurement coverage. Fitted lines are clipped to the observed time range and the first 15 years. Kelp/urchin observations also enter the joint model but are not displayed here.", "",
    "Pathway display: points are MPA mean post-designation log ratios over the same paired years for each predator/urchin comparison. Lines and nominal 95% confidence bands use REML meta-regressions with Knapp-Hartung inference. There are ten MPAs for each density relationship and eight for biomass. These are level comparisons, not designation-induced changes. No mean-only large-predator association survives correction across the four links. The legacy mixed-route alternatives are reported separately in the tables.", "",
    "Large sheephead exceed 30 cm total length; large lobsters are at least 82.55 mm carapace length. Lobster values are a secondary size-pool reconstruction, not reconciled replacement primary biomass. All intervals are conditional on the reconstructed annual measurements. Size-pool/moderator uncertainty and covariance from shared references are not fully propagated. These displays do not establish consumption or a causal pathway.", "",
    "Drivers: code/R/41_predator_size_refits.R and predator_size_refit_display.R. Inputs: tables/table_s_size_refit_*.csv."),
    file.path(root, "docs/size_refit_figure_legends.md"))
}
