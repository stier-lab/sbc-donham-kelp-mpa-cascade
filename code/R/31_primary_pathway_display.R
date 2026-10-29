# Display the six Table 3 regressions without changing primary estimates.
# Run with Rscript code/R/31_primary_pathway_display.R or source this file.
local({
  script_name <- "31_primary_pathway_display.R"
  source_paths <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
  command_paths <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  script_paths <- c(source_paths, command_paths)
  script_paths <- script_paths[basename(script_paths) == script_name]
  root <- if (length(script_paths)) {
    dirname(dirname(dirname(normalizePath(tail(script_paths, 1L), mustWork = TRUE))))
  } else normalizePath(getwd(), mustWork = TRUE)

  run_display <- function() {
    had_rng <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    original_rng <- if (had_rng) get(".Random.seed", envir = globalenv()) else NULL
    original_options <- options()
    on.exit({
      if (had_rng) assign(".Random.seed", original_rng, envir = globalenv())
      else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(list = ".Random.seed", envir = globalenv())
      }
      added_options <- setdiff(names(options()), names(original_options))
      if (length(added_options)) {
        options(stats::setNames(rep(list(NULL), length(added_options)), added_options))
      }
      options(original_options)
    }, add = TRUE)

    packages <- c("dplyr", "ggplot2", "ggrepel", "metafor", "ragg", "scales")
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) stop("Required packages unavailable: ", paste(missing, collapse = ", "))
    input_files <- c(primary = "data/sumstats_final.csv",
                     table3 = "tables/table_03_cross_taxa_meta_regression.csv")
    paths <- file.path(root, input_files)
    if (any(!file.exists(paths))) stop("Missing pathway input: ", paste(paths[!file.exists(paths)], collapse = ", "))
    input_md5 <- unname(tools::md5sum(paths))
    final <- utils::read.csv(paths[1], stringsAsFactors = FALSE)
    table3 <- utils::read.csv(paths[2], stringsAsFactors = FALSE)
    needed <- c("MPA", "Taxa", "Source", "Resp", "Mean", "SE", "Primary", "Location", "AnalysisType", "Model")
    if (!all(needed %in% names(final))) stop("Incomplete primary effect-size schema")
    if (!all(c("Relationship", "Predictor", "Response_Var", "k", "Intercept_Est", "Slope_Est", "Slope_SE",
               "Slope_CI_lower", "Slope_CI_upper", "Slope_pval") %in% names(table3))) {
      stop("Incomplete Table 3 schema")
    }
    keys <- c("MPA", "Taxa", "Source", "Resp")
    if (!nrow(final) || anyNA(final[keys]) || anyDuplicated(final[keys]) ||
        anyNA(final$Primary) || any(final$Primary != "Y") ||
        any(!is.finite(final$Mean)) || any(!is.finite(final$SE) | final$SE <= 0)) {
      stop("Primary rows must have unique complete keys, finite effects, positive SE and Primary=Y")
    }
    if (anyNA(final$Location) || any(!final$Location %in% c("C", "ML")) ||
        any(vapply(split(final$Location, final$MPA), function(x) length(unique(x)), integer(1)) != 1L)) {
      stop("Each MPA must have one valid island/mainland assignment")
    }
    final <- final[do.call(order, final[keys]), ]
    final$vi <- final$SE^2

    # Match script 09: equal source means, with Var(mean) = sum(vi) / n^2.
    pooled <- dplyr::summarise(dplyr::group_by(final, MPA, Taxa, Resp),
      Location = dplyr::first(Location),
      Sources = paste(Source, collapse = "; "),
      primary_keys = paste(paste(MPA, Taxa, Source, Resp, sep = "|"), collapse = "; "),
      AnalysisTypes = paste(AnalysisType, collapse = "; "),
      Models = paste(Model, collapse = "; "),
      n_sources = dplyr::n(), Mean = mean(Mean), vi = sum(vi) / dplyr::n()^2,
      .groups = "drop")
    pooled$SE <- sqrt(pooled$vi)
    definitions <- data.frame(
      panel = letters[1:6],
      x_Taxa = c("P. interruptus", "P. interruptus", "S. pulcher", "S. pulcher", "S. purpuratus", "M. franciscanus"),
      x_Resp = c("Bio", "Den", "Bio", "Den", "Den", "Den"),
      y_Taxa = c(rep("S. purpuratus", 4), rep("M. pyrifera", 2)),
      y_Resp = c("Bio", "Den", "Bio", "Den", "Bio", "Bio"),
      panel_name = c("(a) Lobster biomass /\npurple urchin biomass",
        "(b) Lobster density /\npurple urchin density",
        "(c) Sheephead biomass /\npurple urchin biomass",
        "(d) Sheephead density /\npurple urchin density",
        "(e) Purple urchin density /\nkelp biomass",
        "(f) Red urchin density /\nkelp biomass"), stringsAsFactors = FALSE)
    response_name <- c(Bio = "biomass", Den = "density")
    point_rows <- prediction_rows <- check_rows <- models <- vector("list", nrow(definitions))
    for (i in seq_len(nrow(definitions))) {
      def <- definitions[i, ]
      extract_side <- function(taxon, response, prefix) {
        x <- as.data.frame(pooled[pooled$Taxa == taxon & pooled$Resp == response, ])
        names(x)[names(x) != "MPA"] <- paste0(prefix, "_", names(x)[names(x) != "MPA"])
        x
      }
      pairs <- merge(extract_side(def$x_Taxa, def$x_Resp, "x"),
                     extract_side(def$y_Taxa, def$y_Resp, "y"), by = "MPA", sort = TRUE)
      if (nrow(pairs) < 3L || any(pairs$x_Location != pairs$y_Location)) {
        stop("Insufficient or inconsistent primary membership for panel ", def$panel)
      }
      predictor <- paste(def$x_Taxa, response_name[def$x_Resp])
      response <- paste(def$y_Taxa, response_name[def$y_Resp])
      expected <- table3[table3$Predictor == predictor & table3$Response_Var == response, ]
      if (nrow(expected) != 1L) stop("Table 3 relationship missing or duplicated: ", predictor, " -> ", response)
      model <- metafor::rma(yi = y_Mean, vi = y_vi, mods = ~ x_Mean,
        data = pairs, method = "REML", test = "knha")
      actual <- c(Intercept_Est = unname(model$b[1]), Slope_Est = unname(model$b[2]),
        Slope_SE = model$se[2], Slope_CI_lower = model$ci.lb[2], Slope_CI_upper = model$ci.ub[2])
      target <- unlist(expected[names(actual)], use.names = FALSE)
      if (model$k != expected$k || any(!is.finite(target)) ||
          any(abs(round(actual, 4) - target) > 1e-8) ||
          !is.finite(expected$Slope_pval) || abs(model$pval[2] - expected$Slope_pval) > 1e-8) {
        stop("Primary pathway verification failed before plotting: ", predictor, " -> ", response)
      }
      pairs$panel <- def$panel
      pairs$panel_name <- def$panel_name
      pairs$Relationship <- expected$Relationship
      pairs$Region <- ifelse(pairs$x_Location == "C", "Islands", "Mainland")
      pairs$k <- model$k
      pairs$verified_slope <- actual["Slope_Est"]
      pairs$table3_slope <- expected$Slope_Est
      pairs$primary_input_md5 <- input_md5[1]
      pairs$table3_input_md5 <- input_md5[2]
      point_rows[[i]] <- pairs
      x_grid <- seq(min(pairs$x_Mean), max(pairs$x_Mean), length.out = 200L)
      pred <- stats::predict(model, newmods = x_grid)
      prediction_rows[[i]] <- data.frame(panel_name = def$panel_name, x = x_grid,
        fitted = pred$pred, ci_low = pred$ci.lb, ci_high = pred$ci.ub)
      check_rows[[i]] <- data.frame(panel = def$panel, k = model$k,
        slope = actual["Slope_Est"], table3_slope = expected$Slope_Est,
        rounded_match = TRUE, p = model$pval[2], row.names = NULL)
      models[[i]] <- model
    }
    membership <- dplyr::bind_rows(point_rows)
    predictions <- dplyr::bind_rows(prediction_rows)
    checks <- dplyr::bind_rows(check_rows)
    if (!identical(input_md5, unname(tools::md5sum(paths)))) stop("Primary inputs changed during pathway verification")
    panel_levels <- definitions$panel_name
    membership$panel_name <- factor(membership$panel_name, levels = panel_levels)
    predictions$panel_name <- factor(predictions$panel_name, levels = panel_levels)
    counts <- unique(membership[c("panel_name", "k")])
    counts$label <- paste0("k = ", counts$k)
    focal_mpas <- c("Abalone Cove SMCA", "Campus Point SMCA", "Naples SMCA", "Point Vicente SMCA")
    labels <- membership[membership$panel == "a", ]
    if (!setequal(labels$MPA, focal_mpas)) stop("Focal four-MPA membership changed; update the caption deliberately")
    labels$label <- sub(" (SMCA|SMR)$", "\n\\1", labels$MPA)
    # Start the focal labels away from the fit and whiskers; repel them within the facet.
    label_dx <- c("Abalone Cove SMCA" = -0.15, "Campus Point SMCA" = 0.4,
                  "Naples SMCA" = 0.9, "Point Vicente SMCA" = -0.4)
    label_dy <- c("Abalone Cove SMCA" = 0.9, "Campus Point SMCA" = -0.45,
                  "Naples SMCA" = 1.6, "Point Vicente SMCA" = 0.9)

    width_mm <- 170; height_mm <- 210; dpi <- 600
    panel_width <- (width_mm - 16 - 4) / 2
    base_size <- if (panel_width < 50) 8 else if (panel_width < 80) 9 else 10
    region_colors <- c(Islands = "#0072B2", Mainland = "#D55E00")
    plot <- ggplot2::ggplot(membership, ggplot2::aes(x = x_Mean, y = y_Mean)) +
      ggplot2::geom_ribbon(data = predictions,
        ggplot2::aes(x = x, ymin = ci_low, ymax = ci_high), inherit.aes = FALSE,
        fill = "grey80", alpha = 0.5) +
      ggplot2::geom_hline(yintercept = 0, color = "grey65", linewidth = 0.25, linetype = "dashed") +
      ggplot2::geom_vline(xintercept = 0, color = "grey65", linewidth = 0.25, linetype = "dashed") +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = y_Mean - y_SE, ymax = y_Mean + y_SE),
        width = 0, color = "grey55", linewidth = 0.3) +
      ggplot2::geom_errorbar(ggplot2::aes(xmin = x_Mean - x_SE, xmax = x_Mean + x_SE),
        orientation = "y", width = 0, color = "grey55", linewidth = 0.3) +
      ggplot2::geom_line(data = predictions, ggplot2::aes(x = x, y = fitted),
        inherit.aes = FALSE, color = "grey20", linewidth = 0.55) +
      ggplot2::geom_point(ggplot2::aes(color = Region, shape = Region), size = 1.8) +
      ggrepel::geom_text_repel(data = labels, ggplot2::aes(label = label),
        size = 2.5, family = "sans", lineheight = 0.95, seed = 31L,
        nudge_x = unname(label_dx[labels$MPA]), nudge_y = unname(label_dy[labels$MPA]),
        max.overlaps = Inf, max.iter = 10000, max.time = Inf,
        box.padding = 0.5, point.padding = 0.3, min.segment.length = 0,
        segment.color = "grey45", segment.size = 0.25) +
      ggplot2::geom_text(data = counts, ggplot2::aes(x = Inf, y = Inf, label = label),
        inherit.aes = FALSE, hjust = 1.1, vjust = 1.3, size = 2.6) +
      ggplot2::facet_wrap(~ panel_name, ncol = 2, scales = "free", axes = "all") +
      ggplot2::scale_color_manual(values = region_colors) +
      ggplot2::scale_shape_manual(values = c(Islands = 16, Mainland = 17)) +
      ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 4),
        expand = ggplot2::expansion(mult = c(0.12, 0.18))) +
      ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(n = 4),
        expand = ggplot2::expansion(mult = c(0.12, 0.18))) +
      ggplot2::labs(x = "Consumer MPA effect (log contrast; unitless)",
        y = "Resource MPA effect (log contrast; unitless)", color = NULL, shape = NULL) +
      ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
      ggplot2::theme(strip.background = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(size = 8, face = "bold", hjust = 0),
        axis.text = ggplot2::element_text(size = 8, color = "grey20"),
        axis.title = ggplot2::element_text(size = 9),
        panel.spacing.x = grid::unit(4, "mm"), panel.spacing.y = grid::unit(4, "mm"),
        legend.position = "bottom", legend.key.width = grid::unit(5, "mm"),
        legend.margin = ggplot2::margin(0, 0, 0, 0),
        plot.margin = ggplot2::margin(3, 3, 3, 3, "mm"))

    outputs <- c(pdf = "plots/fig_s_primary_pathway.pdf", png = "plots/fig_s_primary_pathway.png",
      membership = "tables/table_s_primary_pathway_membership.csv", legend = "docs/primary_pathway_legend.md")
    for (directory in unique(dirname(file.path(root, outputs)))) dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
    ggplot2::ggsave(file.path(root, outputs["pdf"]), plot, device = pdf_device,
      width = width_mm, height = height_mm, units = "mm", bg = "white")
    ggplot2::ggsave(file.path(root, outputs["png"]), plot, device = ragg::agg_png,
      width = width_mm, height = height_mm, units = "mm", dpi = dpi, bg = "white")
    utils::write.csv(membership, file.path(root, outputs["membership"]), row.names = FALSE)
    writeLines(c("# Primary Consumer-Resource Pathways", "", "## Caption", "",
      "The six primary cross-taxa meta-regressions. Each point is one MPA, with primary source-specific effect sizes averaged equally within taxon and response; source units are never averaged. Panels relate (a) lobster biomass to purple urchin biomass, (b) lobster density to purple urchin density, (c) sheephead biomass to purple urchin biomass, (d) sheephead density to purple urchin density, (e) purple urchin density to kelp biomass and (f) red urchin density to kelp biomass. Lobster is Panulirus interruptus, sheephead is Semicossyphus pulcher, purple urchin is Strongylocentrotus purpuratus, red urchin is Mesocentrotus franciscanus and giant kelp is Macrocystis pyrifera.", "",
      "Axes show primary MPA effect sizes as unitless log contrasts, not physical biomass or density. Their estimand follows the selected source-specific model, as recorded in the membership export. For Mean CI histories, it is the mean log contrast over retained CI input years, which can include pre-designation and designation-year observations; these histories are not uniformly post-only. Blue circles denote islands; orange triangles denote mainland MPAs. Whiskers are +/-1 pooled sampling SE in both axes. Solid lines are metafor rma fits with REML and Knapp-Hartung tests; grey bands are nominal 95% confidence intervals for the mean relationship, not prediction intervals. Fits are displayed only within the observed predictor range; panel scales differ. k counts MPAs, not source rows or annual observations. The four points in panel (a) are labeled because only four MPAs contribute.", "",
      "For n contributing sources, the effect is their arithmetic mean and its sampling variance is sum(SE^2)/n^2, matching the primary Table 3 implementation. This arithmetic assumes independent source estimates. Predictor uncertainty is displayed but is not modeled by rma. Predictor and response can draw on different source sets; kelp includes satellite canopy as well as diver-derived biomass indices. These are estimated cross-sectional associations, not direct observations of predation or evidence that an MPA caused a trophic cascade.", "",
      "MPAs are not wholly independent: Naples SMCA and Campus Point SMCA share mapped PISCO Arroyo Quemado reference reefs, including in purple-urchin biomass. Programme footprints can overlap, and satellite reference polygon overlap is unverified. Equal-source averaging does not remove these dependencies. The primary regression and its nominal intervals are retained without spatial/shared-reference adjustment; no new primary estimates or significance claims are introduced.", "",
      "## Reproduction And Verification", "",
      "Run Rscript code/R/31_primary_pathway_display.R, or source the script; all working objects and RNG/options changes remain local or are restored. Both run_review_concerns.R and run_all.R include this module. Standalone execution loads no full pipeline, raw data or external workspace.", "",
      "Every fit must match Table 3 k, the intercept, slope, slope SE and slope CI endpoints to their four-decimal export rounding, plus the unrounded slope p-value within 1e-8, before any figure or table is written. A mismatch fails loudly rather than overwriting the primary table.", "",
      sprintf("Export dimensions: %d x %d mm; PNG %d dpi; six facets in two columns.", width_mm, height_mm, dpi), "",
      paste0("- Primary input: ", input_files[1], "; MD5 ", input_md5[1]),
      paste0("- Verification input: ", input_files[2], "; MD5 ", input_md5[2]),
      paste0("- Outputs: ", paste(outputs, collapse = ", ")), "",
      "The membership CSV has one row per plotted MPA/panel, with both taxa/responses, individual primary source keys, source lists, source models/designs, source counts, equal-source means, variances and SEs. verified_slope and table3_slope document concordance. All primary values are read-only.", "",
      "| Panel | k | Verified slope | Table 3 slope |", "|---|---:|---:|---:|",
      sprintf("| %s | %d | %.8f | %.4f |", checks$panel, checks$k, checks$slope, checks$table3_slope), "",
      paste0("Available runtime: ", R.version.string, "; ", paste(
        paste0(packages, " ", vapply(packages, function(p) as.character(utils::packageVersion(p)), character(1))),
        collapse = "; "), ".")), file.path(root, outputs["legend"]))
    if (!identical(input_md5, unname(tools::md5sum(paths)))) stop("Primary inputs changed while rendering")
    print(checks, row.names = FALSE)
    cat(sprintf("Primary pathways: six verified fits; %d MPA/panel points; %d unique MPAs; %d x %d mm at %d dpi.\n",
      nrow(membership), length(unique(membership$MPA)), width_mm, height_mm, dpi))
    cat(paste(outputs, collapse = "\n"), "\n")
    invisible(list(membership = membership, checks = checks, plot = plot, models = models, outputs = outputs))
  }
  run_display()
})
