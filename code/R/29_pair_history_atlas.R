# Annual paired histories for exactly the primary effect-size keys.
# Run with Rscript code/R/29_pair_history_atlas.R or source this file.
# All objects remain local; no pipeline scripts or raw data are loaded.

local({
  script_name <- "29_pair_history_atlas.R"
  source_paths <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
  command_paths <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))
  script_paths <- c(source_paths, command_paths)
  script_paths <- script_paths[basename(script_paths) == script_name]
  if (length(script_paths)) {
    script_path <- normalizePath(tail(script_paths, 1L), mustWork = TRUE)
    root <- dirname(dirname(dirname(script_path)))
  } else {
    root <- normalizePath(getwd(), mustWork = TRUE)
  }

  run_atlas <- function() {
    source(file.path(root, "code/R/taxon_display.R"), local = TRUE)
    # Restore session state, including an initially absent RNG, after rendering.
    had_rng <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    original_rng <- if (had_rng) get(".Random.seed", envir = globalenv()) else NULL
    original_options <- options()
    on.exit({
      if (had_rng) {
        assign(".Random.seed", original_rng, envir = globalenv())
      } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(list = ".Random.seed", envir = globalenv())
      }
      added_options <- setdiff(names(options()), names(original_options))
      if (length(added_options)) {
        options(stats::setNames(rep(list(NULL), length(added_options)), added_options))
      }
      options(original_options)
    }, add = TRUE)
    required_packages <- c("ggplot2", "ggrepel", "ragg")
    missing_packages <- required_packages[!vapply(
      required_packages, requireNamespace, logical(1), quietly = TRUE
    )]
    if (length(missing_packages)) {
      stop("Required packages unavailable: ", paste(missing_packages, collapse = ", "))
    }

    input_paths <- c(
      responses = "data/harmonized/harmonized_response_ratios.csv",
      landsat = "data/harmonized/harmonized_landsat_rr.csv",
      metadata = "data/harmonized/harmonized_site_metadata.csv",
      final = "data/sumstats_final.csv"
    )
    absolute_inputs <- file.path(root, input_paths)
    if (any(!file.exists(absolute_inputs))) {
      stop("Atlas inputs missing: ", paste(input_paths[!file.exists(absolute_inputs)],
                                           collapse = ", "))
    }
    inputs <- lapply(absolute_inputs, utils::read.csv, stringsAsFactors = FALSE)
    names(inputs) <- names(input_paths)
    input_md5 <- unname(tools::md5sum(absolute_inputs))

    require_columns <- function(data, columns, label) {
      missing <- setdiff(columns, names(data))
      if (length(missing)) {
        stop(label, " missing columns: ", paste(missing, collapse = ", "))
      }
    }
    key_columns <- c("MPA", "Taxa", "Source", "Resp")
    make_key <- function(data) do.call(paste, c(data[key_columns], sep = "|"))
    taxa_map <- c(
      "Panulirus interruptus" = "P. interruptus",
      "Semicossyphus pulcher" = "S. pulcher",
      "Strongylocentrotus purpuratus" = "S. purpuratus",
      "Mesocentrotus franciscanus" = "M. franciscanus",
      "Macrocystis pyrifera" = "M. pyrifera"
    )
    response_labels <- c(Bio = "Biomass", Den = "Density")
    source_order <- c("KFM", "LTER", "PISCO", "Landsat")

    standardize_pairs <- function(data, input_file) {
      require_columns(data, c("CA_MPA_Name_Short", "year", "y", "source", "resp",
                              "mpa", "reference", "BA"), input_file)
      out <- data.frame(
        MPA = data$CA_MPA_Name_Short, Taxa = unname(taxa_map[data$y]),
        Source = data$source, Resp = data$resp, year = data$year,
        mpa = data$mpa, reference = data$reference, BA = data$BA,
        input_file = input_file, dataset_year_min = min(data$year),
        dataset_year_max = max(data$year), stringsAsFactors = FALSE
      )
      if (anyNA(out$Taxa)) stop("Unmapped harmonized taxon in ", input_file)
      if (anyNA(out$year) || any(!is.finite(out$year)) ||
          any(out$year != as.integer(out$year))) {
        stop("Annual years must be finite integers in ", input_file)
      }
      out$key <- make_key(out)
      out
    }
    pairs <- rbind(
      standardize_pairs(inputs$responses, input_paths[["responses"]]),
      standardize_pairs(inputs$landsat, input_paths[["landsat"]])
    )

    final <- inputs$final
    require_columns(final, c(key_columns, "Primary", "N", "AnalysisType", "Model",
                             "MPA_Start"), input_paths[["final"]])
    final <- final[!is.na(final$Primary) & final$Primary == "Y", , drop = FALSE]
    if (!nrow(final)) stop("No primary series in the final effect-size table")
    if (anyNA(final[key_columns]) || any(duplicated(make_key(final)))) {
      stop("Final primary keys must be complete and unique")
    }
    if (any(!final$Taxa %in% unname(taxa_map)) ||
        any(!final$Source %in% source_order) ||
        any(!final$Resp %in% names(response_labels))) {
      stop("Unexpected taxon, source or response in final primary keys")
    }
    final$key <- make_key(final)
    unmatched <- final[!final$key %in% pairs$key, key_columns, drop = FALSE]
    if (nrow(unmatched)) {
      stop("Final primary series missing from harmonized paired inputs:\n",
           paste(make_key(unmatched), collapse = "\n"))
    }
    excluded_key_count <- length(setdiff(unique(pairs$key), final$key))
    pairs <- pairs[pairs$key %in% final$key, , drop = FALSE]
    if (any(duplicated(pairs[c("key", "year")]))) {
      stop("Duplicate annual source/taxon/response keys; cannot average source units")
    }
    if (any(is.infinite(pairs$mpa) | is.infinite(pairs$reference)) ||
        any(c(pairs$mpa, pairs$reference) < 0, na.rm = TRUE)) {
      stop("Invalid harmonized proportions: negative or infinite values")
    }

    metadata <- inputs$metadata
    require_columns(metadata, c("CA_MPA_Name_Short", "MPA_Start"),
                    input_paths[["metadata"]])
    metadata <- unique(metadata[c("CA_MPA_Name_Short", "MPA_Start")])
    if (anyDuplicated(metadata$CA_MPA_Name_Short)) {
      stop("Ambiguous metadata designation year")
    }
    final <- final[order(final$MPA, match(final$Source, source_order),
                         match(final$Resp, names(response_labels)),
                         match(final$Taxa, unname(taxa_map))), , drop = FALSE]
    rownames(final) <- NULL
    designation <- metadata$MPA_Start[match(final$MPA, metadata$CA_MPA_Name_Short)]
    if (anyNA(designation) || any(!is.finite(designation)) ||
        any(designation != as.integer(designation))) {
      stop("Missing or noninteger metadata designation year for a final MPA")
    }
    if (anyNA(final$MPA_Start) || any(final$MPA_Start != designation)) {
      stop("Final table and site metadata disagree on MPA_Start")
    }
    final$MPA_Start <- designation
    pairs$MPA_Start <- final$MPA_Start[match(pairs$key, final$key)]
    expected_ba <- ifelse(pairs$year <= pairs$MPA_Start, "Before", "After")
    if (anyNA(pairs$BA) || any(pairs$BA != expected_ba)) {
      stop("Harmonized BA conflicts with metadata (designation year belongs to Before)")
    }

    width_mm <- 170
    height_mm <- 200
    dpi <- 300
    panels_per_page <- 4L
    calc_layout <- function(n_cols, n_rows, width, height) {
      panel_w <- (width - 16 - (n_cols - 1) * 4) / n_cols
      panel_h <- (height - 16 - (n_rows - 1) * 4 - 12) / n_rows
      list(panel_w = panel_w, panel_h = panel_h,
           base_size = if (panel_w < 50) 8 else if (panel_w < 80) 9 else 10,
           point_size = if (panel_w < 50) 1 else 1.5,
           line_width = if (panel_w < 50) 0.4 else 0.6)
    }
    layout <- calc_layout(1, 4, width_mm, height_mm)
    message(sprintf("Atlas layout: four stacked panels, %.0f x %.0f mm each",
                    layout$panel_w, layout$panel_h))
    palette <- c(MPA = "#0072B2", Reference = "#707070")
    font <- "sans"
    if (requireNamespace("systemfonts", quietly = TRUE) &&
        "Arial" %in% systemfonts::system_fonts()$family) font <- "Arial"
    year_range <- range(pairs$year)
    y_upper <- ceiling(max(c(pairs$mpa, pairs$reference), na.rm = TRUE) / 0.25) * 0.25
    y_upper <- max(y_upper, 1)
    x_limits <- c(year_range[1] - 0.6, year_range[2] + 6.5)
    y_limits <- c(-0.025 * y_upper, 1.025 * y_upper)

    final$page <- (seq_len(nrow(final)) - 1L) %/% panels_per_page + 1L
    final$panel <- LETTERS[(seq_len(nrow(final)) - 1L) %% panels_per_page + 1L]
    final$facet_label <- paste0(final$panel, "  ", final$MPA, " | ", final$Source,
                               " | ", taxon_display(final$Taxa), " | ", response_labels[final$Resp])
    final$png <- sprintf("plots/pair_history_atlas/page_%03d.png", final$page)
    final$pdf <- sprintf("plots/pair_history_atlas/page_%03d.pdf", final$page)

    manifest_rows <- lapply(seq_len(nrow(final)), function(i) {
      d <- pairs[pairs$key == final$key[i], , drop = FALSE]
      d <- d[order(d$year), , drop = FALSE]
      valid_pair <- is.finite(d$mpa) & is.finite(d$reference)
      absent_years <- setdiff(seq.int(min(d$year), max(d$year)), d$year)
      data.frame(
        page = final$page[i], panel = final$panel[i], png = final$png[i],
        pdf = final$pdf[i], MPA = final$MPA[i], Source = final$Source[i],
        Taxa = final$Taxa[i], Resp = final$Resp[i],
        response = unname(response_labels[final$Resp[i]]),
        taxon_full = names(taxa_map)[match(final$Taxa[i], taxa_map)],
        MPA_Start = final$MPA_Start[i],
        dataset = unique(d$input_file),
        dataset_year_min = unique(d$dataset_year_min),
        dataset_year_max = unique(d$dataset_year_max),
        year_min = min(d$year), year_max = max(d$year),
        n_annual_rows = nrow(d), n_pairs = sum(valid_pair),
        n_before = sum(valid_pair & d$year <= final$MPA_Start[i]),
        n_after = sum(valid_pair & d$year > final$MPA_Start[i]),
        n_missing_years = length(absent_years),
        missing_years = paste(absent_years, collapse = ";"),
        n_missing_mpa = sum(!is.finite(d$mpa)),
        n_missing_reference = sum(!is.finite(d$reference)),
        final_N = final$N[i], AnalysisType = final$AnalysisType[i],
        Model = final$Model[i], coverage_status = "matched",
        y_axis_min = 0, y_axis_max = y_upper,
        width_mm = width_mm, height_mm = height_mm, dpi = dpi,
        stringsAsFactors = FALSE
      )
    })
    manifest <- do.call(rbind, manifest_rows)
    if (any(manifest$n_pairs == 0L)) {
      stop("Final primary series has no finite annual pair:\n",
           paste(make_key(manifest[manifest$n_pairs == 0L, ]), collapse = "\n"))
    }
    if (anyNA(manifest$final_N) || any(manifest$final_N != manifest$n_pairs)) {
      stop("Final N differs from harmonized annual pairs:\n",
           paste(make_key(manifest[is.na(manifest$final_N) |
                                     manifest$final_N != manifest$n_pairs, ]),
                 collapse = "\n"))
    }
    if (length(unique(manifest$MPA)) < 21L) {
      stop("Coverage below the requested minimum of 21 MPAs")
    }

    # Complete each annual grid before drawing; only adjacent finite years connect.
    long_rows <- lapply(seq_len(nrow(final)), function(i) {
      d <- pairs[pairs$key == final$key[i], , drop = FALSE]
      years <- seq.int(min(d$year), max(d$year))
      d <- d[match(years, d$year), , drop = FALSE]
      do.call(rbind, lapply(c("MPA", "Reference"), function(status) {
        value <- if (status == "MPA") d$mpa else d$reference
        valid <- is.finite(value)
        segment <- cumsum(!valid | c(TRUE, !head(valid, -1)))
        data.frame(key = final$key[i], page = final$page[i],
                   facet_label = final$facet_label[i], year = years,
                   value = value, status = status,
                   segment = paste(i, status, segment, sep = "|"),
                   stringsAsFactors = FALSE)
      }))
    })
    long <- do.call(rbind, long_rows)
    long$status <- factor(long$status, levels = c("MPA", "Reference"))
    observed <- long[is.finite(long$value), , drop = FALSE]
    segment_counts <- table(observed$segment)
    lines <- observed[segment_counts[observed$segment] >= 2L, , drop = FALSE]
    ends <- do.call(rbind, lapply(split(observed, interaction(observed$key,
                                                            observed$status, drop = TRUE)),
                                 function(d) d[which.max(d$year), , drop = FALSE]))

    build_page <- function(page, colours = palette) {
      keys <- final[final$page == page, , drop = FALSE]
      with_facets <- function(d) {
        d <- d[d$page == page, , drop = FALSE]
        d$facet_label <- factor(d$facet_label, levels = keys$facet_label)
        d
      }
      points_page <- with_facets(observed)
      lines_page <- with_facets(lines)
      labels_page <- with_facets(ends)
      designation_page <- data.frame(
        facet_label = factor(keys$facet_label, levels = keys$facet_label),
        MPA_Start = keys$MPA_Start
      )
      ggplot2::ggplot(points_page, ggplot2::aes(year, value)) +
        ggplot2::annotate("rect", xmin = c(1996.5, 2013.5),
                          xmax = c(1998.5, 2016.5), ymin = -Inf, ymax = Inf,
                          fill = "#D5B96A", alpha = 0.13) +
        ggplot2::geom_vline(data = designation_page,
                            ggplot2::aes(xintercept = MPA_Start),
                            colour = "#555555", linetype = "dotdash", linewidth = 0.35) +
        ggplot2::geom_line(data = lines_page,
                           ggplot2::aes(colour = status, linetype = status,
                                        group = segment), linewidth = layout$line_width) +
        ggplot2::geom_point(ggplot2::aes(colour = status, shape = status),
                            size = layout$point_size, stroke = 0.25) +
        ggrepel::geom_label_repel(
          data = labels_page, ggplot2::aes(label = status, colour = status),
          nudge_x = 1.2, direction = "y", hjust = 0, seed = 29,
          size = 2.8, family = font, box.padding = 0.45, point.padding = 0.15,
          force = 2, force_pull = 0.3,
          fill = "white", linewidth = 0,
          label.padding = grid::unit(0.08, "lines"), label.r = grid::unit(0, "lines"),
          max.overlaps = Inf, min.segment.length = 0, segment.size = 0.2,
          segment.linetype = "dotted", segment.alpha = 0.6,
          ylim = c(0.025 * y_upper, 0.975 * y_upper),
          xlim = c(year_range[2] - 0.5, x_limits[2] - 0.2)
        ) +
        ggplot2::facet_wrap(~facet_label, ncol = 1, scales = "fixed") +
        ggplot2::scale_colour_manual(values = colours, guide = "none") +
        ggplot2::scale_linetype_manual(values = c(MPA = "solid", Reference = "dashed"),
                                       guide = "none") +
        ggplot2::scale_shape_manual(values = c(MPA = 16, Reference = 17), guide = "none") +
        ggplot2::scale_x_continuous(limits = x_limits,
                                    breaks = c(seq(1985, year_range[2], by = 10),
                                               year_range[2]),
                                    expand = ggplot2::expansion(mult = 0)) +
        ggplot2::scale_y_continuous(limits = y_limits, breaks = seq(0, y_upper, by = 0.5),
                                    expand = ggplot2::expansion(mult = 0)) +
        ggplot2::labs(x = "Year", y = "Normalized response (shared-scale proportion)") +
        ggplot2::theme_minimal(base_size = layout$base_size, base_family = font) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.major.y = ggplot2::element_line(colour = "#E7E7E7", linewidth = 0.2),
          axis.line.x = ggplot2::element_line(colour = "#777777", linewidth = 0.25),
          axis.ticks = ggplot2::element_line(colour = "#777777", linewidth = 0.25),
          axis.text = ggplot2::element_text(size = 8, colour = "#333333"),
          axis.title = ggplot2::element_text(size = 10, colour = "#222222"),
          axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 7)),
          strip.background = ggplot2::element_blank(),
          strip.text = ggplot2::element_text(size = 9, face = "bold", hjust = 0,
                                             margin = ggplot2::margin(t = 3, b = 4)),
          panel.spacing.y = grid::unit(4, "mm"),
          plot.margin = ggplot2::margin(4, 5, 4, 5, "mm"),
          legend.position = "none", plot.background = ggplot2::element_rect(
            fill = "white", colour = NA
          )
        )
    }

    output_dir <- file.path(root, "plots", "pair_history_atlas")
    qa_dir <- file.path(output_dir, "qa")
    dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(root, "tables"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(root, "docs"), recursive = TRUE, showWarnings = FALSE)
    pages <- sort(unique(final$page))
    save_png <- function(plot, path) {
      ggplot2::ggsave(path, plot = plot, device = ragg::agg_png,
                      width = width_mm, height = height_mm, units = "mm", dpi = dpi,
                      bg = "white")
    }
    open_pdf <- function(path, onefile) {
      if (capabilities("cairo")) {
        grDevices::cairo_pdf(path, width = width_mm / 25.4, height = height_mm / 25.4,
                             family = font, onefile = onefile, bg = "white")
      } else {
        grDevices::pdf(path, width = width_mm / 25.4, height = height_mm / 25.4,
                       family = "Helvetica", onefile = onefile, bg = "white")
      }
    }
    combined_pdf <- file.path(root, "plots", "fig_s_pair_history_atlas.pdf")
    open_pdf(combined_pdf, TRUE)
    combined_device <- grDevices::dev.cur()
    on.exit({
      if (combined_device %in% grDevices::dev.list()) grDevices::dev.off(combined_device)
    }, add = TRUE)
    for (page in pages) {
      plot <- build_page(page)
      if (any(vapply(c("title", "subtitle"), function(label) {
        !is.null(plot$labels[[label]])
      }, logical(1)))) stop("Atlas figures must have no title or subtitle")
      png_path <- file.path(root, unique(final$png[final$page == page]))
      pdf_path <- file.path(root, unique(final$pdf[final$page == page]))
      save_png(plot, png_path)
      ggplot2::ggsave(pdf_path, plot = plot,
                      device = if (capabilities("cairo")) grDevices::cairo_pdf else "pdf",
                      width = width_mm, height = height_mm, units = "mm", bg = "white")
      grDevices::dev.set(combined_device)
      print(plot)
      message(sprintf("Rendered page %03d/%03d", page, length(pages)))
    }
    grDevices::dev.off(combined_device)

    manifest_path <- file.path(root, "tables", "table_s_pair_history_atlas_index.csv")
    utils::write.csv(manifest, manifest_path, row.names = FALSE, na = "")

    # QA assets live in a subdirectory and are deliberately absent from the atlas index.
    if (requireNamespace("magick", quietly = TRUE)) {
      page_paths <- file.path(root, unique(manifest$png))
      for (sheet in seq_len(ceiling(length(pages) / 6))) {
        indices <- seq.int((sheet - 1L) * 6L + 1L, min(sheet * 6L, length(pages)))
        thumbs <- lapply(indices, function(i) {
          magick::image_resize(magick::image_read(page_paths[i]), "480x")
        })
        rows <- lapply(split(seq_along(thumbs), ceiling(seq_along(thumbs) / 3)),
                       function(i) magick::image_append(do.call(c, thumbs[i])))
        sheet_image <- magick::image_append(do.call(c, rows), stack = TRUE)
        magick::image_write(sheet_image,
                            file.path(qa_dir, sprintf("contact_sheet_%02d.png", sheet)))
      }
    } else {
      warning("magick unavailable; contact-sheet QA was not generated")
    }
    if (requireNamespace("colorspace", quietly = TRUE)) {
      for (mode in c("deutan", "protan", "desaturate")) {
        simulate <- getExportedValue("colorspace", mode)
        cvd_palette <- stats::setNames(simulate(unname(palette)), names(palette))
        save_png(build_page(pages[1], cvd_palette),
                 file.path(qa_dir, paste0("page_001_", mode, ".png")))
      }
    }

    source_counts <- as.data.frame(table(manifest$Source, manifest$Resp))
    source_counts <- source_counts[source_counts$Freq > 0, , drop = FALSE]
    names(source_counts) <- c("Source", "Resp", "n")
    count_lines <- sprintf("| %s | %s | %d |", source_counts$Source,
                           response_labels[as.character(source_counts$Resp)], source_counts$n)
    n_mpa <- length(unique(manifest$MPA))
    legend <- c(
      "# Pair-history atlas legend", "",
      sprintf(paste0("Annual paired histories for the %d primary MPA/source/taxon/response ",
                     "series in `data/sumstats_final.csv`, covering all %d final MPAs. ",
                     "Each of %d pages contains four facets (A-D, top to bottom). ",
                     "Facet strips identify the full MPA name, source, abbreviated taxon ",
                     "and biomass or density response."), nrow(manifest), n_mpa, length(pages)),
      "",
      paste0("Blue (#0072B2) solid lines with circles denote MPA observations; grey ",
             "(#707070) dashed lines with triangles denote paired reference observations. ",
             "Direct labels identify the last observed values; thin dotted label leaders ",
             "are annotations, not additional observations. Lines connect consecutive ",
             "observed annual values only. Missing years and missing values break the ",
             "lines; isolated observations remain visible as markers."),
      "",
      paste0("Vertical grey dot-dash lines mark each MPA's actual `MPA_Start` from ",
             "`data/harmonized/harmonized_site_metadata.csv`. Subtle gold shading marks ",
             "the complete 1997-1998 and 2014-2016 calendar-year windows (1996.5-1998.5 ",
             "and 2013.5-2016.5 on the annual axis). These windows provide temporal ",
             "context; they are not estimates of local disturbance exposure."),
      "",
      sprintf(paste0("The vertical axis shows the actual `mpa` and `reference` proportions ",
                     "supplied in the harmonized paired-response tables, on their shared ",
                     "upstream normalization. Values are plotted unchanged, with the same ",
                     "0-%.2f axis across all facets and pages; no MPA/reference rescaling, ",
                     "source averaging, physical-unit reconstruction, smoothing or model ",
                     "fit is applied. Values above 1 in the harmonized files are retained. ",
                     "These are normalized responses, not absolute biomass or density, ",
                     "community fractions, or estimates of among-program absolute abundance. ",
                     "The supplied tables have no annual uncertainty estimates, so no ",
                     "confidence ribbons are inferred."), y_upper),
      "",
      "## Coverage and index", "",
      sprintf(paste0("All %d final primary keys have matched annual paired observations; ",
                     "there are no unmatched final series. The selected histories have ",
                     "%d missing MPA values and %d missing reference values. ",
                     "The atlas displays %d paired annual rows (%d individual values), ",
                     "spanning %d-%d. Annual paired counts match the final effect-size `N` ",
                     "for every series. %d series contain internal missing calendar years."),
              nrow(manifest), sum(manifest$n_missing_mpa), sum(manifest$n_missing_reference),
              sum(manifest$n_pairs), nrow(observed),
              year_range[1], year_range[2], sum(manifest$n_missing_years > 0)),
      "",
      sprintf(paste0("Exactly %d harmonized series whose full key is absent from the final ",
                     "primary set are omitted, including rejected series and other responses. ",
                     "No raw-response or external data files are read."), excluded_key_count),
      "",
      sprintf(paste0("%d selected histories have no Before-period observations (year <= ",
                     "`MPA_Start`) in the supplied paired inputs. Their facets therefore ",
                     "show only the observed later history; no baseline is imputed."),
              sum(manifest$n_before == 0L)),
      "",
      "| Source | Response | Series |", "| --- | --- | ---: |", count_lines, "",
      paste0("`tables/table_s_pair_history_atlas_index.csv` has one row per facet, keyed ",
             "by `page` and `panel`. `png` and `pdf` are repository-relative paths to ",
             "individual atlas pages. Include each distinct `png` once, sorted by `page`; ",
             "QA PNGs in `plots/pair_history_atlas/qa/` are not atlas pages. Pages follow ",
             "MPA name, source (KFM, LTER, PISCO, Landsat), response (Bio, Den) and ",
             "taxon (predators through kelp) order."),
      "",
      paste0("`dataset_year_min/max` describe the full input-file range; `year_min/max` ",
             "describe each selected paired history. `n_before` counts finite annual pairs ",
             "with year <= `MPA_Start`; `n_after` counts pairs with year > `MPA_Start`, ",
             "matching the harmonized BA convention. `missing_years` lists absent internal ",
             "calendar years separated by semicolons. The index also records the input ",
             "dataset, full taxon name, final `N`, analysis type, selected model, ",
             "missing-value counts, common scale and export specifications."),
      "",
      "## Outputs and reproduction", "",
      "- `plots/fig_s_pair_history_atlas.pdf`: combined multipage vector PDF.",
      "- `plots/pair_history_atlas/page_001.png` onward: 170 x 200 mm pages at 300 dpi.",
      "- `plots/pair_history_atlas/page_001.pdf` onward: matching single-page vector PDFs.",
      paste0("- `plots/pair_history_atlas/qa/contact_sheet_01.png` onward: six-page QA contact ",
             "sheets, ordered left to right, then top to bottom. Sheet 01 contains pages ",
             "001-006, sheet 02 pages 007-012, and so on."),
      "- `plots/pair_history_atlas/qa/page_001_{deutan,protan,desaturate}.png`: accessibility checks.",
      "",
      "Run `Rscript code/R/29_pair_history_atlas.R` or `source(\"code/R/29_pair_history_atlas.R\")`.",
      paste0("The script resolves the repository root from its own path, runs in local scope, ",
             "and does not attach packages or load pipeline objects. Required packages are ",
             "ggplot2, ggrepel and ragg. magick and colorspace add reproducible QA assets. ",
             "Figures contain no title or subtitle."),
      "",
      "Input MD5 checksums for this render:", "",
      "| Input | MD5 |", "| --- | --- |",
      sprintf("| `%s` | `%s` |", input_paths, input_md5)
    )
    if (!identical(input_md5, unname(tools::md5sum(absolute_inputs)))) {
      stop("Atlas input changed during rendering; rerun against a consistent final set")
    }
    writeLines(legend, file.path(root, "docs", "pair_history_atlas_legend.md"))
    message(sprintf("Atlas complete: %d series, %d MPAs, %d annual pairs, %d pages",
                    nrow(manifest), n_mpa, sum(manifest$n_pairs), length(pages)))
    invisible(list(manifest = manifest, pdf = combined_pdf, pages = unique(manifest$png)))
  }
  run_atlas()
})
