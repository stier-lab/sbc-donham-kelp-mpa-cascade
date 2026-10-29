# Optional diagnostics: coverage selection and changing within-site sampling.
run_fish_retention_zones <- function(root = here::here()) {
  suppressPackageStartupMessages({library(data.table); library(dplyr); library(tidyr)})
  for (name in c("primary_uncertainty_utils.R", "review_concern_utils.R", "fish_siting_utils.R", "fish_recovery_utils.R"))
    source(file.path(root, "code/R", name), local = TRUE)
  out <- function(d, name) write.csv(d, file.path(root, "tables", paste0("table_s_fish_sampling_", name, ".csv")), row.names = FALSE)
  read_tab <- function(prefix, name) read.csv(file.path(root, "tables", paste0("table_s_fish_", prefix, "_", name, ".csv")))
  raw_path <- file.path(root, "data/PISCO/MLPA_kelpforest_fish_2024.csv")
  if (!file.exists(raw_path) || !file.exists(file.path(root, "tables/table_s_fish_recovery_status.csv"))) {
    out(data.frame(status = "skipped_missing_inputs"), "status"); return(invisible(NULL))
  }
  protected <- c(list.files(file.path(root, "tables"), pattern = "^table_s_fish_(siting|recovery)_", full.names = TRUE),
    list.files(file.path(root, "data/harmonized"), full.names = TRUE), file.path(root, "data/sumstats_final.csv"))
  protected <- protected[!dir.exists(protected)]; hashes <- tools::md5sum(protected)
  input_paths <- c(raw_path, protected, file.path(root, c("code/R/49_fish_retention_zones.R",
    "docs/FISH_RETENTION_ZONE_PLAN.md", "code/R/fish_siting_utils.R", "code/R/fish_recovery_utils.R",
    "code/R/primary_uncertainty_utils.R", "code/R/review_concern_utils.R")))
  input_hashes <- tools::md5sum(input_paths)
  panel <- read_tab("siting", "panel"); map <- read_tab("siting", "roster")
  stopifnot(read_tab("recovery", "status")$status == "complete")
  coverage <- read_tab("recovery", "coverage") %>% filter(scope == "primary_years", frame == "bottom_mid", rule == "core80")
  roster <- read_tab("recovery", "roster") %>% filter(scenario %in% coverage$scenario)
  selection <- read_tab("recovery", "selection") %>% filter(scenario %in% coverage$scenario)
  old_site <- read_tab("siting", "site_annual") %>% filter(frame == "bottom_mid", site %in% map$site, year <= 2023)
  raw <- fread(raw_path)[site %in% map$site & year <= 2023 & campus %in% c("UCSB", "VRG") &
    method %in% c("SBTL_FISH_PISCO", "SBTL_FISH_VRG") & LOOKEDFOR == "yes" & level %in% c("BOT", "MID")]
  cell_keys <- c("campus", "site", "year", "month", "day", "zone", "level", "transect")
  stopifnot(!anyNA(raw[, ..cell_keys]))
  raw[, cell := .GRP, by = cell_keys]
  cells <- raw[, .(n_depth = uniqueN(depth[is.finite(depth)]), depth = if (all(!is.finite(depth))) NA_real_ else median(depth, na.rm = TRUE)),
    by = c("cell", cell_keys)]
  stopifnot(all(cells$n_depth <= 1L))
  out(cells[, .(cells = .N, missing_depth = sum(!is.finite(depth)), conflicting_depth = sum(n_depth > 1)), by = .(campus, level)], "metadata")
  focal <- raw[classcode %in% panel$species]
  focal[, cutoff := panel$cutoff[match(classcode, panel$species)]]
  bad <- unique(focal[!is.finite(count) | count < 0 | (count > 0 & (!is.finite(fish_tl) | fish_tl <= 0)), cell])
  cells <- cells[!cell %in% bad]; focal <- focal[!cell %in% bad]
  summed <- focal[, .(all_sizes = sum(count), size_filtered = sum(count[fish_tl > cutoff])), by = .(cell, species = classcode)]
  grid <- as.data.table(merge(as.data.frame(cells), data.frame(species = panel$species), by = NULL))
  grid <- merge(grid, summed, by = c("cell", "species"), all.x = TRUE)
  grid[is.na(all_sizes), `:=`(all_sizes = 0, size_filtered = 0)]
  # Complete vertical coverage is reconstructed after every restriction.
  zone_means <- function(g) {
    v <- g[, .(all_sizes = mean(all_sizes), size_filtered = mean(size_filtered)),
      by = .(campus, site, year, month, day, zone, level, species)]
    both <- v[, .(n_levels = uniqueN(level)), by = .(campus, site, year, month, day, zone)][n_levels == 2]
    v <- merge(v, both, by = c("campus", "site", "year", "month", "day", "zone"))
    a <- v[, .(all_sizes = mean(all_sizes), size_filtered = mean(size_filtered)), by = .(site, year, zone, level, species)]
    a[, .(all_sizes = sum(all_sizes) / 60, size_filtered = sum(size_filtered) / 60), by = .(site, year, zone, species)]
  }
  zone <- zone_means(grid)
  reconstructed <- zone[, .(all_sizes = mean(all_sizes), size_filtered = mean(size_filtered)), by = .(site, year, species)]
  reconstructed <- melt(reconstructed, id.vars = c("site", "year", "species"), variable.name = "size", value.name = "density")
  check <- inner_join(as.data.frame(reconstructed), old_site, by = c("site", "year", "species", "size"), suffix = c("_new", "_old"))
  stopifnot(nrow(check) == nrow(old_site), max(abs(check$density_new - check$density_old)) < 1e-10)
  out(data.frame(n = nrow(check), maximum_difference = max(abs(check$density_new - check$density_old))), "reconstruction")
  bot <- cells[level == "BOT"]
  visits <- cells[, .(n_levels = uniqueN(level)), by = .(campus, site, year, month, day, zone)][n_levels == 2]
  depth_bot <- merge(bot, visits, by = c("campus", "site", "year", "month", "day", "zone"))
  depth_site <- depth_bot[, .(depth = mean(depth, na.rm = TRUE)), by = .(site, year, month, day, zone)][,
    .(depth = mean(depth, na.rm = TRUE)), by = .(site, year, zone)][,
    .(value = mean(depth, na.rm = TRUE)), by = .(site, year)]
  fish_site <- old_site %>% filter(size == "all_sizes") %>% left_join(panel, by = "species") %>%
    filter(Targeted == "Targeted" | thermal == "Warm") %>%
    mutate(metric = ifelse(Targeted == "Targeted", "Targeted", "Nontargeted_warm")) %>%
    group_by(site, year, metric) %>% summarise(value = mean(log(density * 60 + .5)), .groups = "drop")
  metrics <- bind_rows(fish_site, as.data.frame(depth_site) %>% mutate(metric = "Bottom_depth_m"))
  retention <- list(); retention_annual <- list(); retention_audit <- list()
  for (i in seq_len(nrow(coverage))) {
    c <- coverage[i, ]; s <- selection %>% filter(scenario == c$scenario, MPA == c$MPA) %>%
      inner_join(filter(map, MPA == c$MPA), by = c("MPA", "site"))
    anchor <- as.integer(strsplit(c$anchor_years, ";", fixed = TRUE)[[1]])
    for (excluded in s$site[!s$selected]) for (window in c("early3", "all_overlap")) {
      side <- s$status[match(excluded, s$site)]; retained <- s$site[s$selected & s$status == side]
      fp <- old_site %>% filter(year %in% anchor, site %in% c(excluded, retained)) %>% distinct(site, year)
      yy <- sort(intersect(fp$year[fp$site == excluded], fp$year[fp$site %in% retained]))
      if (window == "early3") yy <- head(yy, 3)
      fixed <- retained[vapply(retained, function(k) all(yy %in% fp$year[fp$site == k]), logical(1))]
      usable <- length(yy) > 0 && length(fixed) > 0
      retention_audit[[length(retention_audit) + 1L]] <- data.frame(endpoint = c$endpoint, MPA = c$MPA, site = excluded,
        status = side, window, usable, main_pair_eligible = c$eligible, n_years = length(yy), n_retained = length(fixed),
        years = paste(yy, collapse = ";"), retained = paste(fixed, collapse = ";"))
      if (!usable) next
      for (metric in unique(metrics$metric)) {
        d <- metrics %>% filter(.data$metric == .env$metric, year %in% yy)
        a <- filter(d, site == excluded) %>% select(year, excluded_value = value)
        b <- filter(d, site %in% fixed) %>% group_by(year) %>%
          summarise(retained_value = if (all(is.finite(value)) && n() == length(fixed)) mean(value) else NA_real_, .groups = "drop")
        v <- inner_join(a, b, by = "year") %>% filter(is.finite(excluded_value), is.finite(retained_value))
        if (!nrow(v)) next
        headrow <- data.frame(endpoint = c$endpoint, MPA = c$MPA, site = excluded, status = side, window, metric, main_pair_eligible = c$eligible)
        retention_annual[[length(retention_annual) + 1L]] <- cbind(headrow, v)
        retention[[length(retention) + 1L]] <- cbind(headrow, data.frame(n_years = nrow(v), n_before = sum(v$year < c$start),
          excluded_mean = mean(v$excluded_value), retained_mean = mean(v$retained_value),
          difference = mean(v$excluded_value - v$retained_value)))
      }
    }
  }
  out(bind_rows(retention_audit), "retention_audit"); out(bind_rows(retention), "retention")
  out(bind_rows(retention_annual), "retention_annual")
  zone_rosters <- list(); zone_coverage <- list(); series <- list(); used_sites <- list(); depth_centers <- list(); depth_audit <- list()
  groups <- read_tab("recovery", "groups")
  make_series <- function(d, r, c, version, restriction) {
    z <- d %>% inner_join(r, by = "site") %>% group_by(MPA, year, size, species, status) %>%
      summarise(density = mean(density), .groups = "drop") %>% pivot_wider(names_from = status, values_from = density) %>%
      mutate(contrast = log(mpa + .5 / 60) - log(reference + .5 / 60)) %>% left_join(panel, by = "species")
    g <- bind_rows(filter(z, Targeted == "Targeted") %>% mutate(series = "Targeted"),
      filter(z, Targeted == "Nontargeted", thermal == "Warm") %>% mutate(series = "Nontargeted_warm")) %>%
      group_by(MPA, year, size, series) %>% summarise(contrast = mean(contrast), .groups = "drop")
    w <- g %>% pivot_wider(names_from = series, values_from = contrast)
    g <- bind_rows(g, w %>% transmute(MPA, year, size, series = "Difference_warm", contrast = Targeted - Nontargeted_warm))
    g$endpoint <- c$endpoint; g$start <- c$start; g$version <- version; g$restriction <- restriction
    g$reference_group <- groups$reference_group[match(paste(c$scenario, c$MPA), paste(groups$scenario, groups$MPA))]
    g
  }
  for (i in which(coverage$eligible)) {
    c <- coverage[i, ]; r <- filter(roster, scenario == c$scenario, MPA == c$MPA)
    yy <- as.integer(strsplit(c$years, ";", fixed = TRUE)[[1]])
    z <- zone[site %in% r$site & year %in% yy]
    zr <- z[, .(n_years = uniqueN(year)), by = .(site, zone)][n_years == length(yy)]
    zone_rosters[[length(zone_rosters) + 1L]] <- cbind(endpoint = c$endpoint, MPA = c$MPA, as.data.frame(zr))
    center <- bot[site %in% r$site & year %in% yy, .(center_depth = median(depth, na.rm = TRUE)), by = .(site, zone)]
    depth_centers[[length(depth_centers) + 1L]] <- cbind(endpoint = c$endpoint, MPA = c$MPA, as.data.frame(center))
    for (restriction in c("fixed_zones", "depth_1m", "depth_2m", "depth_3m")) {
      if (restriction == "fixed_zones") zz <- merge(z, zr[, .(site, zone)], by = c("site", "zone")) else {
        width <- as.numeric(sub("m", "", sub("depth_", "", restriction)))
        g <- merge(grid[site %in% r$site & year %in% yy], zr[, .(site, zone)], by = c("site", "zone"))
        g <- merge(g, center, by = c("site", "zone"))
        g <- g[level == "MID" | (is.finite(depth) & is.finite(center_depth) & abs(depth - center_depth) <= width)]
        gb <- unique(g[level == "BOT", .(cell, depth, center_depth)])
        depth_audit[[length(depth_audit) + 1L]] <- data.frame(endpoint = c$endpoint, MPA = c$MPA, restriction,
          width, n_bottom_cells = nrow(gb), missing_depth = sum(!is.finite(gb$depth)),
          max_deviation = if (nrow(gb)) max(abs(gb$depth - gb$center_depth)) else NA_real_)
        zz <- zone_means(g)
      }
      complete <- zz[, .(n_zones = uniqueN(zone)), by = .(site, year)]
      target <- zr[, .(required_zones = .N), by = site]
      complete <- merge(complete, target, by = "site")[n_zones == required_zones]
      years <- intersect(yy, fish_siting_complete_years(complete, r$site))
      if (!all(r$site %in% zr$site)) years <- integer()
      eligible <- fish_recovery_eligible(years, c$start, c$endpoint)
      zone_coverage[[length(zone_coverage) + 1L]] <- data.frame(endpoint = c$endpoint, MPA = c$MPA, restriction,
        n_sites = nrow(r), sites_no_common_zone = sum(!r$site %in% zr$site), original_years = length(yy),
        n_before = sum(years < c$start), n_after = sum(years > c$start), eligible, years = paste(years, collapse = ";"))
      if (!eligible) next
      zz <- zz[year %in% years]
      ss <- zz[, .(all_sizes = mean(all_sizes), size_filtered = mean(size_filtered)), by = .(site, year, species)]
      ss <- melt(ss, id.vars = c("site", "year", "species"), variable.name = "size", value.name = "density")
      used_sites[[length(used_sites) + 1L]] <- cbind(endpoint = c$endpoint, MPA = c$MPA, restriction, as.data.frame(unique(zz[, .(site, year, zone)])))
      series[[length(series) + 1L]] <- make_series(as.data.frame(ss), r, c, "restricted", restriction)
      series[[length(series) + 1L]] <- make_series(filter(old_site, year %in% years), r, c, "original_matched", restriction)
    }
  }
  out(bind_rows(zone_rosters), "zone_roster"); out(bind_rows(depth_centers), "depth_centers")
  out(bind_rows(depth_audit), "depth_audit")
  out(bind_rows(zone_coverage), "coverage"); out(bind_rows(used_sites), "used_zones")
  annual <- bind_rows(series); out(annual, "annual")
  fits <- list(); keys <- c("endpoint", "MPA", "size", "series", "restriction", "version", "reference_group")
  for (d in split(annual, interaction(annual[keys], drop = TRUE))) {
    f <- fish_siting_fit(d, unique(d$endpoint), unique(d$start), 2)
    stopifnot(!is.null(f)); fits[[length(fits) + 1L]] <- cbind(d[1, setdiff(keys, "endpoint"), drop = FALSE], f$result)
  }
  fits <- bind_rows(fits); out(fits, "fits")
  pooled <- list(); keys <- c("endpoint", "size", "series", "restriction", "version")
  for (d in split(fits, interaction(fits[keys], drop = TRUE)))
    pooled[[length(pooled) + 1L]] <- cbind(d[1, keys, drop = FALSE], review_group_summary(d$estimate, d$reference_group))
  out(bind_rows(pooled), "pooled")
  stopifnot(identical(unname(hashes), unname(tools::md5sum(protected))),
    identical(unname(input_hashes), unname(tools::md5sum(input_paths))))
  out(data.frame(path = protected, before_md5 = unname(hashes), after_md5 = unname(tools::md5sum(protected))), "protected")
  out(data.frame(path = input_paths, md5 = unname(input_hashes)), "inputs")
  out(data.frame(status = "complete", invalid_cells = length(bad), canopy_values_present = sum(!is.na(raw$pctcnpy))), "status")
  print(bind_rows(zone_coverage) %>% group_by(endpoint, restriction) %>% summarise(n_eligible = sum(eligible), .groups = "drop"))
  print(bind_rows(pooled) %>% filter(size == "all_sizes", series == "Difference_warm"))
}
if (sys.nframe() == 0L) run_fish_retention_zones()
