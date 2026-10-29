# Optional raw-data diagnostic. Run from the analysis repository root.
run_fish_siting <- function(root = here::here()) {
  suppressPackageStartupMessages({library(data.table); library(dplyr); library(tidyr)})
  source(file.path(root, "code/R/primary_uncertainty_utils.R"), local = TRUE)
  source(file.path(root, "code/R/review_concern_utils.R"), local = TRUE)
  source(file.path(root, "code/R/fish_siting_utils.R"), local = TRUE)
  out <- function(x, name) write.csv(x, file.path(root, "tables", paste0("table_s_fish_siting_", name, ".csv")), row.names = FALSE)
  paths <- file.path(root, c("data/PISCO/MLPA_kelpforest_fish_2024.csv",
    "data/PISCO/spp_attribute_table_downloaded_9-13-22_SHSPUL.csv",
    "data/PISCO/master_site_table_Emilyedit.csv", "data/harmonized/reference_sites.csv",
    "data/harmonized/reference_primary_keys.csv", "data/Site_List_All.csv"))
  if (any(!file.exists(paths))) {
    out(data.frame(status = "skipped_missing_raw_inputs", path = paths[!file.exists(paths)]), "status")
    return(invisible(NULL))
  }
  protected <- c(list.files(file.path(root, "data/harmonized"), full.names = TRUE),
    file.path(root, "data/sumstats_final.csv"), list.files(file.path(root, "plots"), pattern = "^fig_0", full.names = TRUE))
  protected <- protected[!dir.exists(protected)]
  hashes <- tools::md5sum(protected)
  inputs <- c(paths, file.path(root, c("code/R/47_fish_siting.R", "code/R/fish_siting_utils.R",
    "code/R/primary_uncertainty_utils.R", "code/R/review_concern_utils.R", "docs/FISH_SITING_PLAN.md")))
  input_hash <- tools::md5sum(inputs)
  codes <- c("SPUL", "PCLA", "OCAL", "CPUN", "HSEM", "HRUB", "BFRE", "OPIC", "HROS")
  attr <- as.data.frame(fread(paths[2]))
  panel <- attr[match(codes, attr$pisco_classcode), c("pisco_classcode", "common_name", "Targeted",
    "yoy_size", "juv_size_cutoff", "thermal_composite_class_from_Freedman_et_al_2021")]
  names(panel)[c(1, 6)] <- c("species", "thermal")
  stopifnot(identical(panel$species, codes), all(panel$Targeted == c(rep("Targeted", 2), rep("Nontargeted", 7))))
  panel$cutoff <- ifelse(is.finite(panel$yoy_size), panel$yoy_size, panel$juv_size_cutoff)
  stopifnot(all(is.finite(panel$cutoff)))
  out(panel, "panel")
  keys <- as.data.frame(fread(paths[5])) %>% filter(Source == "PISCO", Taxa == "S. pulcher", Resp == "Den")
  stopifnot(!anyDuplicated(keys$MPA))
  starts <- as.data.frame(fread(paths[6])) %>% distinct(CA_MPA_Name_Short, MPA_Start)
  keys <- left_join(keys, starts, by = c("MPA" = "CA_MPA_Name_Short"))
  stopifnot(!anyNA(keys$MPA_Start))
  registry <- as.data.frame(fread(paths[4])) %>% filter(Source == "PISCO", Taxa == "S. pulcher", Resp == "Den")
  map <- as.data.frame(fread(paths[3])) %>% filter(BaselineRegion == "SOUTH", CA_MPA_Name_Short %in% keys$MPA) %>%
    transmute(MPA = CA_MPA_Name_Short, site, site_new, status = site_status) %>% distinct()
  stopifnot(!anyDuplicated(map[c("MPA", "site")]))
  refs <- registry %>% filter(included, !excluded) %>% distinct(MPA, site = site_id)
  map <- bind_rows(filter(map, status == "mpa"), semi_join(filter(map, status == "reference"), refs, by = c("MPA", "site")))
  # The same exclusions used upstream also apply to inside records.
  excluded <- unique(registry$site_id[registry$excluded])
  map <- filter(map, !site %in% excluded)
  fish <- fread(paths[1])
  fish <- fish[campus %in% c("UCSB", "VRG") & method %in% c("SBTL_FISH_PISCO", "SBTL_FISH_VRG") &
    LOOKEDFOR == "yes" & site %in% map$site & level %in% c("BOT", "MID")]
  stopifnot(!anyNA(fish[, .(site, year, month, day, zone, level, transect)]))
  # Survey cells are reconstructed from all observed fishes, not the focal panel.
  cell_keys <- c("campus", "site", "year", "month", "day", "zone", "level", "transect")
  fish[, cell := .GRP, by = cell_keys]
  cells <- unique(fish[, c("cell", cell_keys), with = FALSE])
  focal <- fish[classcode %in% codes]
  focal[, cutoff := panel$cutoff[match(classcode, panel$species)]]
  bad <- unique(focal[!is.finite(count) | count < 0 | (count > 0 & (!is.finite(fish_tl) | fish_tl <= 0)), cell])
  out(as.data.frame(cells[cell %in% bad]), "invalid_cells")
  cells <- cells[!cell %in% bad]; focal <- focal[!cell %in% bad]
  summed <- focal[, .(all_sizes = sum(count), size_filtered = sum(count[fish_tl > cutoff])), by = .(cell, species = classcode)]
  grid <- as.data.table(merge(as.data.frame(cells), data.frame(species = codes), by = NULL))
  grid <- merge(grid, summed, by = c("cell", "species"), all.x = TRUE)
  grid[is.na(all_sizes), `:=`(all_sizes = 0, size_filtered = 0)]
  # Matched vertical strata in each visit and zone; no canopy-dependent sampling.
  visit <- grid[, .(all_sizes = mean(all_sizes), size_filtered = mean(size_filtered), n_transects = .N),
    by = .(campus, site, year, month, day, zone, level, species)]
  complete <- visit[, .(n_levels = uniqueN(level)), by = .(campus, site, year, month, day, zone)]
  visit <- merge(visit, complete[n_levels == 2], by = c("campus", "site", "year", "month", "day", "zone"))
  annual_zone <- visit[, .(all_sizes = mean(all_sizes), size_filtered = mean(size_filtered)),
    by = .(site, year, zone, level, species)]
  zones <- unique(annual_zone[, .(site, year, zone)])
  out(zones, "zones")
  site_outputs <- list()
  for (frame in c("bottom_mid", "bottom_only")) {
    z <- if (frame == "bottom_only") annual_zone[level == "BOT"] else annual_zone
    z <- z[, .(all_sizes = sum(all_sizes), size_filtered = sum(size_filtered)), by = .(site, year, zone, species)]
    z <- z[, .(all_sizes = mean(all_sizes) / 60, size_filtered = mean(size_filtered) / 60), by = .(site, year, species)]
    z <- melt(z, id.vars = c("site", "year", "species"), variable.name = "size", value.name = "density")
    z[, frame := frame]; site_outputs[[frame]] <- z
  }
  site <- as.data.frame(rbindlist(site_outputs))
  out(site, "site_annual")
  coverage <- list(); rosters <- list(); omitted <- list(); pair_rows <- list()
  for (i in seq_len(nrow(keys))) {
    m <- keys$MPA[i]; years <- as.integer(strsplit(keys$observed_years[i], ",", fixed = TRUE)[[1]])
    sy <- unique(site[site$year %in% years, c("site", "year")])
    roster <- map %>% filter(MPA == m, site %in% sy$site)
    rosters[[m]] <- roster
    common <- fish_siting_complete_years(sy, roster$site)
    if (!all(c("mpa", "reference") %in% roster$status)) common <- integer()
    common <- setdiff(common, keys$MPA_Start[i])
    pre <- sum(common < keys$MPA_Start[i]); post <- sum(common > keys$MPA_Start[i])
    after <- common[common > keys$MPA_Start[i]]
    coverage[[m]] <- data.frame(MPA = m, start = keys$MPA_Start[i], n_sites = nrow(roster),
      n_inside = sum(roster$status == "mpa"), n_reference = sum(roster$status == "reference"),
      n_before = pre, n_after = post, before_after_eligible = pre >= 3 & post >= 5,
      tracking_eligible = pre >= 6, post_trend_eligible = length(after) >= 8 && diff(range(after)) >= 8,
      years = paste(sort(common), collapse = ";"), anchor_years = paste(years, collapse = ";"))
    for (y in setdiff(years, common)) omitted[[length(omitted) + 1L]] <- data.frame(MPA = m, year = y,
      reason = if (y == keys$MPA_Start[i]) "designation year" else "incomplete fixed site/vertical roster",
      missing_sites = paste(setdiff(roster$site, sy$site[sy$year == y]), collapse = ";"))
    if (!length(common)) next
    z <- inner_join(filter(site, year %in% common), roster, by = "site") %>%
      group_by(MPA, frame, size, species, year, status) %>% summarise(density = mean(density), .groups = "drop") %>%
      pivot_wider(names_from = status, values_from = density)
    stopifnot(all(is.finite(z$mpa)), all(is.finite(z$reference)))
    z$start <- keys$MPA_Start[i]; pair_rows[[m]] <- z
  }
  coverage <- bind_rows(coverage); roster <- bind_rows(rosters)
  ids <- roster %>% mutate(identity = ifelse(is.na(site_new) | !nzchar(site_new), site, site_new)) %>%
    group_by(MPA) %>% summarise(ids = paste(sort(unique(c(site, identity))), collapse = ";"), .groups = "drop")
  ids$reference_group <- review_reference_components(ids$MPA, ids$ids)
  out(coverage, "coverage"); out(roster, "roster"); out(bind_rows(omitted), "omissions"); out(ids, "groups")
  pairs <- bind_rows(pair_rows); all_series <- list(); tracking <- list()
  species_coverage <- pairs %>% filter(frame == "bottom_mid", size == "all_sizes") %>%
    group_by(MPA, species) %>% summarise(n_before = sum(year < start), n_after = sum(year > start),
      positive_inside_before = sum(mpa > 0 & year < start), positive_reference_before = sum(reference > 0 & year < start),
      positive_inside_after = sum(mpa > 0 & year > start), positive_reference_after = sum(reference > 0 & year > start),
      both_zero_years = sum(mpa == 0 & reference == 0), .groups = "drop")
  out(species_coverage, "species_coverage")
  for (offset in c(.1, .5, 1)) {
    z <- pairs %>% mutate(log_inside = log(mpa + offset / 60), log_reference = log(reference + offset / 60),
      contrast = log_inside - log_reference, offset_count = offset)
    z <- left_join(z, panel, by = "species")
    track <- z %>% filter(year < start) %>% group_by(MPA, frame, size, species, offset_count) %>%
      summarise(n_before = n(), correlation = fish_siting_tracking(year, log_inside, log_reference), .groups = "drop")
    tracking[[length(tracking) + 1L]] <- track
    all_series[[length(all_series) + 1L]] <- z %>% transmute(MPA, frame, size, year, start, offset_count,
      series = species, contrast)
    groups <- bind_rows(z %>% filter(Targeted == "Targeted") %>% mutate(series = "Targeted"),
      z %>% filter(Targeted == "Nontargeted", thermal == "Warm") %>% mutate(series = "Nontargeted_warm"),
      z %>% filter(Targeted == "Nontargeted") %>% mutate(series = "Nontargeted_all")) %>%
      group_by(MPA, frame, size, year, start, offset_count, series) %>% summarise(contrast = mean(contrast), .groups = "drop")
    all_series[[length(all_series) + 1L]] <- groups
    w <- groups %>% pivot_wider(names_from = series, values_from = contrast)
    all_series[[length(all_series) + 1L]] <- bind_rows(
      w %>% transmute(MPA, frame, size, year, start, offset_count, series = "Difference_warm", contrast = Targeted - Nontargeted_warm),
      w %>% transmute(MPA, frame, size, year, start, offset_count, series = "Difference_all", contrast = Targeted - Nontargeted_all))
  }
  series <- bind_rows(all_series)
  out(pairs, "paired_density"); out(series, "annual"); out(bind_rows(tracking), "tracking")
  fits <- list(); predictions <- list()
  split_keys <- c("MPA", "frame", "size", "offset_count", "series")
  for (d in split(series, interaction(series[split_keys], drop = TRUE))) {
    for (endpoint in c("before_after", "post_trend")) for (lag in 1:3) {
      f <- fish_siting_fit(d, endpoint, unique(d$start), lag)
      if (is.null(f)) next
      header <- d[1, split_keys, drop = FALSE]
      fits[[length(fits) + 1L]] <- cbind(header, f$result)
      if (lag == 2L) predictions[[length(predictions) + 1L]] <- cbind(header, endpoint, f$predictions)
    }
  }
  fits <- left_join(bind_rows(fits), ids[c("MPA", "reference_group")], by = "MPA")
  out(fits, "fits"); out(bind_rows(predictions), "predictions")
  pooled <- list(); loo <- list()
  split_keys <- c("frame", "size", "offset_count", "series", "endpoint")
  for (d in split(filter(fits, hac_lag == 2), interaction(filter(fits, hac_lag == 2)[split_keys], drop = TRUE))) {
    summary <- review_group_summary(d$estimate, d$reference_group)
    pooled[[length(pooled) + 1L]] <- cbind(d[1, split_keys, drop = FALSE], summary)
    for (g in unique(d$reference_group)) loo[[length(loo) + 1L]] <- cbind(d[1, split_keys, drop = FALSE],
      data.frame(omitted_group = g, estimate = mean(d$estimate[d$reference_group != g]), n_remaining = sum(d$reference_group != g)))
  }
  out(bind_rows(pooled), "pooled"); out(bind_rows(loo), "leave_group_out")
  stopifnot(identical(unname(hashes), unname(tools::md5sum(protected))),
    identical(unname(input_hash), unname(tools::md5sum(inputs))))
  out(data.frame(path = inputs, md5 = unname(input_hash)), "inputs")
  out(data.frame(path = protected, before_md5 = unname(hashes), after_md5 = unname(tools::md5sum(protected))), "protected")
  out(data.frame(status = "complete", pairs_screened = nrow(keys), invalid_cells_removed = length(bad),
    before_after_pairs = sum(coverage$before_after_eligible), post_trend_pairs = sum(coverage$post_trend_eligible)), "status")
  capture.output(sessionInfo(), file = file.path(root, "tables/fish_siting_sessionInfo.txt"))
  print(coverage)
  print(bind_rows(pooled) %>% filter(frame == "bottom_mid", size == "all_sizes", offset_count == .5,
    series %in% c("Targeted", "Nontargeted_warm", "Difference_warm")) %>% select(series, endpoint, n_mpa, n_reference_groups, estimate, ci_low, ci_high))
  invisible(list(coverage = coverage, fits = fits, pooled = bind_rows(pooled)))
}

if (sys.nframe() == 0L) run_fish_siting()
