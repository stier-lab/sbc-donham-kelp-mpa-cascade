# Coverage-based recovery; module 47 and manuscript analyses remain unchanged.
run_fish_coverage_recovery <- function(root = here::here()) {
  suppressPackageStartupMessages({library(data.table); library(dplyr); library(tidyr)})
  for (name in c("primary_uncertainty_utils.R", "review_concern_utils.R", "fish_siting_utils.R", "fish_recovery_utils.R"))
    source(file.path(root, "code/R", name), local = TRUE)
  out <- function(d, name) write.csv(d, file.path(root, "tables", paste0("table_s_fish_recovery_", name, ".csv")), row.names = FALSE)
  read_old <- function(name) read.csv(file.path(root, "tables", paste0("table_s_fish_siting_", name, ".csv")))
  raw_paths <- file.path(root, c("data/PISCO/MLPA_kelpforest_fish_2024.csv", "data/PISCO/master_site_table_Emilyedit.csv"))
  if (any(!file.exists(raw_paths)) || !file.exists(file.path(root, "tables/table_s_fish_siting_panel.csv"))) {
    out(data.frame(status = "skipped_missing_inputs"), "status"); return(invisible(NULL))
  }
  protected <- c(list.files(file.path(root, "tables"), pattern = "^table_s_fish_siting_", full.names = TRUE),
    list.files(file.path(root, "data/harmonized"), full.names = TRUE), file.path(root, "data/sumstats_final.csv"))
  protected <- protected[!dir.exists(protected)]; hashes <- tools::md5sum(protected)
  inputs <- c(raw_paths, protected, file.path(root, c("code/R/48_fish_coverage_recovery.R",
    "code/R/fish_recovery_utils.R", "code/R/fish_siting_utils.R", "code/R/primary_uncertainty_utils.R",
    "code/R/review_concern_utils.R", "docs/FISH_COVERAGE_RECOVERY_PLAN.md")))
  input_hash <- tools::md5sum(inputs)
  panel <- read_old("panel"); original <- read_old("coverage"); old_roster <- read_old("roster")
  keys <- read.csv(file.path(root, "data/harmonized/reference_primary_keys.csv")) %>%
    filter(Source == "PISCO", Taxa == "S. pulcher", Resp == "Den")
  registry <- read.csv(file.path(root, "data/harmonized/reference_sites.csv")) %>%
    filter(Source == "PISCO", Taxa == "S. pulcher", Resp == "Den")
  map <- as.data.frame(fread(raw_paths[2])) %>% filter(BaselineRegion == "SOUTH", CA_MPA_Name_Short %in% keys$MPA) %>%
    transmute(MPA = CA_MPA_Name_Short, site, site_new, status = site_status) %>% distinct()
  refs <- registry %>% filter(included, !excluded) %>% distinct(MPA, site = site_id)
  map <- bind_rows(filter(map, status == "mpa"), semi_join(filter(map, status == "reference"), refs, by = c("MPA", "site"))) %>%
    filter(!site %in% registry$site_id[registry$excluded])
  stopifnot(!anyDuplicated(map[c("MPA", "site")]))
  identities <- map %>% group_by(MPA, site_new) %>% summarise(n_raw_sites = n_distinct(site),
    statuses = paste(sort(unique(status)), collapse = ";"), sites = paste(sort(unique(site)), collapse = ";"), .groups = "drop")
  out(identities, "identity_audit")
  raw <- fread(raw_paths[1])[site %in% map$site & year <= 2023]
  raw_sy <- unique(raw[, .(site, year)])
  standard <- raw[campus %in% c("UCSB", "VRG") & method %in% c("SBTL_FISH_PISCO", "SBTL_FISH_VRG") & LOOKEDFOR == "yes"]
  std_sy <- unique(standard[, .(site, year)])
  bottom <- standard[level == "BOT"]
  cell_keys <- c("campus", "site", "year", "month", "day", "zone", "level", "transect")
  stopifnot(!anyNA(bottom[, ..cell_keys]))
  bottom[, cell := .GRP, by = cell_keys]
  cells <- unique(bottom[, c("cell", cell_keys), with = FALSE])
  focal <- bottom[classcode %in% panel$species]
  focal[, cutoff := panel$cutoff[match(classcode, panel$species)]]
  bad <- unique(focal[!is.finite(count) | count < 0 | (count > 0 & (!is.finite(fish_tl) | fish_tl <= 0)), cell])
  cells <- cells[!cell %in% bad]; focal <- focal[!cell %in% bad]
  summed <- focal[, .(all_sizes = sum(count), size_filtered = sum(count[fish_tl > cutoff])), by = .(cell, species = classcode)]
  grid <- as.data.table(merge(as.data.frame(cells), data.frame(species = panel$species), by = NULL))
  grid <- merge(grid, summed, by = c("cell", "species"), all.x = TRUE)
  grid[is.na(all_sizes), `:=`(all_sizes = 0, size_filtered = 0)]
  visit <- grid[, .(all_sizes = mean(all_sizes), size_filtered = mean(size_filtered)), by = .(site, year, month, day, zone, species)]
  zone <- visit[, .(all_sizes = mean(all_sizes), size_filtered = mean(size_filtered)), by = .(site, year, zone, species)]
  bottom_site <- zone[, .(all_sizes = mean(all_sizes) / 60, size_filtered = mean(size_filtered) / 60), by = .(site, year, species)]
  bottom_site <- melt(bottom_site, id.vars = c("site", "year", "species"), variable.name = "size", value.name = "density")
  bottom_site[, frame := "bottom_recovered"]
  site <- bind_rows(filter(read_old("site_annual"), frame == "bottom_mid"), as.data.frame(bottom_site))
  out(site, "site_annual"); out(unique(zone[, .(site, year, zone)]), "bottom_zones")
  footprint <- site %>% distinct(frame, site, year)
  pre_quality <- standard[level %in% c("BOT", "MID"), .(levels = uniqueN(level)), by = .(site, year, month, day, zone)][levels == 2]
  has_sy <- function(d, s, y) any(d$site == s & d$year == y)
  missing <- list(); site_history <- list()
  for (i in seq_len(nrow(original))) {
    m <- original$MPA[i]; years <- as.integer(strsplit(original$anchor_years[i], ";", fixed = TRUE)[[1]])
    r <- old_roster %>% filter(MPA == m)
    for (s in r$site) {
      yy <- footprint$year[footprint$frame == "bottom_mid" & footprint$site == s]
      site_history[[length(site_history) + 1L]] <- data.frame(MPA = m, site = s, first = min(yy), last = max(yy),
        n_primary_years = sum(years %in% yy), n_anchor_years = length(years),
        fraction = mean(years %in% yy), years = paste(sort(yy), collapse = ";"))
      for (y in setdiff(years, yy)) {
        reason <- if (!has_sy(raw_sy, s, y)) "no recorded fish survey" else if (!has_sy(std_sy, s, y))
          "no standard-method survey" else if (!has_sy(pre_quality, s, y)) "no matched BOT+MID visit/zone" else "quality-filter loss"
        missing[[length(missing) + 1L]] <- data.frame(MPA = m, site = s, year = y, reason,
          designation_year = y == original$start[i])
      }
    }
  }
  out(bind_rows(missing), "missing_site_years"); out(bind_rows(site_history), "site_history")
  missing_summary <- bind_rows(missing) %>% filter(!designation_year) %>% count(reason, name = "missing_site_year_links")
  out(missing_summary, "missing_summary")
  cover <- list(); rosters <- list(); selection <- list(); annual <- list(); dens <- list()
  for (scope in c("primary_years", "raw_years")) for (frame in c("bottom_mid", "bottom_recovered"))
    for (rule in c("full", "core70", "core80", "core90")) for (endpoint in c("post_trend", "before_after")) {
      scenario <- paste(scope, frame, rule, endpoint, sep = "|")
      fp <- filter(footprint, .data$frame == .env$frame)
      for (i in seq_len(nrow(original))) {
        m <- original$MPA[i]; start <- original$start[i]; r <- filter(map, MPA == m)
        years <- if (scope == "primary_years") as.integer(strsplit(original$anchor_years[i], ";", fixed = TRUE)[[1]]) else
          sort(unique(std_sy$year[std_sy$site %in% r$site]))
        if (rule == "full" && scope == "primary_years") r <- filter(old_roster, MPA == m)
        if (rule == "full") {
          r <- filter(r, site %in% fp$site[fp$year %in% years])
          scores <- data.frame(site = r$site, min_coverage = NA_real_, selected = TRUE)
        } else {
          scores <- fish_recovery_select(fp, r, years, start, endpoint, as.numeric(sub("core", "", rule)) / 100)
          r <- filter(r, site %in% scores$site[scores$selected])
        }
        selection[[length(selection) + 1L]] <- cbind(scenario, MPA = m, scores)
        common <- intersect(years, fish_siting_complete_years(fp, r$site))
        if (!all(c("mpa", "reference") %in% r$status)) common <- integer()
        common <- setdiff(common, start)
        if (endpoint == "post_trend") common <- common[common > start]
        eligible <- fish_recovery_eligible(common, start, endpoint)
        cover[[length(cover) + 1L]] <- data.frame(scenario, scope, frame, rule, endpoint, MPA = m, start,
          n_sites = nrow(r), n_inside = sum(r$status == "mpa"), n_reference = sum(r$status == "reference"),
          n_before = sum(common < start), n_after = sum(common > start), eligible,
          years = paste(sort(common), collapse = ";"), anchor_years = paste(years, collapse = ";"))
        if (nrow(r)) rosters[[length(rosters) + 1L]] <- cbind(scenario, r)
        if (!eligible) next
        z <- site %>% filter(.data$frame == .env$frame, year %in% common) %>% inner_join(r, by = "site") %>%
          group_by(MPA, frame, size, species, year, status) %>% summarise(density = mean(density), .groups = "drop") %>%
          pivot_wider(names_from = status, values_from = density)
        z$scenario <- scenario; z$start <- start; z$endpoint <- endpoint
        dens[[length(dens) + 1L]] <- z
        for (offset in c(.1, .5, 1)) {
          v <- z %>% mutate(contrast = log(mpa + offset / 60) - log(reference + offset / 60)) %>% left_join(panel, by = "species")
          g <- bind_rows(filter(v, Targeted == "Targeted") %>% mutate(series = "Targeted"),
            filter(v, Targeted == "Nontargeted", thermal == "Warm") %>% mutate(series = "Nontargeted_warm"),
            filter(v, Targeted == "Nontargeted") %>% mutate(series = "Nontargeted_all")) %>%
            group_by(scenario, endpoint, MPA, size, year, start, series) %>% summarise(contrast = mean(contrast), .groups = "drop")
          w <- g %>% pivot_wider(names_from = series, values_from = contrast)
          delta <- bind_rows(w %>% transmute(scenario, endpoint, MPA, size, year, start, series = "Difference_warm", contrast = Targeted - Nontargeted_warm),
            w %>% transmute(scenario, endpoint, MPA, size, year, start, series = "Difference_all", contrast = Targeted - Nontargeted_all))
          g <- bind_rows(g, delta); g$offset_count <- offset
          annual[[length(annual) + 1L]] <- g
        }
      }
    }
  coverage <- bind_rows(cover); roster <- bind_rows(rosters); annual <- bind_rows(annual)
  out(coverage, "coverage"); out(roster, "roster"); out(bind_rows(selection), "selection")
  out(annual, "annual"); out(bind_rows(dens), "density")
  fits <- list()
  split_keys <- c("scenario", "endpoint", "MPA", "size", "offset_count", "series")
  for (d in split(annual, interaction(annual[split_keys], drop = TRUE))) {
    f <- fish_siting_fit(d, unique(d$endpoint), unique(d$start), 2)
    stopifnot(!is.null(f))
    fits[[length(fits) + 1L]] <- cbind(d[1, setdiff(split_keys, "endpoint"), drop = FALSE], f$result)
  }
  fits <- bind_rows(fits); out(fits, "fits")
  groups <- roster %>% mutate(identity = ifelse(is.na(site_new) | !nzchar(site_new), site, site_new)) %>%
    group_by(scenario, MPA) %>% summarise(ids = paste(sort(unique(c(site, identity))), collapse = ";"), .groups = "drop") %>%
    group_by(scenario) %>% mutate(reference_group = review_reference_components(MPA, ids)) %>% ungroup()
  fits <- left_join(fits, groups[c("scenario", "MPA", "reference_group")], by = c("scenario", "MPA"))
  out(groups, "groups")
  pooled <- list()
  split_keys <- c("scenario", "endpoint", "size", "offset_count", "series")
  for (d in split(fits, interaction(fits[split_keys], drop = TRUE)))
    pooled[[length(pooled) + 1L]] <- cbind(d[1, split_keys, drop = FALSE], review_group_summary(d$estimate, d$reference_group))
  pooled <- bind_rows(pooled); out(pooled, "pooled")
  # Compare old/new rosters on identical years at the original trend pairs.
  old <- read_old("annual"); matched <- list()
  scenario <- "primary_years|bottom_mid|core80|post_trend"
  for (m in original$MPA[original$post_trend_eligible]) {
    a <- filter(old, MPA == m, frame == "bottom_mid", size == "all_sizes", offset_count == .5,
      series == "Difference_warm", year > start)
    b <- filter(annual, .data$scenario == .env$scenario, MPA == m, size == "all_sizes", offset_count == .5, series == "Difference_warm")
    yy <- intersect(a$year, b$year)
    for (version in c("strict", "core80")) {
      d <- if (version == "strict") filter(a, year %in% yy) else filter(b, year %in% yy)
      if (!nrow(d)) next
      f <- fish_siting_fit(d, "post_trend", unique(d$start), 2)
      if (!is.null(f)) matched[[length(matched) + 1L]] <- cbind(MPA = m, version, f$result)
    }
  }
  out(bind_rows(matched), "common_year_comparison")
  stopifnot(identical(unname(hashes), unname(tools::md5sum(protected))),
    identical(unname(input_hash), unname(tools::md5sum(inputs))))
  out(data.frame(path = inputs, md5 = unname(input_hash)), "inputs")
  out(data.frame(path = protected, before_md5 = unname(hashes), after_md5 = unname(tools::md5sum(protected))), "protected")
  out(data.frame(status = "complete", scenarios = n_distinct(coverage$scenario)), "status")
  print(missing_summary)
  print(coverage %>% group_by(scope, frame, rule, endpoint) %>% summarise(n_eligible = sum(eligible), .groups = "drop"))
  print(filter(pooled, size == "all_sizes", offset_count == .5, series == "Difference_warm",
    grepl("core80", scenario)) %>% select(scenario, n_mpa, n_reference_groups, estimate, ci_low, ci_high))
  invisible(list(coverage = coverage, pooled = pooled))
}
if (sys.nframe() == 0L) run_fish_coverage_recovery()
