# Optional raw-input diagnostic; never updates primary inputs or selects models.
# Rscript code/R/44_sheephead_yoy_check.R [SOURCE_REPO]

yoy_mask <- function(fish, scenario) {
  stopifnot(scenario %in% c("baseline", "remove_lt10", "remove_le10"))
  valid <- fish$classcode == "SPUL" & is.finite(fish$fish_tl) & fish$fish_tl > 0
  small <- if (scenario == "remove_le10") fish$fish_tl <= 10 else fish$fish_tl < 10
  if (scenario == "baseline") return(rep(FALSE, nrow(fish)))
  which_small <- valid & small
  which_small[is.na(which_small)] <- FALSE
  which_small
}

yoy_filter_counts <- function(fish, scenario) {
  fish$count[yoy_mask(fish, scenario)] <- 0
  fish
}

yoy_key <- function(x, cols) do.call(paste, c(x[cols], sep = "|"))

yoy_equal <- function(actual, expected, label, tolerance = 1e-8) {
  if (length(actual) != length(expected) || any(!is.finite(actual)) ||
      any(!is.finite(expected)) ||
      any(abs(actual - expected) > tolerance * pmax(1, abs(expected)))) {
    stop("YOY baseline guard failed: ", label)
  }
  invisible(TRUE)
}

yoy_environment <- function(root, source_root) {
  env <- new.env(parent = globalenv())
  sys.source(file.path(source_root, "code/R/00c_analysis_constants.R"), env)
  # Only definitions/constants: avoid unrelated calibration reads and any runners.
  for (expr in parse(file.path(source_root, "code/R/01_utils.R"))) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        is.symbol(expr[[2]]) &&
        (is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function")) ||
         as.character(expr[[2]]) %in%
           c("EXCLUDED_REFERENCE_SITES", "EXCLUDED_MPAS"))) eval(expr, env)
  }
  sites <- read.csv(file.path(root, "data/PISCO/master_site_table_Emilyedit.csv"))
  env$sites.short <- sites[!duplicated(sites[c("site", "CA_MPA_Name_Short")]), ]
  env
}

yoy_fish_block <- function(source_root) {
  lines <- readLines(file.path(source_root, "code/R/04_pisco_processing.R"), warn = FALSE)
  start <- grep("^fish_att_path <-", lines)
  end <- grep('^cat\\("  Sheephead records:', lines)
  stopifnot(length(start) == 1L, length(end) == 1L, start < end)
  parse(text = lines[start:end])
}

yoy_aggregate_fish <- function(env, block, root, scenario, raw = NULL, prepared = NULL) {
  e <- new.env(parent = env)
  filtered <- FALSE
  for (expr in block) {
    if (!is.null(prepared) && is.call(expr) &&
        identical(expr[[1]], as.name("for"))) {
      stopifnot(identical(expr[[2]], as.name("i")),
        identical(e$Fish.merged[setdiff(names(e$Fish.merged), c("count", "biomass"))],
                  prepared[setdiff(names(prepared), c("count", "biomass"))]))
      # Counts were zeroed before the unchanged merge. Reuse baseline allometry;
      # zero times a finite positive-length weight is exactly zero.
      e$Fish.merged$biomass <- prepared$biomass
      e$Fish.merged$biomass[yoy_mask(e$Fish.merged, scenario)] <- 0
      next
    }
    if (is.call(expr) && identical(expr[[1]], as.name("<-"))) {
      lhs <- as.character(expr[[2]])
      if (length(lhs) == 1L && lhs == "fish_att_path") {
        e$fish_att_path <- file.path(root, "data/PISCO/spp_attribute_table_downloaded_9-13-22_SHSPUL.csv")
        next
      }
      if (length(lhs) == 1L && lhs == "fish_path") {
        e$fish_path <- file.path(root, "data/PISCO/MLPA_kelpforest_fish_2024.csv")
        next
      }
      if (!is.null(raw) && identical(expr[[2]], as.name("Fish")) &&
          is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("tryCatch"))) {
        e$Fish <- raw
        next
      }
    }
    eval(expr, e)
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        identical(expr[[2]], as.name("Fish")) && is.call(expr[[3]]) &&
        identical(expr[[3]][[1]], as.name("subset"))) {
      stopifnot(!filtered)
      e$raw_fish <- e$Fish
      e$Fish <- yoy_filter_counts(e$Fish, scenario)
      filtered <- TRUE
    }
  }
  stopifnot(filtered)
  list(site = as.data.frame(e$Fish.SPUL.All), raw = e$raw_fish,
       roster = as.data.frame(e$Site.yr.fish), prepared = e$Fish.merged)
}

yoy_eligible_sites <- function(root, env, site) {
  swath <- read.csv(file.path(root, "data/PISCO/MLPA_kelpforest_swath_2024.csv"))
  swath <- subset(swath, (campus == "VRG" |
    (campus == "UCSB" & zone != "DEEP" & zone != "SPECIALO")) &
    month >= env$SURVEY_MONTH_START & month <= env$SURVEY_MONTH_END)
  south <- subset(env$sites.short, BaselineRegion == "SOUTH" & !is.na(CA_MPA_Name_Short))$site
  # Swath's complete species roster ensures its target taxa cover these site-years.
  siteyears <- unique(rbind(site[c("site", "year")],
    subset(swath, site %in% south)[c("site", "year")]))
  years <- table(siteyears$site)
  names(years)[years >= env$MINIMUM_YEARS_OF_DATA]
}

yoy_annual_values <- function(site, eligible, env) {
  site <- subset(site, site %in% eligible & !site %in% env$EXCLUDED_REFERENCE_SITES &
    !CA_MPA_Name_Short %in% env$EXCLUDED_MPAS)
  x <- site %>% dplyr::group_by(CA_MPA_Name_Short, year, y, site_status, site_designation) %>%
    dplyr::summarise(biomass = mean(biomass, na.rm = TRUE),
      density = mean(density, na.rm = TRUE), .groups = "drop")
  long <- rbind(data.frame(x[c("CA_MPA_Name_Short", "year", "site_status")],
    resp = "Bio", value = x$biomass / env$PISCO_SWATH_AREA_M2),
    data.frame(x[c("CA_MPA_Name_Short", "year", "site_status")],
    resp = "Den", value = x$density))
  names(long)[names(long) == "site_status"] <- "status"
  list(mpa = x, long = long)
}

yoy_ratios <- function(mpa, env) {
  do.call(rbind, lapply(c("Bio", "Den"), function(resp) {
    col <- if (resp == "Bio") "biomass" else "density"
    d <- subset(mpa, site_status != "")
    if (resp == "Bio") d <- subset(d, site_designation != "N/A")
    d <- d[is.finite(d[[col]]), ]
    d <- env$calculate_proportions(d, col)
    d <- tidyr::spread(d[c("CA_MPA_Name_Short", "year", "y", "site_designation",
                          "site_status", "PropCorr")], site_status, PropCorr)
    d <- env$calculate_log_response_ratio(d)
    data.frame(MPA = d$CA_MPA_Name_Short, Resp = resp, year = d$year, lnDiff = d$lnDiff)
  }))
}

yoy_fixed_fit <- function(row, rebuilt, lnDiff) {
  stopifnot(row$Model %in% c("Mean", "Linear"), length(lnDiff) == nrow(rebuilt$data),
            all(is.finite(lnDiff)))
  d <- rebuilt$data
  d$lnDiff <- lnDiff
  fit <- if (row$Model == "Mean") lm(lnDiff ~ 1, data = d) else lm(lnDiff ~ time.model, data = d)
  g <- rebuilt$gradient
  c(estimate = sum(g * coef(fit)), se = sqrt(drop(t(g) %*% vcov(fit) %*% g)))
}

run_sheephead_yoy_check <- function(root = here::here(),
    source_root = file.path(dirname(root), "Donham-Stier-CA-MPA-Data-2026"), write = TRUE) {
  suppressPackageStartupMessages({library(dplyr); library(tidyr)})
  raw_paths <- file.path(root, "data/PISCO", c("MLPA_kelpforest_fish_2024.csv",
    "MLPA_kelpforest_swath_2024.csv", "master_site_table_Emilyedit.csv",
    "spp_attribute_table_downloaded_9-13-22_SHSPUL.csv"))
  source_paths <- file.path(source_root, "code/R", c("04_pisco_processing.R", "01_utils.R",
                                                               "00c_analysis_constants.R"))
  unavailable <- c(raw_paths, source_paths)[!file.exists(c(raw_paths, source_paths))]
  if (length(unavailable)) {
    status <- data.frame(status = "skipped_raw_unavailable", path = unavailable)
    if (write) write.csv(status, file.path(root, "tables/table_s_sheephead_yoy_status.csv"), row.names = FALSE)
    message("Optional YOY diagnostic skipped: raw inputs/source unavailable")
    return(invisible(list(status = status)))
  }
  source(file.path(root, "code/R/primary_uncertainty_utils.R"), local = TRUE)
  protected <- unique(c(list.files(file.path(root, "data/harmonized"), full.names = TRUE),
    list.files(file.path(root, "data/cache"), full.names = TRUE),
    list.files(file.path(source_root, "data/cache"), full.names = TRUE),
    file.path(root, c("data/sumstats_final.csv", "tables/table_02_meta_analysis.csv"))))
  protected <- protected[!dir.exists(protected)]
  protected_hash <- tools::md5sum(protected)
  inputs <- unique(c(raw_paths, source_paths, file.path(root, c(
    "data/sumstats_final.csv", "data/harmonized/harmonized_raw_responses.csv",
    "data/harmonized/harmonized_response_ratios.csv", "data/harmonized/harmonized_landsat_rr.csv",
    "code/R/44_sheephead_yoy_check.R", "code/R/primary_uncertainty_utils.R",
    "code/R/08_effect_sizes.R", "code/R/02_pBACIPS_function.R", "tables/table_02_meta_analysis.csv"))))
  fingerprints <- data.frame(path = inputs, md5 = unname(tools::md5sum(inputs)))
  env <- yoy_environment(root, source_root)
  block <- yoy_fish_block(source_root)
  attr <- read.csv(raw_paths[4])
  attr <- subset(attr, pisco_classcode == "SPUL")
  stopifnot(nrow(attr) == 1L, attr$yoy_size == 10, attr$juv_size_cutoff == 10,
            attr$yoy_code == "SPULYOY")
  baseline <- yoy_aggregate_fish(env, block, root, "baseline")
  eligible <- yoy_eligible_sites(root, env, baseline$site)
  base_values <- yoy_annual_values(baseline$site, eligible, env)
  retained <- read.csv(file.path(root, "data/harmonized/harmonized_raw_responses.csv"))
  retained <- subset(retained, source == "PISCO" & taxon_name == "Semicossyphus pulcher")
  rawkeys <- c("CA_MPA_Name_Short", "year", "status", "resp")
  stopifnot(!anyDuplicated(retained[rawkeys]), !anyDuplicated(base_values$long[rawkeys]))
  j <- match(yoy_key(retained, rawkeys), yoy_key(base_values$long, rawkeys))
  stopifnot(!anyNA(j))
  yoy_equal(base_values$long$value[j], retained$value, "exact retained annual raw values")
  primary <- read.csv(file.path(root, "data/sumstats_final.csv"))
  annual <- uncertainty_annual(root, primary)
  definitions <- uncertainty_definitions(root)
  ix <- which(primary$Taxa == "S. pulcher" & primary$Source == "PISCO")
  stopifnot(length(ix) > 0, all(primary$Model[ix] %in% c("Mean", "Linear")))
  rebuilt <- lapply(ix, function(i) uncertainty_reconstruct(primary[i, ], annual, definitions))
  base_ratios <- yoy_ratios(base_values$mpa, env)
  rrkeys <- c("MPA", "Resp", "year")
  for (k in seq_along(ix)) {
    d <- rebuilt[[k]]$data
    j <- match(yoy_key(d, rrkeys), yoy_key(base_ratios, rrkeys))
    stopifnot(!anyNA(j))
    yoy_equal(base_ratios$lnDiff[j], d$lnDiff, "retained annual lnRR")
    f <- yoy_fixed_fit(primary[ix[k], ], rebuilt[[k]], base_ratios$lnDiff[j])
    yoy_equal(f, c(primary$Mean[ix[k]], primary$SE[ix[k]]), "primary contrast and SE")
  }
  message("PASS: baseline annual values, lnRR, all PISCO sheephead contrasts and SEs")
  pooled_baseline <- uncertainty_joint(primary, primary$SE, "baseline")
  published <- read.csv(file.path(root, "tables/table_02_meta_analysis.csv"))
  published$Resp <- ifelse(published$Response == "Biomass", "Bio", "Den")
  match_pooled <- merge(pooled_baseline, published, by = c("Taxa", "Resp"))
  stopifnot(nrow(match_pooled) == 9L)
  yoy_equal(match_pooled$estimate, match_pooled$Estimate, "joint baseline estimate")
  yoy_equal(match_pooled$se, match_pooled$SE, "joint baseline SE")
  yoy_equal(match_pooled$q, match_pooled$pval_fdr, "joint baseline nine-test BH")
  contrasts <- list(); contributions <- list(); pooled <- list(pooled_baseline)
  for (scenario in c("baseline", "remove_lt10", "remove_le10")) {
    agg <- if (scenario == "baseline") baseline else
      yoy_aggregate_fish(env, block, root, scenario, baseline$raw, baseline$prepared)
    stopifnot(identical(agg$roster, baseline$roster),
      identical(agg$site[setdiff(names(agg$site), c("count", "density", "biomass"))],
                baseline$site[setdiff(names(baseline$site), c("count", "density", "biomass"))]))
    values <- yoy_annual_values(agg$site, eligible, env)
    j <- match(yoy_key(retained, rawkeys), yoy_key(values$long, rawkeys))
    stopifnot(!anyNA(j), identical(is.na(values$long$value[j]), is.na(retained$value)))
    filtered <- values$long$value[j]
    stopifnot(all(filtered >= 0), all(filtered <= retained$value + 1e-8))
    contribution <- retained[rawkeys]
    contribution$scenario <- scenario
    contribution$baseline <- retained$value
    contribution$filtered <- filtered
    contribution$removed <- retained$value - filtered
    contribution$fraction_removed <- ifelse(retained$value > 0,
                                           contribution$removed / retained$value, NA_real_)
    contributions[[scenario]] <- contribution
    rr <- yoy_ratios(values$mpa, env)
    updated <- primary
    contrasts[[scenario]] <- do.call(rbind, lapply(seq_along(ix), function(k) {
      i <- ix[k]; row <- primary[i, ]; d <- rebuilt[[k]]$data
      j <- match(yoy_key(d, rrkeys), yoy_key(rr, rrkeys))
      stopifnot(!anyNA(j))
      fitted <- yoy_fixed_fit(row, rebuilt[[k]], rr$lnDiff[j])
      updated$Mean[i] <<- fitted["estimate"]
      updated$SE[i] <<- fitted["se"]
      data.frame(row[c("MPA", "Taxa", "Source", "Resp", "Model", "AnalysisType", "N")],
        scenario, baseline_estimate = row$Mean, baseline_se = row$SE,
        estimate = unname(fitted["estimate"]), se = unname(fitted["se"]),
        delta = unname(fitted["estimate"]) - row$Mean,
        years = paste(d$year, collapse = ","), input_time = paste(d$input_time, collapse = ","))
    }))
    if (scenario != "baseline") pooled[[scenario]] <- uncertainty_joint(updated, updated$SE, scenario)
    stopifnot(identical(updated[-ix, ], primary[-ix, ]))
  }
  # Physical raw-record shares: semi-join to retained sampling cells, then unique row IDs.
  fish <- baseline$raw
  fish$raw_row <- seq_len(nrow(fish))
  raw <- subset(fish, classcode == "SPUL")
  raw <- merge(raw, env$sites.short[c("site", "CA_MPA_Name_Short", "site_status", "BaselineRegion")], by = "site")
  raw <- subset(raw, site %in% eligible & !site %in% env$EXCLUDED_REFERENCE_SITES &
    !CA_MPA_Name_Short %in% env$EXCLUDED_MPAS & BaselineRegion == "SOUTH" & level != "CAN")
  cells <- unique(retained[c("CA_MPA_Name_Short", "year", "status")])
  names(cells)[3] <- "site_status"
  raw <- merge(raw, cells, by = c("CA_MPA_Name_Short", "year", "site_status"))
  raw <- raw[!duplicated(raw$raw_row), ]
  raw_summary <- do.call(rbind, lapply(list(all_campus_SPUL = subset(fish, classcode == "SPUL"),
                                           retained_unique_SPUL = raw), function(d) {
    valid <- is.finite(d$fish_tl) & d$fish_tl > 0
    mass <- d$count * env$bio_sheephead(d$fish_tl)
    do.call(rbind, lapply(c("remove_lt10", "remove_le10"), function(scenario) {
      small <- yoy_mask(d, scenario)
      data.frame(scenario, rows = nrow(d), count = sum(d$count, na.rm = TRUE),
        removed_rows = sum(small), removed_count = sum(d$count[small], na.rm = TRUE),
        count_fraction = sum(d$count[small], na.rm = TRUE) / sum(d$count, na.rm = TRUE),
        biomass_fraction = sum(mass[small], na.rm = TRUE) / sum(mass, na.rm = TRUE),
        at10_count = sum(d$count[valid & d$fish_tl == 10], na.rm = TRUE),
        positive_missing_size_rows = sum(d$count > 0 & !valid, na.rm = TRUE),
        positive_missing_size_count = sum(d$count[d$count > 0 & !valid], na.rm = TRUE),
        missing_count_rows = sum(is.na(d$count)), zero_count_rows = sum(d$count == 0, na.rm = TRUE))
    }))
  }))
  raw_summary$population <- rep(c("all_campus_SPUL", "retained_unique_SPUL"), each = 2)
  contrasts <- do.call(rbind, contrasts)
  contributions <- do.call(rbind, contributions)
  a <- subset(contributions, scenario == "remove_lt10")
  b <- subset(contributions, scenario == "remove_le10")
  stopifnot(identical(yoy_key(a, rawkeys), yoy_key(b, rawkeys)),
            all(b$filtered <= a$filtered + 1e-8))
  pooled <- do.call(rbind, pooled)
  pooled <- merge(pooled, pooled_baseline[c("Taxa", "Resp", "estimate", "q")],
                  by = c("Taxa", "Resp"), suffixes = c("", "_baseline"))
  pooled$delta <- pooled$estimate - pooled$estimate_baseline
  pooled$bh_detection_changed <- (pooled$q < .05) != (pooled$q_baseline < .05)
  protection <- data.frame(path = protected, before_md5 = unname(protected_hash),
                           after_md5 = unname(tools::md5sum(protected)))
  stopifnot(identical(protection$before_md5, protection$after_md5),
    identical(fingerprints$md5, unname(tools::md5sum(inputs))))
  status <- data.frame(status = "complete_baseline_verified", primary_contrasts = length(ix),
    raw_annual_cells = nrow(retained), raw_SPUL_rows = sum(fish$classcode == "SPUL"),
    raw_SPULYOY_rows = sum(fish$classcode == "SPULYOY"), protected_files = length(protected))
  outputs <- list(status = status, annual = contributions, raw_shares = raw_summary,
    contrasts = contrasts, pooled = pooled, inputs = fingerprints, protected = protection)
  if (write) for (name in names(outputs)) write.csv(outputs[[name]],
    file.path(root, "tables", paste0("table_s_sheephead_yoy_", name, ".csv")), row.names = FALSE)
  print(status); print(raw_summary)
  print(subset(pooled, Taxa == "S. pulcher")[c("Resp", "scenario", "estimate", "ci_low", "ci_high", "q")])
  invisible(outputs)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) <= 1L)
  if (length(args)) run_sheephead_yoy_check(source_root = args[1]) else run_sheephead_yoy_check()
}
