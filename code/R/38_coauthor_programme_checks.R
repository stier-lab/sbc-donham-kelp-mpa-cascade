# Secondary checks requested during coauthor review. Never overwrite primary inputs.
local({
  suppressPackageStartupMessages({library(dplyr); library(tidyr)})
  root <- here::here()
  source(file.path(root, "code/R/coauthor_sensitivity_utils.R"), local = TRUE)
  source(file.path(root, "code/R/00c_analysis_constants.R"), local = TRUE)
  source(file.path(root, "code/R/01_utils.R"), local = TRUE)
  source(file.path(root, "code/R/lobster_biomass_utils.R"), local = TRUE)
  primary <- read.csv(file.path(root, "data/sumstats_final.csv")) %>% filter(Primary == "Y")
  primary$Mean <- as.numeric(primary$Mean)
  primary$vi <- as.numeric(primary$SE)^2
  fingerprint <- tools::md5sum(file.path(root, "data/sumstats_final.csv"))
  write <- function(d, name) write.csv(d, file.path(root, "tables", paste0(name, ".csv")), row.names = FALSE)

  # Refit the unchanged joint model after omissions; one-source models need no Source RE.
  joint <- function(d, scenario) bind_rows(lapply(c("Bio", "Den"), function(response) {
    x <- d[d$Resp == response, ]
    x$Taxa <- factor(x$Taxa); x$Source <- factor(x$Source)
    re <- if (nlevels(x$Source) > 1) list(~1 | MPA, ~1 | Source) else ~1 | MPA
    fit <- metafor::rma.mv(Mean, V = vi, mods = ~Taxa - 1, random = re,
                           method = "REML", test = "t", data = x)
    taxa <- sub("^Taxa", "", names(coef(fit)))
    data.frame(scenario, Taxa = taxa, Resp = response,
      k = vapply(taxa, function(t) sum(x$Taxa == t), 1L),
      n_mpa = vapply(taxa, function(t) n_distinct(x$MPA[x$Taxa == t]), 1L),
      estimate = as.numeric(coef(fit)), se = fit$se, ci_low = fit$ci.lb,
      ci_high = fit$ci.ub, p = fit$pval)
  })) %>% mutate(q = p.adjust(p, "BH"))
  scenarios <- list(retained_primary = primary,
    omit_kfm_sheephead = filter(primary, !(Source == "KFM" & Taxa == "S. pulcher")),
    omit_skunk_point = filter(primary, !grepl("Skunk", MPA)),
    omit_gull_kfm = filter(primary, !(Source == "KFM" & MPA == "Gull Island SMR")),
    pisco_only = filter(primary, Source == "PISCO"))
  omissions <- bind_rows(lapply(names(scenarios), function(n) joint(scenarios[[n]], n)))
  baseline <- omissions %>% filter(scenario == "retained_primary") %>%
    inner_join(read.csv(file.path(root, "tables/table_02_meta_analysis.csv")) %>%
      mutate(Resp = ifelse(Response == "Biomass", "Bio", "Den")), by = c("Taxa", "Resp"))
  stopifnot(nrow(baseline) == 9L, max(abs(baseline$estimate - baseline$Estimate)) < 1e-7)
  write(omissions, "table_s_programme_omissions")
  write(primary %>% group_by(MPA, Taxa, Resp) %>% filter(n_distinct(Source) > 1) %>%
    mutate(mpa_mean = mean(Mean), deviation_from_mean = Mean - mpa_mean) %>% ungroup(),
    "table_s_programme_overlap")

  # Additional regressions use the same aggregation as Table 3, not extra MPA replication.
  pooled <- primary %>% group_by(MPA, Taxa, Resp) %>%
    summarise(estimate = mean(Mean), variance = sum(vi) / n()^2,
      sources = paste(sort(unique(Source)), collapse = ";"), .groups = "drop")
  links <- bind_rows(expand.grid(predator = c("P. interruptus", "S. pulcher"),
      prey = c("S. purpuratus", "M. franciscanus"), predictor_response = "Bio",
      outcome_response = "Bio", stringsAsFactors = FALSE),
    expand.grid(predator = c("P. interruptus", "S. pulcher"), prey = "M. pyrifera",
      predictor_response = c("Bio", "Den"), outcome_response = "Bio", stringsAsFactors = FALSE))
  regression_pairs <- list()
  regressions <- bind_rows(lapply(seq_len(nrow(links)), function(i) {
    l <- links[i, ]
    x <- pooled %>% filter(Taxa == l$predator, Resp == l$predictor_response) %>%
      transmute(MPA, x = estimate, predictor_sources = sources)
    y <- pooled %>% filter(Taxa == l$prey, Resp == l$outcome_response) %>%
      transmute(MPA, y = estimate, vi = variance, outcome_sources = sources)
    d <- inner_join(x, y, by = "MPA")
    regression_pairs[[i]] <<- cbind(link = i, d)
    fit <- if (nrow(d) >= 3 && sd(d$x) > 0) metafor::rma(y, vi = vi, mods = ~x,
      method = "REML", test = "knha", data = d) else NULL
    cbind(l, n_mpa = nrow(d), data.frame(slope = if (is.null(fit)) NA_real_ else coef(fit)[2],
      ci_low = if (is.null(fit)) NA_real_ else fit$ci.lb[2],
      ci_high = if (is.null(fit)) NA_real_ else fit$ci.ub[2],
      p = if (is.null(fit)) NA_real_ else fit$pval[2]))
  })) %>% mutate(q = p.adjust(p, "BH"),
    limitation = "Exploratory; mixed primary contrasts, uncertain predictors, shared controls; not a mediation test")
  write(regressions, "table_s_extended_pathways")
  write(bind_rows(regression_pairs), "table_s_extended_pathway_pairs")

  # Compare PISCO/KFM sheephead on identical calendar years, not different time spans.
  raw <- read.csv(file.path(root, "data/harmonized/harmonized_raw_responses.csv"))
  fish_keys <- primary %>% filter(Taxa == "S. pulcher", Resp == "Den") %>% distinct(MPA, Source)
  annual <- raw %>% filter(taxon_name == "Semicossyphus pulcher", resp == "Den", BA == "After") %>%
    inner_join(fish_keys, by = c("CA_MPA_Name_Short" = "MPA", "source" = "Source")) %>%
    pivot_wider(id_cols = c(CA_MPA_Name_Short, source, year), names_from = status, values_from = value) %>%
    filter(is.finite(mpa), is.finite(reference)) %>% group_by(CA_MPA_Name_Short, source) %>%
    mutate(lnRR = coauthor_offset_ratio(mpa, reference)) %>% ungroup()
  common <- annual %>% filter(source %in% c("PISCO", "KFM")) %>%
    group_by(CA_MPA_Name_Short, year) %>% filter(n_distinct(source) == 2) %>% ungroup()
  write(common, "table_s_sheephead_programme_common_years")
  write(common %>% group_by(CA_MPA_Name_Short, source) %>% summarise(n_years = n(),
    years = paste(sort(year), collapse = ","), mean_lnRR = mean(lnRR), .groups = "drop"),
    "table_s_sheephead_programme_comparison")
  island_mpas <- primary %>% filter(Location == "C", !grepl("Onshore|Long Point", MPA)) %>% pull(MPA) %>% unique()
  write(annual %>% filter(source == "PISCO", CA_MPA_Name_Short %in% island_mpas) %>%
    group_by(CA_MPA_Name_Short) %>% summarise(n_years = n(), mean_lnRR = mean(lnRR),
      year_min = min(year), year_max = max(year), .groups = "drop"), "table_s_pisco_island_sheephead")

  # Raw-size checks are optional in a public derived-data checkout; outputs say so explicitly.
  fish_path <- file.path(root, "data/PISCO/MLPA_kelpforest_fish_2024.csv")
  site_path <- file.path(root, "data/PISCO/master_site_table_Emilyedit.csv")
  if (file.exists(fish_path) && file.exists(site_path)) {
    fish <- read.csv(fish_path)
    sites <- read.csv(site_path) %>% filter(site_status %in% c("mpa", "reference")) %>%
      distinct(site, CA_MPA_Name_Short, site_status, BaselineRegion)
    stopifnot(!anyDuplicated(sites[c("site", "CA_MPA_Name_Short")]))
    fish <- fish %>% filter(level != "CAN", !is.na(count), count >= 0,
      !site %in% EXCLUDED_REFERENCE_SITES)
    # All surveyed transects supply zeros, including those with no sheephead.
    roster <- fish %>% distinct(site, year, zone, transect)
    sheep <- fish %>% filter(classcode == "SPUL")
    results <- lapply(c(0, 30), function(threshold) {
      counts <- sheep %>% group_by(site, year, zone, transect) %>%
        summarise(missing_size = any(count > 0 & (!is.finite(fish_tl) | fish_tl <= 0)),
          density = sum(ifelse(is.finite(fish_tl) & fish_tl > threshold, count, 0)),
          biomass = sum(ifelse(is.finite(fish_tl) & fish_tl > threshold,
            count * 0.0144 * fish_tl^3.04, 0)), .groups = "drop")
      roster %>% left_join(counts, by = c("site", "year", "zone", "transect")) %>%
        mutate(missing_size = replace_na(missing_size, FALSE),
          density = replace_na(density, 0), biomass = replace_na(biomass, 0)) %>%
        group_by(site, year) %>% filter(!any(missing_size)) %>% ungroup() %>%
        group_by(site, year, zone) %>% summarise(density = mean(density), biomass = mean(biomass), .groups = "drop") %>%
        group_by(site, year) %>% summarise(density = mean(density) / PISCO_SWATH_AREA_M2,
          biomass = mean(biomass) / PISCO_SWATH_AREA_M2, .groups = "drop") %>%
        inner_join(sites, by = "site", relationship = "many-to-many") %>%
        filter(BaselineRegion == "SOUTH", site_status %in% c("mpa", "reference")) %>%
        inner_join(fish_keys %>% filter(Source == "PISCO"), by = c("CA_MPA_Name_Short" = "MPA")) %>%
        group_by(CA_MPA_Name_Short, year, site_status) %>%
        summarise(density = mean(density), biomass = mean(biomass), .groups = "drop") %>%
        mutate(threshold_cm = threshold)
    })
    sized <- bind_rows(results) %>% inner_join(read.csv(file.path(root, "data/harmonized/harmonized_site_metadata.csv")) %>%
      distinct(CA_MPA_Name_Short, MPA_Start), by = "CA_MPA_Name_Short") %>% filter(year > MPA_Start)
    size_pairs <- sized %>% pivot_longer(c(density, biomass), names_to = "measure", values_to = "value") %>%
      pivot_wider(id_cols = c(CA_MPA_Name_Short, year, threshold_cm, measure), names_from = site_status, values_from = value) %>%
      filter(is.finite(mpa), is.finite(reference)) %>% group_by(CA_MPA_Name_Short, measure, year) %>%
      filter(n_distinct(threshold_cm) == 2) %>% ungroup() %>%
      group_by(CA_MPA_Name_Short, measure) %>% mutate(lnRR = coauthor_offset_ratio(mpa, reference)) %>% ungroup()
    check_size <- size_pairs %>% select(CA_MPA_Name_Short, year, measure, threshold_cm, mpa, reference) %>%
      pivot_wider(names_from = threshold_cm, values_from = c(mpa, reference))
    stopifnot(nrow(check_size) > 0, all(check_size$mpa_30 <= check_size$mpa_0 + 1e-10),
      all(check_size$reference_30 <= check_size$reference_0 + 1e-10),
      all(size_pairs$mpa >= 0), all(size_pairs$reference >= 0))
    write(size_pairs, "table_s_sheephead_size_annual")
    write(size_pairs %>% group_by(CA_MPA_Name_Short, measure, threshold_cm) %>%
      summarise(n_years = n(), mean_lnRR = mean(lnRR), years = paste(sort(year), collapse = ","), .groups = "drop"),
      "table_s_sheephead_size_comparison")
    write(data.frame(path = c(fish_path, site_path), md5 = unname(tools::md5sum(c(fish_path, site_path)))),
      "table_s_predator_size_inputs")
    write(data.frame(check = "PISCO sheephead >30 cm TL", status = "completed",
      limitation = "Post-only descriptive matched-year ratios, not a refit of primary designation contrasts; positive-count missing-size site-years omitted; threshold not a universal feeding or legal cutoff"),
      "table_s_predator_size_status")
  } else {
    write(data.frame(check = "PISCO sheephead >30 cm TL", status = "raw inputs unavailable",
      limitation = "Obtain contributing-programme raw survey release to reproduce size check"), "table_s_predator_size_status")
  }
  size_path <- file.path(root, "data/ALL_sizefreq_2024.csv")
  if (file.exists(size_path) && file.exists(site_path)) {
    sizes <- read.csv(size_path) %>% filter(classcode == "PANINT", is.finite(size), size > 0,
      is.finite(count), count > 0, month >= 5, month <= 10)
    lobster_sites <- read.csv(site_path) %>% filter(site_status %in% c("mpa", "reference"),
      !site %in% EXCLUDED_REFERENCE_SITES) %>% distinct(site, CA_MPA_Name_Short, site_status)
    stopifnot(!anyDuplicated(lobster_sites[c("site", "CA_MPA_Name_Short")]))
    # Expected large-animal fraction/mass from measured size classes, not a new primary biomass export.
    pools <- sizes %>% inner_join(lobster_sites, by = "site", relationship = "many-to-many") %>%
      group_by(CA_MPA_Name_Short, year, status = site_status) %>% summarise(n_sized = sum(count),
        large_fraction = sum(count[size >= 82.55]) / sum(count),
        all_mass = sum(count * lobster_mass_g(size)) / sum(count),
        large_mass = sum(count * lobster_mass_g(size) * (size >= 82.55)) / sum(count), .groups = "drop")
    lobster_keys <- primary %>% filter(Taxa == "P. interruptus", Resp == "Den", Source == "PISCO") %>%
      distinct(MPA)
    lobster <- raw %>% filter(source == "PISCO", taxon_name == "Panulirus interruptus",
      resp == "Den", BA == "After") %>% inner_join(lobster_keys, by = c("CA_MPA_Name_Short" = "MPA")) %>%
      left_join(pools, by = c("CA_MPA_Name_Short", "year", "status")) %>%
      filter(is.finite(large_fraction) | value == 0) %>%
      mutate(all_density = value, large_density = ifelse(value == 0, 0, value * large_fraction),
        all_biomass = ifelse(value == 0, 0, value * all_mass),
        large_biomass = ifelse(value == 0, 0, value * large_mass))
    write(lobster, "table_s_lobster_size_coverage")
    lobster_pairs <- lobster %>% pivot_longer(c(all_density, large_density, all_biomass, large_biomass),
      names_to = "measure", values_to = "value_sized") %>%
      pivot_wider(id_cols = c(CA_MPA_Name_Short, year, measure), names_from = status, values_from = value_sized) %>%
      filter(is.finite(mpa), is.finite(reference)) %>%
      mutate(response = ifelse(grepl("density", measure), "density", "biomass")) %>%
      group_by(CA_MPA_Name_Short, response) %>%
      mutate(lnRR = coauthor_offset_ratio(mpa, reference)) %>% ungroup()
    write(lobster_pairs, "table_s_lobster_size_annual")
    write(lobster_pairs %>% group_by(CA_MPA_Name_Short, measure) %>%
      summarise(n_years = n(), mean_lnRR = mean(lnRR), years = paste(sort(year), collapse = ","), .groups = "drop"),
      "table_s_lobster_size_comparison")
    write(data.frame(path = c(size_path, site_path), md5 = unname(tools::md5sum(c(size_path, site_path)))),
      "table_s_lobster_size_inputs")
    write(data.frame(check = "PISCO lobster >=82.55 mm CL", status = "completed",
      limitation = "Descriptive size-pool reconstruction; no size-pool uncertainty or primary biomass substitution"),
      "table_s_lobster_size_status")
  } else {
    write(data.frame(check = "PISCO lobster >=82.55 mm CL", status = "raw inputs unavailable",
      limitation = "Do not reuse previously generated size outputs as a fresh calculation"),
      "table_s_lobster_size_status")
  }
  stopifnot(identical(fingerprint, tools::md5sum(file.path(root, "data/sumstats_final.csv"))))
  print(regressions)
  cat("Programme and exploratory pathway checks complete; primary contrasts unchanged.\n")
})
