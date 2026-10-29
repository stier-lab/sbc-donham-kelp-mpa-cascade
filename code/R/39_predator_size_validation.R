# Secondary size-pool validation and legacy-protection checks. Primary inputs stay fixed.
local({
  suppressPackageStartupMessages({library(dplyr); library(tidyr)})
  root <- here::here()
  source(file.path(root, "code/R/lobster_biomass_utils.R"), local = TRUE)
  primary_path <- file.path(root, "data/sumstats_final.csv")
  fingerprint <- tools::md5sum(primary_path)
  write <- function(d, stem) write.csv(d, file.path(root, "tables", paste0(stem, ".csv")), row.names = FALSE)
  primary <- read.csv(primary_path) %>% filter(Primary == "Y")
  primary$Mean <- as.numeric(primary$Mean)
  legacy <- c("Blue Cavern Onshore SMCA", "Farnsworth Onshore SMCA", "Long Point SMR")
  refit <- function(d, scenario) bind_rows(lapply(c("Bio", "Den"), function(response) {
    x <- filter(d, Resp == response)
    x$Taxa <- factor(x$Taxa)
    f <- metafor::rma.mv(Mean, V = SE^2, mods = ~Taxa - 1,
      random = list(~1 | MPA, ~1 | Source), method = "REML", test = "t", data = x)
    data.frame(scenario, Taxa = sub("^Taxa", "", names(coef(f))), Resp = response,
      estimate = as.numeric(coef(f)), ci_low = f$ci.lb, ci_high = f$ci.ub, p = f$pval)
  })) %>% mutate(q = p.adjust(p, "BH"))
  omit <- refit(filter(primary, !MPA %in% legacy), "omit_catalina_primary")
  base <- refit(primary, "retained_primary")
  joined <- omit %>% inner_join(base %>% select(Taxa, Resp, primary_estimate = estimate, primary_q = q), by = c("Taxa", "Resp")) %>%
    mutate(direction_changed = sign(estimate) != sign(primary_estimate), bh_changed = (q < 0.05) != (primary_q < 0.05))
  write(joined, "table_s_catalina_omission")
  write(primary %>% filter(MPA %in% legacy) %>% select(MPA, Taxa, Resp, Source), "table_s_catalina_membership")

  # Calendar slopes are descriptive summaries, without model-based P values or causal labels.
  slopes <- function(d, keys, outcome) {
    d %>% group_by(across(all_of(keys))) %>% group_modify(~ {
      stopifnot(!anyDuplicated(.x$year), all(is.finite(.x[[outcome]])))
      n <- nrow(.x)
      slope <- if (n >= 3 && diff(range(.x$year)) > 0) unname(coef(lm(.x[[outcome]] ~ .x$year))[2]) else NA_real_
      data.frame(n_years = n, first_year = min(.x$year), last_year = max(.x$year),
        slope_per_calendar_year = slope, inference = "Descriptive post-designation slope; no consumption or designation-effect inference")
    }) %>% ungroup()
  }
  for (species in c("sheephead", "lobster")) {
    path <- file.path(root, "tables", paste0("table_s_", species, "_size_annual.csv"))
    status_path <- file.path(root, "tables", if (species == "sheephead") "table_s_predator_size_status.csv" else "table_s_lobster_size_status.csv")
    if (file.exists(path) && file.exists(status_path) && all(read.csv(status_path)$status == "completed")) {
      d <- read.csv(path)
      keys <- c("CA_MPA_Name_Short", "measure", if (species == "sheephead") "threshold_cm")
      write(slopes(d, keys, "lnRR"), paste0("table_s_", species, "_size_temporal"))
    }
  }

  size_path <- file.path(root, "data/ALL_sizefreq_2024.csv")
  site_path <- file.path(root, "data/PISCO/master_site_table_Emilyedit.csv")
  if (file.exists(size_path) && file.exists(site_path)) {
    source(file.path(root, "code/R/00c_analysis_constants.R"), local = TRUE)
    source(file.path(root, "code/R/01_utils.R"), local = TRUE)
    sites <- read.csv(site_path) %>% filter(site_status %in% c("mpa", "reference"), !site %in% EXCLUDED_REFERENCE_SITES) %>%
      distinct(site, CA_MPA_Name_Short, site_status)
    stopifnot(!anyDuplicated(sites[c("site", "CA_MPA_Name_Short")]))
    sizes <- read.csv(size_path) %>% filter(classcode == "PANINT", is.finite(size), size > 0,
      is.finite(count), count > 0, month >= 5, month <= 10) %>%
      inner_join(sites, by = "site", relationship = "many-to-many")
    stopifnot(all(abs(sizes$count - round(sizes$count)) < 1e-10))
    set.seed(20261016)
    pools <- sizes %>% group_by(CA_MPA_Name_Short, year, status = site_status) %>% group_modify(~ {
      x <- .x
      n <- sum(x$count)
      draws <- rmultinom(1000, size = n, prob = x$count / n)
      fraction <- colSums(draws[x$size >= 82.55, , drop = FALSE]) / n
      mass <- lobster_mass_g(x$size)
      all_mass <- as.numeric(crossprod(mass, draws)) / n
      large_mass <- as.numeric(crossprod(mass * (x$size >= 82.55), draws)) / n
      interval <- function(z) unname(quantile(z, c(0.025, 0.975)))
      data.frame(n_sized = n, n_sites = n_distinct(x$site), n_size_classes = n_distinct(x$size),
        fraction = sum(x$count[x$size >= 82.55]) / n,
        fraction_low = interval(fraction)[1], fraction_high = interval(fraction)[2],
        mass = sum(x$count * mass) / n, mass_low = interval(all_mass)[1], mass_high = interval(all_mass)[2],
        large_mass = sum(x$count * mass * (x$size >= 82.55)) / n,
        large_mass_low = interval(large_mass)[1], large_mass_high = interval(large_mass)[2],
        uncertainty = "Size-pool-only multinomial resampling; conditional on measured individuals; excludes count uncertainty, survey clustering and detectability")
    }) %>% ungroup()
    write(pools, "table_s_lobster_size_pool_uncertainty")
    raw <- read.csv(file.path(root, "data/harmonized/harmonized_raw_responses.csv"))
    existing <- raw %>% filter(source == "PISCO", taxon_name == "Panulirus interruptus", resp == "Bio") %>%
      select(CA_MPA_Name_Short, year, status, primary_biomass = value)
    rebuilt <- read.csv(file.path(root, "tables/table_s_lobster_size_coverage.csv"))
    overlap <- inner_join(rebuilt, existing, by = c("CA_MPA_Name_Short", "year", "status")) %>%
      mutate(difference = all_biomass - primary_biomass,
        ratio = ifelse(primary_biomass > 0, all_biomass / primary_biomass, NA_real_))
    write(overlap, "table_s_lobster_biomass_reconciliation")
    write(data.frame(check = "Expanded lobster biomass", overlap_rows = nrow(overlap),
      median_absolute_difference = if (nrow(overlap)) median(abs(overlap$difference)) else NA_real_,
      status = "secondary_only", reason = "Conditional size-pool validation does not establish equivalence with transect-resampled primary biomass or propagate all survey uncertainty"),
      "table_s_lobster_biomass_validation")
    write(data.frame(path = c(size_path, site_path), md5 = unname(tools::md5sum(c(size_path, site_path)))), "table_s_size_validation_inputs")
  }
  stopifnot(identical(fingerprint, tools::md5sum(primary_path)))
  cat("Secondary size validation and Catalina omission complete; primary cohort unchanged.\n")
})
