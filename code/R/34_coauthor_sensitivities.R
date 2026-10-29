# Reproduce with Rscript code/R/34_coauthor_sensitivities.R; no raw checkout needed.
local({
  library(dplyr)
  library(tidyr)
  root <- here::here()
  source(file.path(root, "code/R/coauthor_sensitivity_utils.R"), local = TRUE)
  primary <- read.csv(file.path(root, "data/sumstats_final.csv"))
  primary <- primary[primary$Primary == "Y", ]
  primary$Mean <- as.numeric(primary$Mean)
  primary$vi <- as.numeric(primary$SE)^2
  naples <- primary$MPA == "Naples SMCA" & primary$Source == "LTER" &
    primary$Taxa == "S. pulcher" & primary$Resp == "Bio"
  stopifnot(nrow(primary) == 145L, sum(naples) == 1L,
    primary$Model[naples] == "Linear", primary$AnalysisType[naples] == "pBACIPS",
    primary$N[naples] == 23L, all(is.finite(primary$vi)), all(primary$vi > 0))
  write_table <- function(x, name) write.csv(x, file.path(root, "tables", paste0(name, ".csv")), row.names = FALSE)
  input_files <- c("data/sumstats_final.csv", "data/harmonized/harmonized_raw_responses.csv",
                   "data/harmonized/reference_sites.csv")
  fingerprints <- unname(tools::md5sum(file.path(root, input_files)))

  fit_joint <- function(d, scenario) {
    bind_rows(lapply(c("Bio", "Den"), function(response) {
      x <- d[d$Resp == response, ]
      x$Taxa <- factor(x$Taxa); x$Source <- factor(x$Source)
      warnings <- character()
      fit <- withCallingHandlers(metafor::rma.mv(Mean, V = vi, mods = ~ Taxa - 1,
        random = list(~1 | MPA, ~1 | Source), data = x, method = "REML", test = "t"),
        warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
      taxa <- sub("^Taxa", "", names(stats::coef(fit)))
      data.frame(scenario, Taxa = taxa, Resp = response,
        k = vapply(taxa, function(t) sum(x$Taxa == t), integer(1)),
        n_mpa = vapply(taxa, function(t) length(unique(x$MPA[x$Taxa == t])), integer(1)),
        estimate = as.numeric(stats::coef(fit)), se = fit$se, ci_low = fit$ci.lb,
        ci_high = fit$ci.ub, p = fit$pval, warnings = paste(unique(warnings), collapse = "; "))
    })) %>% mutate(q = p.adjust(p, "BH"))
  }
  gull <- primary$MPA == "Gull Island SMR" & primary$Source == "KFM"
  stopifnot(sum(gull) == 4L, sum(gull & primary$Taxa == "P. interruptus") == 1L)
  scenarios <- list(retained_primary = primary, omit_gull_kfm = primary[!gull, ],
                    omit_gull_kfm_lobster = primary[!(gull & primary$Taxa == "P. interruptus"), ],
                    omit_naples_lter_sheephead_biomass = primary[!naples, ])
  omitted <- bind_rows(lapply(names(scenarios), function(name) fit_joint(scenarios[[name]], name)))
  published <- read.csv(file.path(root, "tables/table_02_meta_analysis.csv")) %>%
    mutate(Resp = ifelse(Response == "Biomass", "Bio", "Den"))
  baseline <- omitted %>% filter(scenario == "retained_primary") %>%
    inner_join(published, by = c("Taxa", "Resp"))
  stopifnot(nrow(baseline) == 9L,
    max(abs(baseline$estimate - baseline$Estimate)) < 1e-7,
    max(abs(baseline$q - baseline$pval_fdr)) < 1e-7)
  write_table(omitted, "table_s_coauthor_gull_omission")
  write_table(primary[gull, ], "table_s_coauthor_gull_membership")
  write_table(omitted[omitted$scenario %in% c("retained_primary", "omit_naples_lter_sheephead_biomass"), ],
              "table_s_naples_biomass_inclusion")
  write_table(primary[naples, ], "table_s_naples_biomass_history")

  raw <- read.csv(file.path(root, "data/harmonized/harmonized_raw_responses.csv"))
  annual <- bind_rows(lapply(c("Bio", "Den"), function(response) coauthor_predator_cohort(raw, primary, response)))
  write_table(annual, "table_s_coauthor_combined_predator_annual")
  source_pairs <- annual %>% filter(is.finite(predator_lnRR), is.finite(urchin_lnRR)) %>%
    group_by(response, MPA, source) %>% summarise(n_common_years = n(),
      year_min = min(year), year_max = max(year), years = paste(sort(year), collapse = ","),
      predator_lnRR = mean(predator_lnRR), urchin_lnRR = mean(urchin_lnRR), .groups = "drop")
  write_table(source_pairs, "table_s_coauthor_combined_predator_coverage")
  pairs <- source_pairs %>% group_by(response, MPA) %>% summarise(n_sources = n(),
    sources = paste(sort(source), collapse = ";"), predator_lnRR = mean(predator_lnRR),
    urchin_lnRR = mean(urchin_lnRR), .groups = "drop")
  write_table(pairs, "table_s_coauthor_combined_predator_pairs")
  tests <- bind_rows(lapply(c("Bio", "Den"), function(response) {
    d <- pairs[pairs$response == response, ]
    cbind(response, coauthor_spearman(d$predator_lnRR, d$urchin_lnRR))
  })) %>% mutate(q = p.adjust(p, "BH"),
    limitation = "Post-only common-cohort levels, not designation changes; shared references and estimated predictors; not a causal or consumptive-capacity test")
  write_table(tests, "table_s_coauthor_combined_predator_tests")

  kelp_keys <- primary %>% filter(Taxa == "M. pyrifera", Source != "Landsat") %>% distinct(MPA, Source)
  kelp <- raw %>% filter(taxon_name == "Macrocystis pyrifera", resp == "Den") %>%
    inner_join(kelp_keys, by = c("CA_MPA_Name_Short" = "MPA", "source" = "Source")) %>%
    group_by(source, CA_MPA_Name_Short) %>% summarise(n_years = n_distinct(year),
      year_min = min(year), year_max = max(year), .groups = "drop") %>%
    mutate(measure = "harmonized adult stipe density, not stage-specific plant density",
           juvenile_density_available = FALSE, adult_plant_density_available = FALSE)
  write_table(kelp, "table_s_coauthor_kelp_density_coverage")
  stopifnot(identical(fingerprints, unname(tools::md5sum(file.path(root, input_files)))))
  write_table(data.frame(path = input_files, md5 = fingerprints), "table_s_coauthor_sensitivity_inputs")
  cat("Coauthor checks: 145 primary keys retained; Naples biomass exclusion, Gull Island omissions and common-cohort predator checks exported.\n")
  print(tests)
})
