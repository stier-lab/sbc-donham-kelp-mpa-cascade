# Programme-separated pathway sensitivity; no primary estimate or figure is replaced.
(function() {
  suppressPackageStartupMessages({library(dplyr); library(tidyr)})
  root <- here::here()
  for (file in c("primary_uncertainty_utils.R", "programme_pathway_utils.R"))
    source(file.path(root, "code/R", file), local = TRUE)
  write <- function(d, stem) write.csv(d, file.path(root, "tables", paste0("table_s_programme_pathway_", stem, ".csv")), row.names = FALSE)
  completed <- FALSE
  on.exit(if (!completed) write(data.frame(status = "failed", reason = "Do not reuse stale outputs as a completed run"), "status"), add = TRUE)
  protected <- file.path(root, c("data/sumstats_final.csv", "tables/table_02_meta_analysis.csv",
    "tables/table_03_cross_taxa_meta_regression.csv", paste0("plots/fig_0", 1:4, c("_mpa_map.png", "_cascade_case_studies.png", "_mean_effects.png", "_recovery_curves.png"))))
  before <- tools::md5sum(protected)
  input_paths <- c("data/harmonized/harmonized_response_ratios.csv", "data/sumstats_final.csv",
    "tables/table_s_design_reference_sites.csv", "code/R/primary_uncertainty_utils.R",
    "code/R/programme_pathway_utils.R", "code/R/42_programme_pathway_refits.R")
  primary <- read.csv(file.path(root, "data/sumstats_final.csv")) %>% filter(Primary == "Y")
  annual <- programme_pathway_annual(read.csv(file.path(root, input_paths[1])), primary)
  registry <- read.csv(file.path(root, input_paths[3]))
  links <- programme_pathway_links()
  source_pairs <- list(KFM_PISCO = c("KFM", "PISCO"), LTER_PISCO = c("LTER", "PISCO"), KFM_LTER = c("KFM", "LTER"))
  summaries <- effects <- observations <- coverage <- omissions <- comparisons <- groups <- list()
  models <- list()
  for (scope in c("retained_primary_histories", "available_diver_histories")) {
    for (pair in names(source_pairs)) for (i in seq_len(nrow(links))) {
      link <- links[i, ]; sources <- source_pairs[[pair]]
      tag <- cbind(scope, programme_pair = pair, link)
      rownames(tag) <- NULL
      matched <- programme_pathway_match(annual, link, sources, scope)
      key <- paste(scope, pair, link$link, sep = "|")
      if (nrow(matched$coverage)) coverage[[key]] <- cbind(tag, matched$coverage)
      if (nrow(matched$annual)) observations[[key]] <- cbind(tag, matched$annual)
      refs <- programme_pathway_reference_groups(sort(unique(matched$annual$MPA)), sources, registry)
      if (nrow(refs)) groups[[key]] <- cbind(tag, refs)
      for (lag in c(1L, 3L)) {
        e <- programme_pathway_effects(matched$annual, lag)
        if (nrow(e)) effects[[paste(key, lag)]] <- cbind(tag, lag_years = lag, e)
        rowset <- list()
        for (programme in c(sources, "equal_programme_mean")) {
          d <- if (nrow(e)) e[e$Source == programme, ] else data.frame()
          fitted <- programme_pathway_fit(d)
          row <- cbind(tag, lag_years = lag, Source = programme, fitted$result)
          summaries[[paste(key, lag, programme)]] <- row
          rowset[[programme]] <- row
          if (!is.null(fitted$model)) {
            models[[paste(key, lag, programme)]] <- fitted$model
            omit_sets <- c(setNames(as.list(d$MPA), paste0("MPA|", d$MPA)),
              setNames(lapply(unique(refs$reference_group), function(g) refs$MPA[refs$reference_group == g]),
                paste0("reference_group|", unique(refs$reference_group))))
            for (name in names(omit_sets)) {
              refit <- programme_pathway_fit(d[!d$MPA %in% omit_sets[[name]], ])
              omissions[[paste(key, lag, programme, name)]] <- cbind(tag, lag_years = lag, Source = programme,
                omitted = name, refit$result)
            }
          }
        }
        comparisons[[paste(key, lag)]] <- cbind(tag, lag_years = lag,
          n_mpa = rowset[[1]]$n_mpa,
          slope_second_minus_first = rowset[[sources[2]]]$slope - rowset[[sources[1]]]$slope,
          quantity = "Descriptive paired slope difference; no independent-programme test; not a survey-method effect")
      }
    }
  }
  summary <- bind_rows(summaries) %>% group_by(scope, programme_pair, lag_years, Source) %>%
    mutate(q = p.adjust(p, "BH", n = 6L)) %>% ungroup() %>%
    mutate(quantity = "Mean annual post-designation log MPA/reference ratio on identical paired years",
      uncertainty = "Calendar-gap HAC for annual means; conditional Knapp-Hartung REML; predictor error and unverified spatial covariance not fitted")
  write(summary, "fits")
  write(bind_rows(effects), "effects")
  write(bind_rows(observations), "annual")
  write(bind_rows(coverage), "coverage")
  omission_table <- bind_rows(omissions) %>% group_by(scope, programme_pair, lag_years, Source, omitted) %>%
    mutate(q = p.adjust(p, "BH", n = 6L)) %>% ungroup()
  write(omission_table, "omissions")
  write(bind_rows(comparisons), "differences")
  write(bind_rows(groups), "reference_groups")
  write(data.frame(path = input_paths, md5 = unname(tools::md5sum(file.path(root, input_paths)))), "inputs")
  packages <- c("dplyr", "tidyr", "metafor", "sandwich")
  write(data.frame(package = packages, version = vapply(packages, function(p) as.character(packageVersion(p)), "")), "versions")
  saveRDS(models, file.path(root, "data/cache/programme_pathway_refits.rds"))
  stopifnot(identical(before, tools::md5sum(protected)))
  write(data.frame(status = "completed", primary_changed = FALSE,
    limitation = "Matched MPA/year sensitivity, not identical transects or a causal programme comparison; sparse links explicitly unestimable"), "status")
  completed <- TRUE
  print(summary %>% filter(lag_years == 3, programme_pair == "KFM_PISCO") %>%
    select(scope, link, Source, n_mpa, status, slope, ci_low, ci_high, p, q))
})()
