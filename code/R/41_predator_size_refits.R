# Actual matched size-specific temporal and predator-urchin refits; primary inputs fixed.
(function() {
  suppressPackageStartupMessages({library(dplyr); library(tidyr)})
  root <- here::here()
  for (file in c("coauthor_sensitivity_utils.R", "primary_uncertainty_utils.R",
                 "temporal_uncertainty_utils.R", "predator_size_refit_utils.R"))
    source(file.path(root, "code/R", file), local = TRUE)
  write <- function(d, stem) write.csv(d, file.path(root, "tables", paste0(stem, ".csv")), row.names = FALSE)
  status_recorded <- FALSE
  on.exit(if (!status_recorded) write(data.frame(status = "failed",
    reason = "Refit did not complete; do not treat earlier outputs as a fresh run"), "table_s_size_refit_status"), add = TRUE)
  protected <- file.path(root, c("data/sumstats_final.csv", "tables/table_02_meta_analysis.csv",
    "tables/table_03_cross_taxa_meta_regression.csv", "plots/fig_04_recovery_curves.png"))
  initial <- tools::md5sum(protected)
  inputs <- file.path(root, c("tables/table_s_sheephead_size_annual.csv",
    "tables/table_s_lobster_size_annual.csv", "data/harmonized/harmonized_response_ratios.csv",
    "data/harmonized/harmonized_site_metadata.csv", "data/sumstats_final.csv",
    "code/R/41_predator_size_refits.R", "code/R/predator_size_refit_utils.R",
    "code/R/predator_size_refit_display.R"))
  status_paths <- file.path(root, paste0("tables/table_s_", c("predator", "lobster"), "_size_status.csv"))
  if (!all(file.exists(inputs)) || !all(file.exists(status_paths)) ||
      !all(vapply(status_paths, function(p) all(read.csv(p)$status == "completed"), TRUE))) {
    write(data.frame(status = "unavailable", reason = "Current verified size inputs required; no primary changes or fresh refit claim"),
      "table_s_size_refit_status")
    message("Size refits skipped: regenerate verified size inputs with modules 38/39.")
    status_recorded <- TRUE
    return(invisible(NULL))
  }
  for (stem in c("predator", "lobster")) {
    status <- read.csv(file.path(root, paste0("tables/table_s_", stem, "_size_status.csv")))
    stopifnot(all(status$status == "completed"))
    saved <- read.csv(file.path(root, paste0("tables/table_s_", stem, "_size_inputs.csv")))
    if (!all(file.exists(saved$path))) {
      write(data.frame(status = "unavailable", reason = "Upstream raw-source paths unavailable; run modules 38/39 in this checkout"),
        "table_s_size_refit_status")
      status_recorded <- TRUE
      return(invisible(NULL))
    }
    stopifnot(identical(unname(tools::md5sum(saved$path)), as.character(saved$md5)))
  }
  primary <- read.csv(file.path(root, "data/sumstats_final.csv")) %>% filter(Primary == "Y")
  ratios <- read.csv(file.path(root, "data/harmonized/harmonized_response_ratios.csv"))
  resource <- temporal_primary_data(ratios, primary) %>% filter(source == "PISCO")
  taxa <- c("Panulirus interruptus" = "P. interruptus", "Semicossyphus pulcher" = "S. pulcher",
    "Strongylocentrotus purpuratus" = "S. purpuratus", "Mesocentrotus franciscanus" = "M. franciscanus",
    "Macrocystis pyrifera" = "M. pyrifera")
  resource$Species <- unname(taxa[resource$y])
  predators <- size_refit_annual(root)
  slopes <- diagnostics <- attempts <- annual <- predictions <- ar1_audits <- shared <- list()
  models <- list()
  for (scope in c("retained_membership", "expanded_size_cohort")) {
    for (response in c("Bio", "Den")) {
      cohorts <- lapply(c("all", "large"), function(threshold)
        size_refit_temporal_data(predators, resource, primary, scope, threshold, response))
      stopifnot(identical(cohorts[[1]][c("history", "year", "time")],
        cohorts[[2]][c("history", "year", "time")]))
      for (j in 1:2) {
        threshold <- c("all", "large")[j]; d <- cohorts[[j]]
        tag <- paste(scope, response, threshold, sep = "|")
        message("Size refit: ", tag, "; n=", nrow(d))
        fitted <- size_refit_temporal(d)
        m <- fitted$model; models[[tag]] <- m
        attempts[[tag]] <- cbind(scope, Resp = response, size = threshold, fitted$audit)
        annual[[tag]] <- cbind(scope, d)
        slopes[[tag]] <- bind_rows(lapply(levels(d$Species), function(sp)
          cbind(scope, Resp = response, size = threshold, Species = sp, model = "independent_errors",
            n_obs = sum(d$Species == sp), n_mpa = n_distinct(d$MPA[d$Species == sp]),
            temporal_contrast(m, sp, slope = TRUE))))
        predictions[[tag]] <- cbind(scope, size = threshold, model = "independent_errors",
          temporal_prediction_grid(m, response, times = seq(0, 15, by = .5)))
        residual <- residuals(m)
        diagnostics[[tag]] <- data.frame(scope, Resp = response, size = threshold,
          n_obs = nrow(d), n_history = n_distinct(d$history), n_mpa = n_distinct(d$MPA),
          candidate = fitted$candidate, converged = TRUE, singular = lme4::isSingular(m),
          residual_sd = sigma(m), residual_shapiro_p = shapiro.test(residual)$p.value,
          residual_spread_rho = suppressWarnings(cor(abs(residual), fitted(m), method = "spearman")))
        {
          correlated <- size_refit_ar1(d)
          ar1_audits[[tag]] <- cbind(scope, Resp = response, size = threshold, correlated$audit)
          if (!is.null(correlated$model)) {
            models[[paste0(tag, "|AR1")]] <- correlated$model
            slopes[[paste0(tag, "|AR1")]] <- bind_rows(lapply(levels(d$Species), function(sp)
              cbind(scope, Resp = response, size = threshold, Species = sp, model = "calendar_AR1",
                n_obs = sum(d$Species == sp), n_mpa = n_distinct(d$MPA[d$Species == sp]),
                size_refit_ar1_slope(correlated$model, sp))))
            predictions[[paste0(tag, "|AR1")]] <- cbind(scope, size = threshold, model = "calendar_AR1",
              size_refit_ar1_predictions(correlated$model, response))
          }
        }
        slopes[[paste0(tag, "|CR2")]] <- bind_rows(lapply(levels(d$Species), function(sp)
          cbind(scope, Resp = response, size = threshold, Species = sp, model = "MPA_cluster_CR2",
            n_obs = sum(d$Species == sp), n_mpa = n_distinct(d$MPA[d$Species == sp]),
            size_refit_cluster_slope(m, sp))))
        sub <- filter(d, !MPA %in% c("Naples SMCA", "Campus Point SMCA"))
        without_shared <- size_refit_temporal(sub)
        shared[[tag]] <- bind_rows(lapply(levels(sub$Species), function(sp)
          cbind(scope, Resp = response, size = threshold, Species = sp,
            n_mpa = n_distinct(sub$MPA[sub$Species == sp]),
            size_refit_cluster_slope(without_shared$model, sp))))
      }
      # Paired threshold difference models preserve all/large dependence within a year.
      delta <- cohorts[[1]] %>% filter(Species %in% c("P. interruptus", "S. pulcher")) %>%
        select(-size) %>% inner_join(cohorts[[2]] %>% select(history, year, large_lnDiff = lnDiff),
          by = c("history", "year")) %>% mutate(lnDiff = large_lnDiff - lnDiff, Species = droplevels(Species))
      tag <- paste(scope, response, "paired_threshold_difference", sep = "|")
      fitted <- size_refit_temporal(delta)
      models[[tag]] <- fitted$model
      attempts[[tag]] <- cbind(scope, Resp = response, size = "paired_threshold_difference", fitted$audit)
      slopes[[tag]] <- bind_rows(lapply(levels(delta$Species), function(sp)
        cbind(scope, Resp = response, size = "paired_threshold_difference", Species = sp,
          model = "independent_errors", n_obs = sum(delta$Species == sp),
          n_mpa = n_distinct(delta$MPA[delta$Species == sp]), temporal_contrast(fitted$model, sp, slope = TRUE))))
      correlated <- size_refit_ar1(delta)
      ar1_audits[[tag]] <- cbind(scope, Resp = response, size = "paired_threshold_difference", correlated$audit)
      if (!is.null(correlated$model)) {
        models[[paste0(tag, "|AR1")]] <- correlated$model
        slopes[[paste0(tag, "|AR1")]] <- bind_rows(lapply(levels(delta$Species), function(sp)
          cbind(scope, Resp = response, size = "paired_threshold_difference", Species = sp, model = "calendar_AR1",
            n_obs = sum(delta$Species == sp), n_mpa = n_distinct(delta$MPA[delta$Species == sp]),
            size_refit_ar1_slope(correlated$model, sp))))
      }
      slopes[[paste0(tag, "|CR2")]] <- bind_rows(lapply(levels(delta$Species), function(sp)
        cbind(scope, Resp = response, size = "paired_threshold_difference", Species = sp,
          model = "MPA_cluster_CR2", n_obs = sum(delta$Species == sp),
          n_mpa = n_distinct(delta$MPA[delta$Species == sp]), size_refit_cluster_slope(fitted$model, sp))))
    }
  }
  slope_table <- bind_rows(slopes) %>% group_by(scope, Resp, size, model) %>%
    mutate(q = p.adjust(p, "BH")) %>% ungroup() %>%
    mutate(inference = "Conditional on annual reconstructed measurements; shared references and full size/count/detection uncertainty unresolved")
  write(slope_table, "table_s_size_refit_temporal_slopes")
  write(bind_rows(attempts), "table_s_size_refit_temporal_candidates")
  write(bind_rows(diagnostics), "table_s_size_refit_temporal_diagnostics")
  write(bind_rows(ar1_audits), "table_s_size_refit_ar1_candidates")
  write(bind_rows(annual), "table_s_size_refit_temporal_annual")
  write(bind_rows(predictions), "table_s_size_refit_temporal_predictions")
  write(bind_rows(shared), "table_s_size_refit_shared_reference_slopes")
  paths <- size_refit_pathways(predators, resource, primary)
  write(paths$summary, "table_s_size_refit_pathways")
  write(paths$effects, "table_s_size_refit_pathway_effects")
  write(paths$annual, "table_s_size_refit_pathway_annual")
  write(paths$omissions, "table_s_size_refit_pathway_omissions")
  write(data.frame(path = inputs, md5 = unname(tools::md5sum(inputs))), "table_s_size_refit_inputs")
  packages <- c("lme4", "lmerTest", "glmmTMB", "clubSandwich", "metafor", "sandwich", "dplyr", "tidyr")
  write(data.frame(package = packages, version = vapply(packages, function(p) as.character(packageVersion(p)), "")),
    "table_s_size_refit_versions")
  dir.create(file.path(root, "data/cache"), showWarnings = FALSE)
  saveRDS(models, file.path(root, "data/cache/predator_size_refits.rds"))
  stopifnot(identical(initial, tools::md5sum(protected)))
  source(file.path(root, "code/R/predator_size_refit_display.R"), local = TRUE)
  size_refit_display(root)
  write(data.frame(status = "completed", primary_changed = FALSE,
    uncertainty = "Conditional on reconstructed annual measurements; no full count/size/detection propagation; shared references remain"),
    "table_s_size_refit_status")
  status_recorded <- TRUE
  print(filter(slope_table, Species %in% c("P. interruptus", "S. pulcher")))
  print(paths$summary)
  cat("Actual matched size refits complete; primary fits and Figure 4 unchanged.\n")
})()
