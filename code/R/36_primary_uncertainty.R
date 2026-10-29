# Conditional uncertainty on every retained primary model; no selection or input changes.
local({
  root <- here::here()
  source(file.path(root, "code/R/primary_uncertainty_utils.R"), local = TRUE)
  primary <- read.csv(file.path(root, "data/sumstats_final.csv"))
  files <- c("data/sumstats_final.csv", "data/harmonized/harmonized_response_ratios.csv",
    "data/harmonized/harmonized_landsat_rr.csv", "code/R/08_effect_sizes.R",
    "code/R/02_pBACIPS_function.R", "code/R/primary_uncertainty_utils.R")
  hashes <- unname(tools::md5sum(file.path(root, files)))
  stopifnot(nrow(primary) == 145L, !anyDuplicated(uncertainty_key(primary)))
  annual <- uncertainty_annual(root, primary)
  definitions <- uncertainty_definitions(root)
  rows <- vector("list", nrow(primary))
  for (i in seq_len(nrow(primary))) {
    row <- primary[i, ]; warnings <- character()
    rebuilt <- withCallingHandlers(uncertainty_reconstruct(row, annual, definitions),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
    lag_auto <- max(1L, floor(4 * (row$N / 100)^(2 / 9)))
    scenarios <- c(lag1 = 1L, lag3 = 3L, automatic = lag_auto)
    rows[[i]] <- do.call(rbind, lapply(names(scenarios), function(scenario) {
      v <- uncertainty_hac(rebuilt$fit, rebuilt$data$year, scenarios[[scenario]])
      se_hac <- sqrt(drop(t(rebuilt$gradient) %*% v %*% rebuilt$gradient))
      cbind(row[c("MPA", "Taxa", "Source", "Resp", "Model", "AnalysisType", "N")],
        scenario, lag_years = scenarios[[scenario]], clock = rebuilt$clock,
        primary_estimate = row$Mean, reconstructed_estimate = rebuilt$estimate,
        primary_se = row$SE, reconstructed_se = rebuilt$iid_se, hac_se = se_hac,
        synthesis_se = max(row$SE, se_hac),
        se_ratio = max(row$SE, se_hac) / row$SE,
        variance_floor_applied = se_hac < row$SE,
        year_min = min(rebuilt$data$year), year_max = max(rebuilt$data$year),
        missing_years = length(setdiff(seq.int(min(rebuilt$data$year), max(rebuilt$data$year)), rebuilt$data$year)),
        reconstruction_warnings = paste(unique(warnings), collapse = "; "))
    }))
  }
  histories <- do.call(rbind, rows)
  baseline <- uncertainty_joint(primary, primary$SE, "retained_primary")
  published <- read.csv(file.path(root, "tables/table_02_meta_analysis.csv"))
  published$Resp <- ifelse(published$Response == "Biomass", "Bio", "Den")
  j <- merge(baseline, published, by = c("Taxa", "Resp"))
  stopifnot(nrow(j) == 9L, max(abs(j$estimate - j$Estimate)) < 1e-7,
    max(abs(j$se - j$SE)) < 1e-7, max(abs(j$q - j$pval_fdr)) < 1e-7)
  pooled <- do.call(rbind, c(list(baseline), lapply(c("lag1", "lag3", "automatic"), function(scenario) {
    x <- histories[histories$scenario == scenario, ]
    se <- x$synthesis_se[match(uncertainty_key(primary), uncertainty_key(x))]
    stopifnot(length(se) == nrow(primary), all(is.finite(se)), all(se >= primary$SE))
    uncertainty_joint(primary, se, paste0("calendar_HAC_", scenario))
  })))
  pooled <- merge(pooled, baseline[c("Taxa", "Resp", "estimate", "q", "ci_low", "ci_high")],
    by = c("Taxa", "Resp"), suffixes = c("", "_primary"))
  pooled$direction_changed <- sign(pooled$estimate) != sign(pooled$estimate_primary)
  pooled$bh_detection_changed <- (pooled$q < .05) != (pooled$q_primary < .05)
  pooled$ci_width_ratio <- (pooled$ci_high - pooled$ci_low) / (pooled$ci_high_primary - pooled$ci_low_primary)
  stopifnot(identical(hashes, unname(tools::md5sum(file.path(root, files)))))
  write.csv(histories, file.path(root, "tables/table_s_primary_hac_histories.csv"), row.names = FALSE)
  write.csv(pooled, file.path(root, "tables/table_s_primary_hac_pooled.csv"), row.names = FALSE)
  write.csv(data.frame(path = files, md5 = hashes), file.path(root, "tables/table_s_primary_hac_inputs.csv"), row.names = FALSE)
  cat("Exact primary uncertainty: all 145 contrasts and SEs reconstructed; joint baseline reproduces Table 2.\n")
  print(pooled[pooled$scenario == "calendar_HAC_automatic", c("Taxa", "Resp", "estimate", "ci_low", "ci_high", "q", "bh_detection_changed")])
})
