# Separate diagnostic bridge; never writes manuscript inputs or primary results.
local({
  suppressPackageStartupMessages({library(data.table); library(nlme); library(ggplot2)})
  source("code/replications/malakhoff2021/helpers.R", local = TRUE)
  source("code/R/review_concern_utils.R", local = TRUE)
  out <- "outputs/malakhoff2021/bridge"
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  input_paths <- c("outputs/malakhoff2021/baci_input.csv",
                   "outputs/malakhoff2021/common_support_membership.csv",
                   "data/sumstats_final.csv", "data/harmonized/harmonized_response_ratios.csv",
                   "tables/table_s_design_matched_summary.csv",
                   "tables/table_s_design_pisco_common_calendar_summary.csv",
                   "data/harmonized/harmonized_raw_responses.csv",
                   "data/harmonized/reference_sites.csv", "tables/table_02_meta_analysis.csv")
  hashes <- unname(tools::md5sum(input_paths))
  fwrite(data.table(path = input_paths, md5 = hashes), file.path(out, "input_manifest.csv"))
  d <- fread(input_paths[1])
  membership <- fread(input_paths[2])
  sites <- membership[included == TRUE, site]
  reserves <- intersect(sites, sprintf("a-k-%02d", c(2, 6, 9, 14)))
  controls <- setdiff(sites, reserves)
  x <- d[site %in% sites & year >= 1991 & year <= 2022]
  calendar <- x[, .(n_rows = .N,
                    complete_biomass = all(is.finite(biomass)),
                    complete_density = all(is.finite(density)),
                    n_zero = sum(density == 0)), by = year]
  calendar[, included := n_rows == 2L * length(sites) & complete_biomass & complete_density]
  fwrite(calendar, file.path(out, "calendar_support.csv"))
  x <- x[year %in% calendar[included == TRUE, year]]
  # Match the allometry comparison's zero policy: a measured zero count has
  # zero biomass under either allometry, even without sampled sizes.
  x[density == 0, biomass_code_era := 0]
  x[, status := ifelse(site %in% reserves, "Reserve", "Control")]
  stopifnot(!anyDuplicated(x[, .(site, year, species)]),
            all(x[, .N, by = .(year, species)]$N == length(sites)))
  fwrite(x, file.path(out, "balanced_input.csv"))

  # The model grid changes one response/covariance/boundary/window at a time.
  endpoints <- list()
  for (sp in c("red", "purple")) {
    for (quantity in c("density", "biomass", "biomass_code_era")) {
      endpoints[[paste(sp, quantity, sep = "_")]] <- x[species == sp, .(site, year, value = get(quantity))]
    }
  }
  fits <- list()
  for (endpoint in names(endpoints)) for (last in c(2019L, 2022L)) {
    for (first in c(2004L, 2005L)) for (form in c("site_calendar_AR1", "site_independent")) {
      offsets <- if (grepl("_density$", endpoint)) c(.01, .1, 1) else 1
      for (offset in offsets) {
        key <- paste(endpoint, last, first, form, offset, sep = "__")
        m <- fit_baci(endpoints[[endpoint]], endpoint, last, first, form, offset)
        m$row$offset <- offset
        fits[[key]] <- m
      }
    }
  }
  results <- rbindlist(lapply(fits, `[[`, "row"))
  results[, `:=`(eligible_for_manuscript = FALSE,
                 scope = "identical balanced site-year records; modified BACI sensitivity")]
  results[, p_fdr := p.adjust(p, "BH"), by = .(end_year, first_post, formulation, offset)]
  fwrite(results, file.path(out, "balanced_model_results.csv"))
  saveRDS(lapply(fits, `[[`, "fit"), file.path(out, "balanced_models.rds"))

  # Missing sizes must not discard density-only evidence. The manuscript's
  # lone purple-density BACIPS series includes reserve site 14, which is absent
  # from the biomass-supported cohort above.
  density_all <- d[year %in% calendar[included == TRUE, year]]
  stopifnot(all(density_all[, .N, by = .(year, species)]$N == 14L),
            all(is.finite(density_all$density)))
  density_fits <- list()
  for (sp in c("red", "purple")) for (last in c(2019L, 2022L)) {
    for (first in c(2004L, 2005L)) for (form in c("site_calendar_AR1", "site_independent")) {
      for (offset in c(.01, .1, 1)) {
        key <- paste(sp, last, first, form, offset, sep = "__")
        m <- fit_baci(density_all[species == sp, .(site, year, value = density)],
                      paste0(sp, "_density"), last, first, form, offset)
        m$row$offset <- offset
        density_fits[[key]] <- m
      }
    }
  }
  density_results <- rbindlist(lapply(density_fits, `[[`, "row"))[,
    `:=`(cohort = "Paper14 density-only", eligible_for_manuscript = FALSE)]
  density_results <- rbindlist(list(density_results,
    copy(results[grepl("_density$", endpoint)])[, cohort := "Biomass-supported11"]), fill = TRUE)
  fwrite(density_results, file.path(out, "density_cohort_results.csv"))
  saveRDS(lapply(density_fits, `[[`, "fit"), file.path(out, "density_cohort_models.rds"))
  density_site_changes <- list()
  paper_reserves <- sprintf("a-k-%02d", c(2, 6, 9, 14))
  paper_controls <- setdiff(unique(density_all$site), paper_reserves)
  for (sp in c("red", "purple")) for (offset in c(.01, .1, 1)) {
    z <- density_all[species == sp]
    changes <- z[, .(before = mean(log(density[year < 2005] + offset)),
                      after = mean(log(density[year >= 2005] + offset))), by = site]
    for (reserve_id in paper_reserves) for (ref in c("All ten controls", "SBI mapped references")) {
      if (ref == "SBI mapped references" && reserve_id != "a-k-14") next
      selected <- if (ref == "All ten controls") paper_controls else c("a-k-15", "a-k-16")
      a <- changes[site == reserve_id]
      b <- changes[site %in% selected]
      density_site_changes[[paste(sp, offset, reserve_id, ref)]] <- data.table(species = sp, offset = offset,
        reserve = reserve_id, reference_choice = ref, before = a$before - mean(b$before),
        after = a$after - mean(b$after), step = (a$after - a$before) - mean(b$after - b$before),
        n_controls = nrow(b), eligible_for_manuscript = FALSE)
    }
  }
  fwrite(rbindlist(density_site_changes), file.path(out, "density_reserve_contrasts.csv"))

  # log(biomass) = log(density) + log(mean individual mass), on positive rows.
  # Enforce identical positive years for every site within each species; zeros
  # are kept in the model grid above, not assigned invented individual masses.
  positive_calendar <- x[, .(all_positive = all(density > 0 & mean_mass > 0 &
                                                is.finite(mean_mass))), by = .(species, year)]
  positive <- merge(x, positive_calendar[all_positive == TRUE, .(species, year)], by = c("species", "year"))
  fwrite(positive_calendar, file.path(out, "positive_calendar_support.csv"))
  fwrite(positive, file.path(out, "positive_decomposition_input.csv"))
  site_changes <- decomposition <- list()
  for (sp in c("red", "purple")) for (last in c(2019L, 2022L)) for (first in c(2004L, 2005L)) {
    z <- positive[species == sp & year <= last]
    for (quantity in c("density", "mean_mass", "biomass")) {
      s <- z[, .(before = mean(log(get(quantity)[year < first])),
                  after = mean(log(get(quantity)[year >= first])),
                  n_before = sum(year < first), n_after = sum(year >= first)), by = .(site, status)]
      s[, change := after - before]
      s[, `:=`(species = sp, end_year = last, first_post = first, quantity = quantity)]
      key <- paste(sp, last, first, quantity)
      site_changes[[key]] <- s
      a <- s[status == "Reserve", change]; b <- s[status == "Control", change]
      uncertainty <- t.test(a, b)
      decomposition[[key]] <- data.table(species = sp, end_year = last, first_post = first,
        quantity = quantity, estimate = mean(a) - mean(b), lower = uncertainty$conf.int[1],
        upper = uncertainty$conf.int[2], df = unname(uncertainty$parameter),
        n_reserves = length(a), n_controls = length(b),
        n_before_years = unique(s$n_before), n_after_years = unique(s$n_after),
        uncertainty = "exploratory Welch interval among site-level changes; only three reserves",
        eligible_for_manuscript = FALSE)
    }
  }
  site_changes <- rbindlist(site_changes)
  decomposition <- rbindlist(decomposition)
  fwrite(site_changes, file.path(out, "size_density_site_changes.csv"))
  fwrite(decomposition, file.path(out, "size_density_decomposition.csv"))
  identity <- dcast(decomposition, species + end_year + first_post ~ quantity, value.var = "estimate")
  identity[, error := biomass - density - mean_mass]
  stopifnot(max(abs(identity$error)) < 1e-10)
  fwrite(identity, file.path(out, "decomposition_identity.csv"))

  # A reference choice experiment, NOT validation of ecological matching.
  # Gull is deliberately absent: its mapped KFM reference is disputed.
  pairing <- data.table(reserve = c("a-k-02", "a-k-09"),
                        reference = c("a-k-01", "a-k-08"),
                        label = c("Harris Point", "Scorpion"),
                        basis = "declared legacy mapping; not a new ecological matching claim")
  fwrite(pairing, file.path(out, "reference_experiment_pairs.csv"))
  ref_history <- list()
  for (sp in c("red", "purple")) for (quantity in c("density", "biomass")) {
    z <- positive[species == sp]
    for (yr in sort(unique(z$year))) {
      a <- z[year == yr & site %in% pairing$reserve, get(quantity)]
      all_b <- z[year == yr & site %in% controls, get(quantity)]
      b <- z[year == yr & site %in% pairing$reference, get(quantity)]
      stopifnot(length(a) == 2L, length(b) == 2L, length(all_b) == length(controls))
      pairs <- merge(pairing, z[year == yr, .(reserve = site, a = get(quantity))], by = "reserve")
      pairs <- merge(pairs, z[year == yr, .(reference = site, b = get(quantity))], by = "reference")
      ref_history[[paste(sp, quantity, yr)]] <- data.table(species = sp, quantity = quantity, year = yr,
        reference_choice = c("All eight controls: mean logs", "Two mapped controls: mean logs",
                             "Two mapped controls: paired logs", "Two mapped controls: log arithmetic means"),
        value = c(mean(log(a)) - mean(log(all_b)), mean(log(a)) - mean(log(b)),
                  mean(log(pairs$a / pairs$b)), log(mean(a) / mean(b))))
    }
  }
  ref_history <- rbindlist(ref_history)
  fwrite(ref_history, file.path(out, "reference_annual_contrasts.csv"))
  ref_summary <- rbindlist(lapply(c(2019L, 2022L), function(last) {
    rbindlist(lapply(c(2004L, 2005L), function(first) {
      ref_history[year <= last, .(before = mean(value[year < first]), after = mean(value[year >= first]),
                                  step = mean(value[year >= first]) - mean(value[year < first]),
                                  n_before = sum(year < first), n_after = sum(year >= first)),
                  by = .(species, quantity, reference_choice)][, `:=`(end_year = last, first_post = first)]
    }))
  }))
  ref_summary[, `:=`(n_reserves = 2L, eligible_for_manuscript = FALSE,
                     scope = "descriptive reference/aggregation experiment; no network inference")]
  fwrite(ref_summary, file.path(out, "reference_comparison.csv"))

  # Piecewise eras distinguish a period average from changes within that period.
  eras <- list()
  era_fits <- list()
  for (sp in c("red", "purple")) for (quantity in c("density", "biomass")) {
    z <- copy(x[species == sp])
    offset <- if (quantity == "density") .01 else 1
    z[, response := log(get(quantity) + offset)]
    z[, status := factor(status, c("Control", "Reserve"))]
    z[, era := factor(fcase(year <= 2004, "Before", year <= 2013, "2005-2013",
                           year <= 2016, "2014-2016", year <= 2019, "2017-2019", default = "2021-2022"),
                       c("Before", "2005-2013", "2014-2016", "2017-2019", "2021-2022"))]
    warnings <- character()
    fit <- tryCatch(withCallingHandlers(lme(response ~ status * era, random = ~1 | site,
      correlation = corAR1(form = ~year | site), data = z,
      control = lmeControl(msMaxIter = 200)), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
      }), error = function(e) e)
    if (inherits(fit, "error")) stop("Era model failed: ", conditionMessage(fit))
    tt <- summary(fit)$tTable
    terms <- grep("statusReserve:era", rownames(tt), value = TRUE)
    tab <- as.data.table(tt[terms, , drop = FALSE], keep.rownames = "term")
    setnames(tab, c("Value", "Std.Error", "DF", "t-value", "p-value"), c("estimate", "se", "df", "t", "p"))
    tab[, `:=`(species = sp, quantity = quantity, offset = offset,
                lower = estimate - qt(.975, df) * se, upper = estimate + qt(.975, df) * se,
                warnings = paste(unique(warnings), collapse = "; "), eligible_for_manuscript = FALSE)]
    eras[[paste(sp, quantity)]] <- tab
    era_fits[[paste(sp, quantity)]] <- fit
  }
  era_results <- rbindlist(eras)
  era_results[, p_fdr := p.adjust(p, "BH")]
  fwrite(era_results, file.path(out, "era_model_results.csv"))
  saveRDS(era_fits, file.path(out, "era_models.rds"))

  # Inventory primary estimands without conflating sources with independent MPAs.
  primary <- fread("data/sumstats_final.csv")
  inventory <- primary[Taxa %in% c("S. purpuratus", "M. franciscanus"),
    .(n_estimates = .N, n_mpa = uniqueN(MPA),
      MPA_source_keys = paste(paste(MPA, Source, sep = "|"), collapse = ";")),
    by = .(Taxa, Resp, AnalysisType, Model)]
  fwrite(inventory, file.path(out, "primary_estimand_inventory.csv"))
  rr <- fread("data/harmonized/harmonized_response_ratios.csv")
  rr[, Taxa := fcase(y == "Strongylocentrotus purpuratus", "S. purpuratus",
                     y == "Mesocentrotus franciscanus", "M. franciscanus", default = NA_character_)]
  keys <- primary[Taxa %in% c("S. purpuratus", "M. franciscanus"),
    .(CA_MPA_Name_Short = MPA, Taxa, source = Source, resp = Resp, AnalysisType, Model,
      primary_estimate = Mean, primary_se = SE)]
  annual <- merge(rr, keys, by = c("CA_MPA_Name_Short", "Taxa", "source", "resp"))
  stopifnot(!anyDuplicated(annual[, .(CA_MPA_Name_Short, Taxa, source, resp, year)]))
  annual <- annual[is.finite(lnDiff)]
  primary_audit <- annual[, {
    lm_fit <- lm(lnDiff ~ time)
    before <- if (any(BA == "Before")) mean(lnDiff[BA == "Before"]) else NA_real_
    after <- if (any(BA == "After")) mean(lnDiff[BA == "After"]) else NA_real_
    estimated <- if (AnalysisType[1] == "CI" && Model[1] == "Linear")
      unname(coef(lm_fit)["time"] * 11) else if (AnalysisType[1] == "CI" && Model[1] == "Mean")
        mean(lnDiff) else NA_real_
    list(primary_estimate = primary_estimate[1], recomputed_primary = estimated,
      reconstruction_error = estimated - primary_estimate[1],
      n_before = sum(BA == "Before"), n_after = sum(BA == "After"),
      before_level = before, after_level = after, all_year_level = mean(lnDiff),
      descriptive_step = after - before, linear_11_model_intervals = unname(coef(lm_fit)["time"] * 11),
      linear_11_calendar_years = unname(coef(lm(lnDiff ~ year))["year"] * 11),
      extrapolates_model_endpoint = 11 > max(time),
      primary_quantity = if (AnalysisType[1] == "CI" && Model[1] == "Linear") "Fitted change over 11 model-time intervals" else
        if (AnalysisType[1] == "CI") "Mean log contrast over retained input years" else "Designation-associated fitted change")
  }, by = .(CA_MPA_Name_Short, Taxa, source, resp, AnalysisType, Model)]
  stopifnot(all(abs(primary_audit[AnalysisType == "CI", reconstruction_error]) < 1e-8))
  fwrite(primary_audit, file.path(out, "primary_estimate_reconstruction.csv"))
  fwrite(annual, file.path(out, "primary_annual_input.csv"))
  ci_linear <- keys[AnalysisType == "CI" & Model == "Linear"]
  trend_results <- trend_attempts <- list()
  for (i in seq_len(nrow(ci_linear))) {
    k <- ci_linear[i]
    z <- annual[CA_MPA_Name_Short == k$CA_MPA_Name_Short & Taxa == k$Taxa &
                  source == k$source & resp == k$resp]
    setorder(z, year)
    for (correlation in c("IID", "AR1")) {
      result <- review_fit_gls(as.data.frame(z), lnDiff ~ time, correlation)
      if (is.null(result$model)) stop("Primary-trend audit failed: ", result$error)
      contrast <- c("(Intercept)" = 0, time = 11)
      values <- review_gls_contrast(result$model, contrast)
      trend_results[[paste(i, correlation)]] <- cbind(k, correlation = correlation, values,
        review_gls_fit_audit(result$model), result$selection,
        scope = "conditional uncertainty for already-selected linear fits; not selection-adjusted",
        eligible_for_manuscript = FALSE)
      trend_attempts[[paste(i, correlation)]] <- cbind(k, correlation = correlation, result$attempts)
    }
  }
  trend_results <- rbindlist(trend_results)
  fwrite(trend_results, file.path(out, "primary_linear_covariance_check.csv"))
  fwrite(rbindlist(trend_attempts), file.path(out, "primary_linear_fit_attempts.csv"))
  # Ask whether changing uncertainty in these selected urchin fits changes the
  # existing joint synthesis; other taxa/model classes are deliberately fixed.
  if (!requireNamespace("metafor", quietly = TRUE)) stop("metafor is required for the joint-model sensitivity")
  meta_inputs <- meta_results <- meta_models <- list()
  replacements <- trend_results[correlation == "AR1"]
  for (mode in c("Current primary", "Urchin linear AR1 substitution", "Urchin linear variance inflation only")) {
    z <- copy(primary)
    z[, updated := FALSE]
    for (i in seq_len(nrow(replacements))) {
      a <- replacements[i]
      idx <- which(z$MPA == a$CA_MPA_Name_Short & z$Taxa == a$Taxa & z$Source == a$source & z$Resp == a$resp)
      stopifnot(length(idx) == 1L)
      if (mode == "Urchin linear AR1 substitution") {
        z[idx, `:=`(Mean = a$estimate, SE = a$se, updated = TRUE)]
      } else if (mode == "Urchin linear variance inflation only") {
        z[idx, `:=`(SE = max(SE, a$se), updated = TRUE)]
      }
    }
    z[, `:=`(vi = SE^2, scenario = mode)]
    meta_inputs[[mode]] <- z
    for (response in c("Bio", "Den")) {
      zz <- as.data.frame(z[Resp == response])
      zz$Taxa <- factor(zz$Taxa)
      warnings <- character()
      fit <- withCallingHandlers(metafor::rma.mv(yi = Mean, V = vi, mods = ~ Taxa - 1,
        random = list(~1 | MPA, ~1 | Source), data = zz, method = "REML", test = "t"),
        warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
      tab <- data.table(Taxa = sub("^Taxa", "", rownames(fit$b)), response = response,
        estimate = as.numeric(fit$b), se = fit$se, lower = fit$ci.lb, upper = fit$ci.ub, p = fit$pval,
        scenario = mode, total_k = nrow(zz), updated_k = sum(zz$updated),
        warnings = paste(unique(warnings), collapse = "; "), eligible_for_manuscript = FALSE,
        scope = "partial covariance substitution; mixed estimands and selected-model inference unchanged")
      meta_results[[paste(mode, response)]] <- tab
      meta_models[[paste(mode, response)]] <- fit
    }
  }
  meta_results <- rbindlist(meta_results)
  meta_results[, p_fdr := p.adjust(p, "BH"), by = scenario]
  baseline <- merge(meta_results[scenario == "Current primary"],
    fread("tables/table_02_meta_analysis.csv")[, .(Taxa, response = ifelse(Response == "Biomass", "Bio", "Den"),
                                                   recorded = Estimate)], by = c("Taxa", "response"))
  stopifnot(nrow(baseline) == 9L, max(abs(baseline$estimate - baseline$recorded)) < 1e-6)
  fwrite(rbindlist(meta_inputs), file.path(out, "joint_sensitivity_inputs.csv"))
  fwrite(meta_results, file.path(out, "joint_sensitivity_results.csv"))
  saveRDS(meta_models, file.path(out, "joint_sensitivity_models.rds"))

  # Geography context with one programme and a completely observed calendar.
  # This separates the calendar/programme axis, not nonrandom reef selection.
  regional_input <- annual[source == "PISCO" & year >= 2013 & year <= 2022 & BA == "After"]
  geographic_coverage <- regional_input[, .(n_years = uniqueN(year),
    years = paste(sort(year), collapse = ","), Region = ifelse(Location[1] == "C", "Islands", "Mainland")),
    by = .(CA_MPA_Name_Short, Taxa, source, resp)]
  geographic_coverage[, included := n_years == 10L]
  fwrite(geographic_coverage, file.path(out, "geography_complete_calendar_coverage.csv"))
  regional_input <- merge(regional_input, geographic_coverage[included == TRUE,
    .(CA_MPA_Name_Short, Taxa, source, resp, Region)], by = c("CA_MPA_Name_Short", "Taxa", "source", "resp"))
  geo_pairs <- regional_input[, .(after_level = mean(lnDiff),
    change_over_9_calendar_years = unname(coef(lm(lnDiff ~ year))["year"] * 9),
    n_years = .N), by = .(MPA = CA_MPA_Name_Short, Taxa, Source = source, Resp = resp, Region)]
  registry <- fread("data/harmonized/reference_sites.csv")
  refs <- registry[included == TRUE, .(reference_ids = paste(sort(unique(reference_id)), collapse = ";"),
    reference_mapping_complete = all(identity_verified)), by = .(MPA, Taxa, Source, Resp)]
  geo_pairs <- merge(geo_pairs, refs, by = c("MPA", "Taxa", "Source", "Resp"), all.x = TRUE)
  stopifnot(!anyNA(geo_pairs$reference_ids))
  geo_summary <- list()
  for (tx in unique(geo_pairs$Taxa)) for (response in unique(geo_pairs$Resp)) {
    for (region in c("All", "Islands", "Mainland")) for (quantity in c("after_level", "change_over_9_calendar_years")) {
      z <- geo_pairs[Taxa == tx & Resp == response & (Region == region | region == "All")]
      if (!nrow(z)) next
      groups <- review_reference_components(z$MPA, z$reference_ids)
      a <- review_group_summary(z[[quantity]], groups, z$reference_mapping_complete)
      geo_summary[[paste(tx, response, region, quantity)]] <- cbind(Taxa = tx, Resp = response, Region = region,
        quantity = quantity, first_year = 2013L, last_year = 2022L, a,
        eligible_for_manuscript = FALSE, scope = "fixed calendar and programme; descriptive geographic context")
    }
  }
  fwrite(geo_pairs, file.path(out, "geography_complete_calendar_pairs.csv"))
  fwrite(rbindlist(geo_summary), file.path(out, "geography_complete_calendar_summary.csv"))
  matched <- fread("tables/table_s_design_matched_summary.csv")
  matched <- matched[correlation == "Annual means" & source_subset == "All programmes" & Region == "All"]
  fwrite(matched, file.path(out, "manuscript_matched_design_context.csv"))
  regional <- fread("tables/table_s_design_pisco_common_calendar_summary.csv")
  fwrite(regional, file.path(out, "manuscript_common_calendar_geography.csv"))

  # Plot measured site histories, with arithmetic group means, not coefficient forests.
  metrics <- c(density = "Density (individuals / m2)", mean_mass = "Mean individual wet mass (g)",
               biomass = "Wet biomass (g / m2)")
  display <- melt(x, id.vars = c("site", "year", "species", "status"),
                  measure.vars = names(metrics), variable.name = "quantity", value.name = "value")
  display[, metric := factor(metrics[as.character(quantity)], unname(metrics))]
  setorder(display, site, species, metric, year)
  display[, segment := cumsum(c(TRUE, diff(year) != 1)), by = .(site, species, metric)]
  means <- display[, .(value = mean(value, na.rm = TRUE)), by = .(year, species, status, metric)]
  setorder(means, species, status, metric, year)
  means[, segment := cumsum(c(TRUE, diff(year) != 1)), by = .(species, status, metric)]
  fwrite(display, file.path(out, "observed_metric_histories.csv"))
  for (sp in c("red", "purple")) {
    p <- ggplot(display[species == sp], aes(year, value, color = status)) +
      geom_vline(xintercept = 2004.5, linetype = "dashed", color = "grey60") +
      geom_line(aes(group = interaction(site, segment)), alpha = .25, linewidth = .3, na.rm = TRUE) +
      geom_point(size = .6, alpha = .35, na.rm = TRUE) +
      geom_line(data = means[species == sp], aes(group = interaction(status, segment)), linewidth = .8) +
      facet_wrap(~metric, scales = "free_y", nrow = 1) +
      scale_color_manual(values = c(Control = "#D55E00", Reserve = "#0072B2")) +
      labs(x = "Year", y = NULL, color = NULL) + theme_classic(base_size = 9) +
      theme(legend.position = "bottom", strip.background = element_blank(), plot.margin = margin(8, 10, 8, 8))
    ggsave(file.path(out, paste0(sp, "_size_density_histories.png")), p, width = 200, height = 85, units = "mm", dpi = 300)
    ggsave(file.path(out, paste0(sp, "_size_density_histories.pdf")), p, width = 200, height = 85, units = "mm", device = cairo_pdf)
  }
  pdf(file.path(out, "diagnostics.pdf"), width = 9, height = 6)
  all_models <- c(lapply(fits, `[[`, "fit"), lapply(density_fits, `[[`, "fit"), era_fits)
  fwrite(data.table(page = seq_along(all_models), model_key = names(all_models)),
         file.path(out, "diagnostic_index.csv"))
  model_audit <- rbindlist(lapply(names(all_models), function(key) {
    m <- all_models[[key]]
    if (is.null(m)) return(data.table(model_key = key, fit_status = "failed"))
    random_sd <- as.numeric(VarCorr(m)[1, "StdDev"])
    residual <- residuals(m, type = "normalized")
    cbind(data.table(model_key = key, fit_status = "fitted", random_sd = random_sd,
      random_variance_fraction = random_sd^2 / (random_sd^2 + m$sigma^2),
      near_zero_random_fraction = random_sd^2 / (random_sd^2 + m$sigma^2) < 1e-6,
      max_abs_normalized_residual = max(abs(residual)), n_abs_residual_gt3 = sum(abs(residual) > 3),
      eligible_for_manuscript = FALSE), review_gls_fit_audit(m))
  }), fill = TRUE)
  fwrite(model_audit, file.path(out, "model_audit.csv"))
  for (key in names(all_models)) {
    m <- all_models[[key]]
    if (is.null(m)) next
    par(mfrow = c(1, 2), oma = c(0, 0, 3, 0))
    residual <- residuals(m, type = "normalized")
    plot(fitted(m), residual, xlab = "Fitted transformed response", ylab = "Normalized residual")
    abline(h = 0, col = "grey60")
    qqnorm(residual); qqline(residual)
    mtext(key, outer = TRUE, cex = .7)
  }
  dev.off()
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
  stopifnot(identical(hashes, unname(tools::md5sum(input_paths))))
  print(decomposition[end_year == 2022 & first_post == 2005,
                      .(species, quantity, estimate, lower, upper, n_before_years, n_after_years)])
  print(results[end_year == 2022 & first_post == 2005 & formulation == "site_calendar_AR1" &
                  (offset == 1 | offset == .01), .(endpoint, offset, estimate, lower, upper)])
  message("Method bridge written to ", out, "; manuscript unchanged.")
})
