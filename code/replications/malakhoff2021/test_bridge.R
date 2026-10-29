suppressPackageStartupMessages(library(data.table))
out <- "outputs/malakhoff2021/bridge"
x <- fread(file.path(out, "balanced_input.csv"))
calendar <- fread(file.path(out, "calendar_support.csv"))
stopifnot(uniqueN(x$site) == 11L, uniqueN(x$year) == 29L,
          !anyDuplicated(x[, .(site, year, species)]),
          all(x[, .N, by = .(year, species)]$N == 11L),
          all(is.finite(x$biomass)), all(is.finite(x$biomass_code_era)),
          all(x[density == 0, biomass_code_era] == 0),
          !any(c(1993L, 2003L, 2020L) %in% x$year),
          setequal(unique(x$year), calendar[included == TRUE, year]))
r <- fread(file.path(out, "balanced_model_results.csv"))
stopifnot(nrow(r) == 80L, all(r$status == "fitted"), !any(r$eligible_for_manuscript),
          all(r$n_sites == 11L), all(r$n_reserve_sites == 3L))
for (last in unique(r$end_year)) {
  stopifnot(all(r[end_year == last, n] == 11L * sum(calendar$included & calendar$year <= last)))
}
identity <- fread(file.path(out, "decomposition_identity.csv"))
stopifnot(max(abs(identity$error)) < 1e-10)
p <- fread(file.path(out, "positive_decomposition_input.csv"))
stopifnot(all(p$density > 0), all(p$mean_mass > 0),
          max(abs(log(p$biomass) - log(p$density) - log(p$mean_mass))) < 1e-10,
          all(p[, .N, by = .(species, year)]$N == 11L))
refs <- fread(file.path(out, "reference_comparison.csv"))
pair <- refs[reference_choice == "Two mapped controls: paired logs"]
pool <- refs[reference_choice == "Two mapped controls: mean logs"]
check <- merge(pair, pool, by = c("species", "quantity", "end_year", "first_post"))
stopifnot(max(abs(check$step.x - check$step.y)) < 1e-10,
          max(abs(check$before.x - check$before.y)) < 1e-10,
          max(abs(check$after.x - check$after.y)) < 1e-10)
pairing <- fread(file.path(out, "reference_experiment_pairs.csv"))
stopifnot(!"a-k-06" %in% pairing$reserve)
density <- fread(file.path(out, "density_cohort_results.csv"))
stopifnot(all(density$status == "fitted"),
          all(density[cohort == "Paper14 density-only", n_sites] == 14L),
          all(density[cohort == "Paper14 density-only", n_reserve_sites] == 4L),
          all(density[cohort == "Paper14 density-only", n] ==
                ifelse(density[cohort == "Paper14 density-only", end_year] == 2022, 14L * 29L, 14L * 27L)))
primary <- fread(file.path(out, "primary_estimate_reconstruction.csv"))
stopifnot(all(abs(primary[AnalysisType == "CI", reconstruction_error]) < 1e-8))
counts <- fread(file.path(out, "primary_estimand_inventory.csv"))
purple <- counts[Taxa == "S. purpuratus" & Resp == "Den"]
stopifnot(purple[AnalysisType == "CI" & Model == "Linear", n_estimates] == 9L,
          purple[AnalysisType == "CI" & Model == "Mean", n_estimates] == 3L,
          purple[AnalysisType == "pBACIPS" & Model == "Step", n_estimates] == 1L)
trend <- fread(file.path(out, "primary_linear_covariance_check.csv"))
stopifnot(!any(trend$eligible_for_manuscript), nrow(trend) == 46L,
          all(trend$ci_low <= trend$estimate & trend$ci_high >= trend$estimate))
joint <- fread(file.path(out, "joint_sensitivity_inputs.csv"))
baseline <- joint[scenario == "Current primary"]
for (scenario_id in setdiff(unique(joint$scenario), "Current primary")) {
  z <- merge(baseline, joint[scenario == scenario_id], by = c("MPA", "Taxa", "Source", "Resp"))
  stopifnot(nrow(z) == nrow(baseline), all(z[updated.y == FALSE, Mean.x == Mean.y]),
            all(z[updated.y == FALSE, SE.x == SE.y]))
  if (scenario_id == "Urchin linear variance inflation only")
    stopifnot(all(z$Mean.x == z$Mean.y), all(z$SE.y >= z$SE.x))
}
geo <- fread(file.path(out, "geography_complete_calendar_pairs.csv"))
stopifnot(all(geo$n_years == 10L), all(geo$Source == "PISCO"))
manifest <- fread(file.path(out, "input_manifest.csv"))
stopifnot(identical(manifest$md5, unname(tools::md5sum(manifest$path))))
audit <- fread(file.path(out, "model_audit.csv"))
pages <- fread(file.path(out, "diagnostic_index.csv"))
stopifnot(nrow(audit) == 132L, nrow(pages) == 132L,
          setequal(audit$model_key, pages$model_key), !any(audit$eligible_for_manuscript),
          !anyDuplicated(audit$model_key))
cat("PASS: common calendars, endpoint support, log identities, reference algebra, density-only cohort, primary reconstruction, covariance substitutions, geography and unchanged inputs.\n")
