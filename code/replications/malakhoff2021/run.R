# Standalone reconstruction: does not overwrite the manuscript pipeline outputs.
suppressPackageStartupMessages({
  library(data.table)
  library(nlme)
  library(ggplot2)
})
source("code/replications/malakhoff2021/helpers.R")
data_root <- Sys.getenv("MALAKHOFF_DATA_ROOT",
                        "/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data")
out <- "outputs/malakhoff2021"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
counts_path <- file.path(data_root, "MBON/SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv")
sizes_path <- file.path(data_root, "ALL_sizefreq_2024.csv")
stopifnot(file.exists(counts_path), file.exists(sizes_path))
fwrite(data.table(path = c(counts_path, sizes_path),
                  md5 = unname(tools::md5sum(c(counts_path, sizes_path))),
                  source = c("MBON integrated KFM counts, 20231022 export",
                             "Combined size-frequency file, KFM-only subset")),
       file.path(out, "input_manifest.csv"))
counts <- fread(counts_path,
                select = c("data_source", "sample_method", "date", "site_id", "transect_id",
                           "replicate_id", "proj_taxon_id", "area", "count", "taxon_name"),
                na.strings = c("NA", "."))
counts <- counts[data_source == "kfm"]
counts[, year := as.integer(format(as.Date(date), "%Y"))]
counts[, species := fifelse(taxon_name %in% c("Mesocentrotus franciscanus", "Strongylocentrotus franciscanus"),
                             "red", fifelse(taxon_name == "Strongylocentrotus purpuratus", "purple", NA_character_))]
urchins <- counts[!is.na(species) & sample_method == "1mquad"]
stopifnot(all(is.finite(urchins$count) & urchins$count >= 0),
          all(is.finite(urchins$area) & urchins$area > 0))
keys <- c("site_id", "date", "transect_id", "replicate_id", "proj_taxon_id")
stopifnot(!anyDuplicated(urchins[, ..keys]))
annual <- urchins[, .(density = mean(count / area), count_rows = .N),
                  by = .(site = site_id, year, species)]
sizes <- fread(sizes_path, na.strings = c("NA", "."))
sizes <- sizes[campus == "KFM" & classcode %in%
                 c("Strongylocentrotus franciscanus", "Strongylocentrotus purpuratus")]
sizes[, species := fifelse(classcode == "Strongylocentrotus franciscanus", "red", "purple")]
stopifnot(all(is.finite(sizes$size) & sizes$size > 0),
          all(is.finite(sizes$count) & sizes$count > 0))
sizes[, wet_mass := urchin_mass(size, species)]
size_annual <- sizes[, .(n_sized = sum(count), mean_size = weighted.mean(size, count),
                         mean_mass = weighted.mean(wet_mass, count),
                         fraction_small_code = sum(count[size < 36]) / sum(count),
                         fraction_small_paper = sum(count[size < 35]) / sum(count),
                         small_mass_per_sampled_individual = sum(count[size < 35] * wet_mass[size < 35]) / sum(count)),
                     by = .(site, year, species)]
stopifnot(!anyDuplicated(annual[, .(site, year, species)]),
          !anyDuplicated(size_annual[, .(site, year, species)]))
annual <- merge(annual, size_annual, by = c("site", "year", "species"), all.x = TRUE)
annual[, biomass := density * mean_mass]
# Preserve observed count zeros in the corrected reconstruction, never invent sizes.
annual[density == 0, biomass := 0]
annual[, missing_positive_sizes := density > 0 & !is.finite(mean_mass)]
annual[, biomass_2017_mean_size := density * urchin_mass(
  fifelse(is.finite(mean_size), mean_size, 1), species)]
annual[!is.finite(mean_size), biomass_2017_mean_size := NA_real_]
annual[, biomass_code_era := fifelse(year >= 2017,
                                    biomass_2017_mean_size, density * mean_mass)]
annual[, small_count_fraction_code := fraction_small_code * biomass_code_era]
annual[, small_count_fraction_paper := fraction_small_paper * biomass]
annual[, small_mass_corrected := density * small_mass_per_sampled_individual]
annual[density == 0, `:=`(small_count_fraction_paper = 0, small_mass_corrected = 0)]
# The archived dplyr summarise drops the last (species) grouping, leaving
# size-class/year/site groups. Its denominator is therefore across species.
literal_small <- sizes[size < 36, .(n_small = sum(count)), by = .(site, year, species)]
literal_small[, literal_fraction := n_small / sum(n_small), by = .(site, year)]
annual <- merge(annual, literal_small[, .(site, year, species, literal_fraction)],
                by = c("site", "year", "species"), all.x = TRUE)
annual[, small_literal_grouping := literal_fraction * biomass_code_era]
fwrite(annual, file.path(out, "annual_urchins_all_sites.csv"))
cohort <- sprintf("a-k-%02d", c(1:11, 14:16))
baci <- annual[site %in% cohort & year >= 1991 & year <= 2023]
fwrite(baci, file.path(out, "baci_input.csv"))
membership <- baci[, .(first_count_year = min(year), last_count_year = max(year),
                       count_years = .N, biomass_years = sum(is.finite(biomass)),
                       missing_positive_sizes = sum(missing_positive_sizes)), by = .(site, species)]
fwrite(membership, file.path(out, "site_membership.csv"))
fwrite(baci[, .(n_sites = uniqueN(site), n_count_rows = sum(count_rows),
                missing_positive_sizes = sum(missing_positive_sizes),
                zero_density_sites = sum(density == 0), max_year = max(year)),
            by = .(species, year)], file.path(out, "coverage.csv"))
endpoints <- list()
for (sp in c("red", "purple")) {
  for (quantity in c("biomass", "biomass_code_era", "small_count_fraction_code",
                      "small_count_fraction_paper", "small_mass_corrected", "small_literal_grouping")) {
    x <- baci[species == sp, .(site, year, value = get(quantity))]
    endpoints[[paste(sp, quantity, sep = "_")]] <- x
  }
}
total <- dcast(baci, site + year ~ species, value.var = "biomass")
endpoints$total_biomass_sum <- total[, .(site, year, value = red + purple)]
models <- list()
for (name in names(endpoints)) {
  for (last in c(2019L, 2022L, 2023L)) {
    for (first in c(2004L, 2005L)) {
      # A scientifically interpretable sensitivity, not a literal execution of
      # the archived script's year/site singleton covariance grouping.
      offset <- if (name == "red_small_count_fraction_code") 1.5 else 1
      key <- paste(name, last, first, "site_calendar_AR1", sep = "__")
      message("Fitting ", key)
      models[[key]] <- fit_baci(endpoints[[name]], name, last, first,
                                "site_calendar_AR1", offset)
    }
  }
}
for (name in c("red_biomass_code_era", "purple_biomass_code_era")) {
  for (last in c(2019L, 2023L)) {
    key <- paste(name, last, "original_year_site", sep = "__")
    models[[key]] <- fit_baci(endpoints[[name]], name, last, 2005L, "original_year_site")
  }
}
results <- rbindlist(lapply(models, `[[`, "row"), fill = TRUE)
results[, `:=`(replication_status = "exploratory reconstruction; original CSVs unavailable",
               eligible_for_manuscript = FALSE,
               input_caveat = "missing positive-count sizes; 2020 counts absent; 2023 only three original sites")]
fwrite(results, file.path(out, "model_results.csv"))
saveRDS(lapply(models, `[[`, "fit"), file.path(out, "models.rds"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
# Anacapa is a separate, post-only comparison, not an extension of the BACI cohort.
anacapa <- annual[site %in% sprintf("a-k-%02d", 11:13) & is.finite(biomass)]
anacapa[, status := factor(ifelse(site == "a-k-11", "Control", "Reserve"))]
anacapa_results <- list()
for (sp in c("red", "purple")) for (last in c(2019L, 2023L)) {
  x <- anacapa[species == sp & year <= last]
  fit <- lm(biomass ~ year * status, data = x)
  tab <- as.data.table(coef(summary(fit)), keep.rownames = "term")
  tab[, `:=`(species = sp, end_year = last, first_year = min(x$year), n = nrow(x),
              scope = "ANCOVA reconstruction; size support starts in 1985, not 1981")]
  anacapa_results[[paste(sp, last)]] <- tab
}
fwrite(rbindlist(anacapa_results), file.path(out, "anacapa_ancova.csv"))

display <- baci[, .(value = mean(biomass, na.rm = TRUE),
                    se = sd(biomass, na.rm = TRUE) / sqrt(sum(is.finite(biomass))),
                    n_sites = sum(is.finite(biomass))),
                by = .(year, species, status = ifelse(site %in% sprintf("a-k-%02d", c(2, 6, 9, 14)),
                                                       "Reserve", "Control"))]
fwrite(display, file.path(out, "trajectory_summary.csv"))
p <- ggplot(display, aes(year, value, color = status, fill = status)) +
  geom_vline(xintercept = 2004, linetype = "dashed", color = "grey55") +
  geom_ribbon(data = display[is.finite(se)],
              aes(ymin = pmax(0, value - se), ymax = value + se), alpha = .15, color = NA) +
  geom_line(linewidth = .5) + geom_point(size = 1) +
  facet_wrap(~species, scales = "free_y", labeller = as_labeller(c(red = "Red urchin", purple = "Purple urchin"))) +
  scale_color_manual(values = c(Reserve = "#0072B2", Control = "#D55E00")) +
  scale_fill_manual(values = c(Reserve = "#0072B2", Control = "#D55E00")) +
  labs(x = "Year", y = expression("Urchin wet biomass (g " * m^{-2} * ")"), color = NULL, fill = NULL) +
  theme_classic(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(out, "urchin_trajectories.png"), p, width = 170, height = 85, units = "mm", dpi = 300)
ggsave(file.path(out, "urchin_trajectories.pdf"), p, width = 170, height = 85, units = "mm", device = cairo_pdf)

# Predator results were descriptive in the paper. Do not invent a predator BACI test.
predators <- counts[taxon_name %in% c("Panulirus interruptus", "Pycnopodia helianthoides") &
                      sample_method == "band" & site_id %in% cohort & year >= 1991,
                    .(site = site_id, year, taxon = taxon_name, count, area)]
fish_path <- file.path(data_root, "MBON/SBCMBON_kelp_forest_integrated_fish_20231022.csv")
fish <- fread(fish_path, select = c("data_source", "sample_method", "date", "site_id", "taxon_name", "count", "area"),
               na.strings = c("NA", "."))
fish <- fish[data_source == "kfm" & sample_method == "visualfish" & site_id %in% cohort &
               taxon_name %in% c("Semicossyphus pulcher", "Bodianus pulcher")]
fish[, year := as.integer(format(as.Date(date), "%Y"))]
predators <- rbindlist(list(predators, fish[year >= 1991, .(site = site_id, year,
                            taxon = "Bodianus pulcher", count, area)]))
predator_annual <- predators[, .(density = if (any(is.finite(count))) mean(count / area, na.rm = TRUE) else NA_real_,
                                  n_sample_rows = .N, n_missing_count = sum(!is.finite(count))),
                              by = .(site, year, taxon)]
fwrite(predator_annual, file.path(out, "predator_annual.csv"))
manifest <- fread(file.path(out, "input_manifest.csv"))
manifest <- rbind(manifest, data.table(path = fish_path, md5 = unname(tools::md5sum(fish_path)),
                                      source = "MBON integrated KFM visual fish, 20231022 export"))
fwrite(manifest, file.path(out, "input_manifest.csv"))
predator_annual[, status := ifelse(site %in% sprintf("a-k-%02d", c(2, 6, 9, 14)), "Reserve", "Control")]
pred <- ggplot(predator_annual[is.finite(density)], aes(year, density, group = site, color = status)) +
  geom_vline(xintercept = 2004, linetype = "dashed", color = "grey55") +
  geom_line(linewidth = .35, alpha = .5) + geom_point(size = .6, alpha = .7) +
  facet_wrap(~taxon, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c(Reserve = "#0072B2", Control = "#D55E00")) +
  labs(x = "Year", y = expression("Predator density (no. " * m^{-2} * ")"), color = NULL) +
  theme_classic(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(out, "predator_trajectories.png"), pred, width = 150, height = 170, units = "mm", dpi = 300)
ggsave(file.path(out, "predator_trajectories.pdf"), pred, width = 150, height = 170, units = "mm", device = cairo_pdf)
pdf(file.path(out, "model_diagnostics.pdf"), width = 8, height = 6)
for (key in names(models)) {
  m <- models[[key]]
  if (is.null(m$fit)) next
  par(mfrow = c(1, 2))
  plot(fitted(m$fit), residuals(m$fit, type = "normalized"), main = key,
       xlab = "Fitted log biomass", ylab = "Normalized residual")
  abline(h = 0, col = "grey60")
  qqnorm(residuals(m$fit, type = "normalized")); qqline(residuals(m$fit, type = "normalized"))
}
dev.off()
source("code/replications/malakhoff2021/common_support.R")
source("code/replications/malakhoff2021/method_bridge.R")
message("Reconstruction outputs written to ", out,
        "; exact original CSVs, kelp stipe-density input and 2025 restricted package remain unavailable.")
