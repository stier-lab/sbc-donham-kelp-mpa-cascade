# Best-available extension, deliberately separate from the exact-replication audit.
suppressPackageStartupMessages({library(data.table); library(nlme); library(ggplot2)})
source("code/replications/malakhoff2021/helpers.R")
out <- "outputs/malakhoff2021"
d <- fread(file.path(out, "baci_input.csv"))
reserve <- sprintf("a-k-%02d", c(2, 6, 9, 14))
support <- rbindlist(lapply(c(2004L, 2005L), function(first) {
  d[year <= 2019, .(before = sum(is.finite(biomass) & year < first),
                    after = sum(is.finite(biomass) & year >= first)),
    by = .(site, species)][, first_post := first]
}))
# Freeze membership using the original window, both taxa and both period boundaries.
membership <- support[, .(min_before = min(before), min_after = min(after)), by = site]
membership[, included := min_before >= 2 & min_after >= 2]
membership[, status := ifelse(site %in% reserve, "Reserve", "Control")]
membership[, reason := fifelse(included, "At least two before and after years for both species and boundaries",
                              "Insufficient original-window before/after size-supported biomass")]
fwrite(support, file.path(out, "common_support_counts.csv"))
fwrite(membership, file.path(out, "common_support_membership.csv"))
keep <- membership[included == TRUE, site]
x <- d[site %in% keep]
endpoints <- list()
for (sp in c("red", "purple")) {
  endpoints[[paste0(sp, "_biomass")]] <- x[species == sp, .(site, year, value = biomass)]
  endpoints[[paste0(sp, "_small_mass")]] <- x[species == sp, .(site, year, value = small_mass_corrected)]
}
wide <- dcast(x, site + year ~ species, value.var = "biomass")
endpoints$total_biomass <- wide[, .(site, year, value = red + purple)]
models <- list()
for (endpoint in names(endpoints)) for (last in c(2019L, 2022L, 2023L)) {
  for (first in c(2004L, 2005L)) for (form in c("site_calendar_AR1", "site_independent")) {
    key <- paste(endpoint, last, first, form, sep = "__")
    m <- fit_baci(endpoints[[endpoint]], endpoint, last, first, form)
    m$row$omitted_reserve <- "none"
    models[[key]] <- m
  }
}
# These are influence checks, not extra replicates or independent hypothesis tests.
for (endpoint in names(endpoints)) for (omit in intersect(keep, reserve)) {
  key <- paste(endpoint, "omit", omit, sep = "__")
  m <- fit_baci(endpoints[[endpoint]][site != omit], endpoint, 2022L, 2005L, "site_calendar_AR1")
  m$row$omitted_reserve <- omit
  models[[key]] <- m
}
results <- rbindlist(lapply(models, `[[`, "row"), fill = TRUE)
results[, `:=`(eligible_for_manuscript = FALSE,
               scope = "fixed common-support cohort; modified covariance; exploratory",
               units = "BACI contrast in log(wet biomass in g/m2 + 1)")]
fwrite(results, file.path(out, "common_support_results.csv"))
saveRDS(lapply(models, `[[`, "fit"), file.path(out, "common_support_models.rds"))
inputs <- rbindlist(lapply(names(endpoints), function(n) copy(endpoints[[n]])[, endpoint := n]))
fwrite(inputs, file.path(out, "common_support_inputs.csv"))
# Before-period slope contrast checks the parallel-trend assumption, not equivalence.
pretrends <- rbindlist(lapply(names(endpoints), function(n) {
  z <- as.data.frame(endpoints[[n]][year < 2005 & is.finite(value)])
  z$status <- factor(ifelse(z$site %in% reserve, "Reserve", "Control"), c("Control", "Reserve"))
  z$time <- z$year - 1991
  fit <- tryCatch(lme(log1p(value) ~ status * time, random = ~1 | site,
                     correlation = corAR1(form = ~year | site), data = z), error = function(e) e)
  if (inherits(fit, "error")) return(data.table(endpoint = n, diagnostic = conditionMessage(fit)))
  tt <- summary(fit)$tTable["statusReserve:time", ]
  data.table(endpoint = n, estimate = tt["Value"], se = tt["Std.Error"],
             lower = tt["Value"] - qt(.975, tt["DF"]) * tt["Std.Error"],
             upper = tt["Value"] + qt(.975, tt["DF"]) * tt["Std.Error"],
             p = tt["p-value"], diagnostic = "Exploratory before-period slope contrast")
}))
fwrite(pretrends, file.path(out, "common_support_pretrends.csv"))
shown <- results[omitted_reserve == "none" & first_post == 2005 &
                   formulation == "site_calendar_AR1" & !grepl("small", endpoint)]
shown[, endpoint_label := factor(endpoint, c("red_biomass", "purple_biomass", "total_biomass"),
                                 c("Red urchin", "Purple urchin", "Both species (sum)"))]
p <- ggplot(shown, aes(estimate, endpoint_label, color = factor(end_year))) +
  geom_vline(xintercept = 0, color = "grey60", linetype = "dashed") +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", width = .3,
                position = position_dodge(width = .6)) +
  geom_point(position = position_dodge(width = .6), size = 2) +
  scale_color_manual(values = c("2019" = "#333333", "2022" = "#0072B2", "2023" = "#D55E00"),
                     labels = c("Through 2019", "Through 2022", "Through 2023 (partial)")) +
  labs(x = "Before-after reserve-control contrast in log(biomass + 1)", y = NULL, color = NULL) +
  theme_classic(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(out, "common_support_comparison.png"), p, width = 180, height = 95, units = "mm", dpi = 300)
ggsave(file.path(out, "common_support_comparison.pdf"), p, width = 180, height = 95, units = "mm", device = cairo_pdf)
observed <- x[year <= 2022, .(site, year, species, biomass,
                             status = ifelse(site %in% reserve, "Reserve", "Control"))]
setorder(observed, site, species, year)
# Break trajectories at unsampled years rather than drawing across the 2020 gap.
observed[, segment := cumsum(c(TRUE, diff(year) != 1)), by = .(site, species)]
means <- observed[, .(biomass = mean(biomass, na.rm = TRUE), n_sites = sum(is.finite(biomass))),
                  by = .(year, species, status)]
setorder(means, species, status, year)
means[, segment := cumsum(c(TRUE, diff(year) != 1)), by = .(species, status)]
fwrite(means, file.path(out, "common_support_trajectory_summary.csv"))
raw <- ggplot(observed, aes(year, biomass, color = status)) +
  geom_vline(xintercept = 2004.5, linetype = "dashed", color = "grey60") +
  geom_line(aes(group = interaction(site, segment)), alpha = .25, linewidth = .35) +
  geom_point(alpha = .4, size = .7) +
  geom_line(data = means, aes(group = interaction(status, segment)), linewidth = .9) +
  facet_wrap(~species, scales = "free_y", labeller = as_labeller(c(red = "Red urchin", purple = "Purple urchin"))) +
  scale_color_manual(values = c(Reserve = "#0072B2", Control = "#D55E00")) +
  labs(x = "Year", y = expression("Urchin wet biomass (g " * m^{-2} * ")"), color = NULL) +
  theme_classic(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(out, "common_support_trajectories.png"), raw, width = 180, height = 100, units = "mm", dpi = 300)
ggsave(file.path(out, "common_support_trajectories.pdf"), raw, width = 180, height = 100, units = "mm", device = cairo_pdf)
pdf(file.path(out, "common_support_diagnostics.pdf"), width = 9, height = 6)
for (key in names(models)) {
  m <- models[[key]]
  if (is.null(m$fit)) next
  par(mfrow = c(1, 2), oma = c(0, 0, 3, 0))
  plot(fitted(m$fit), residuals(m$fit, type = "normalized"), xlab = "Fitted log(biomass + 1)", ylab = "Normalized residual")
  abline(h = 0, col = "grey60")
  qqnorm(residuals(m$fit, type = "normalized")); qqline(residuals(m$fit, type = "normalized"))
  mtext(key, outer = TRUE, cex = .7)
}
dev.off()
capture.output(sessionInfo(), file = file.path(out, "common_support_sessionInfo.txt"))
print(membership)
print(shown[, .(endpoint, end_year, n, n_sites, n_reserve_sites, estimate, lower, upper, p, phi, residual_lag1)])
