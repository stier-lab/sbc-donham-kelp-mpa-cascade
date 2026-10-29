library(dplyr)
read <- function(stem) read.csv(here::here("tables", paste0(stem, ".csv")))
x <- read("table_s_catalina_omission")
stopifnot(nrow(x) == 9, !anyDuplicated(x[c("Taxa", "Resp")]),
  !any(x$direction_changed), !any(x$bh_changed),
  max(abs(x$q - p.adjust(x$p, "BH"))) < 1e-12)
for (species in c("sheephead", "lobster")) {
  x <- read(paste0("table_s_", species, "_size_temporal"))
  stopifnot(nrow(x) > 0, all(x$n_years > 0), all(x$first_year <= x$last_year),
    all(is.finite(x$slope_per_calendar_year[x$n_years >= 3])))
}
x <- read("table_s_lobster_size_pool_uncertainty")
stopifnot(nrow(x) > 0, !anyDuplicated(x[c("CA_MPA_Name_Short", "year", "status")]),
  all(x$n_sized > 0), all(x$fraction_low >= 0 & x$fraction_high <= 1),
  all(x$fraction_low <= x$fraction_high), all(x$mass_low <= x$mass_high),
  all(x$large_mass_low <= x$large_mass_high), all(x$large_mass <= x$mass + 1e-10),
  all(x$large_mass_high <= x$mass_high + 1e-10),
  all(x$fraction_low[x$fraction == 0] == 0), all(x$fraction_high[x$fraction == 1] == 1))
r <- read("table_s_lobster_biomass_reconciliation")
v <- read("table_s_lobster_biomass_validation")
stopifnot(nrow(r) == v$overlap_rows, v$status == "secondary_only",
  abs(median(abs(r$difference)) - v$median_absolute_difference) < 1e-10)
cat("PASS: Catalina omission, temporal coverage, conditional size-pool bounds and biomass reconciliation\n")
