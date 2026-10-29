source(here::here("code/R/lobster_biomass_utils.R"))
stopifnot(abs(lobster_mass_g(60)-205.450635613253)<1e-9)
one <- lobster_pool_moments(data.frame(size=60,count=1))
repeated <- lobster_pool_moments(data.frame(size=60,count=7))
stopifnot(one$mass_variance==0,repeated$mass_variance==0,
          one$exact_mass==repeated$exact_mass,one$rounded_mass<one$exact_mass)
x <- read.csv(here::here("tables/table_s_lobster_reconciliation_annual.csv"))
s <- read.csv(here::here("tables/table_s_lobster_reconciliation_status.csv"))
stopifnot(s$status=="completed_secondary_only",!s$primary_changed,
  s$max_cached_primary_difference<1e-10,s$max_zero_variance_pool_error<1e-8,
  nrow(x)==s$overlap_rows,!anyDuplicated(x[c("CA_MPA_Name_Short","year","status")]),
  all(abs(x$secondary_biomass-x$seasonal_exact_biomass)<1e-10))
cat("PASS: full-precision masses, singleton/repeated pools, authoritative primary reconstruction and secondary arithmetic\n")
