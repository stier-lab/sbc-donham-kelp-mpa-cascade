# SBC-LTER biomass coefficients, PAIN row, version 3 (Santa Barbara Coastal LTER et al. 2021).
# Carapace length in mm; wet mass in g. Keep precision identical to primary processing.
lobster_mass_g <- function(length_mm) 0.001352821 * length_mm^2.913963113

lobster_pool_moments <- function(d) {
  stopifnot(all(is.finite(d$size)), all(d$size > 0),
            all(is.finite(d$count)), all(d$count > 0))
  w <- d$count / sum(d$count)
  mass <- lobster_mass_g(d$size)
  data.frame(n_sized=sum(d$count), exact_mass=sum(w*mass),
    rounded_mass=sum(w*0.00135*d$size^2.91),
    mass_variance=sum(w*(mass-sum(w*mass))^2))
}
