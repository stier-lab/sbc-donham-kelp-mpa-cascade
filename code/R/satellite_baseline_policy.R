# Satellite before/after eligibility. A failed baseline is excluded, not routed to CI.
# The year labels are the supplied, audited BA labels; time is calendar year.
satellite_baseline_screen <- function(dat, minimum_years = 3L, alpha = 0.05) {
  required <- c('CA_MPA_Name_Short', 'year', 'BA', 'lnDiff')
  stopifnot(all(required %in% names(dat)))
  do.call(rbind, lapply(split(dat, dat$CA_MPA_Name_Short), function(d) {
    b <- d[d$BA == 'Before' & is.finite(d$year) & is.finite(d$lnDiff), ]
    if (anyDuplicated(b$year)) stop('Duplicate satellite baseline years: ', d$CA_MPA_Name_Short[1])
    n <- nrow(b)
    slope <- p <- NA_real_
    if (n >= minimum_years) {
      fit <- lm(lnDiff ~ year, data = b)
      co <- summary(fit)$coefficients
      if ('year' %in% rownames(co)) {
        slope <- co['year', 'Estimate']; p <- co['year', 'Pr(>|t|)']
      }
    }
    eligible <- n >= minimum_years && is.finite(p) && p >= alpha
    reason <- if (n < minimum_years) 'Insufficient baseline years' else if (!is.finite(p))
      'Baseline slope not estimable' else if (p < alpha) 'Detectable baseline trend' else 'Passes baseline screen'
    data.frame(MPA = d$CA_MPA_Name_Short[1], baseline_years = n,
      first_baseline_year = if (n) min(b$year) else NA_integer_,
      last_baseline_year = if (n) max(b$year) else NA_integer_,
      slope_per_year = slope, baseline_p = p, baseline_eligible = eligible, reason = reason)
  }))
}
