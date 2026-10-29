# Comparisons use complete site rosters and shared calendar years.
fish_siting_fit <- function(d, endpoint, start, lag = 2L) {
  stopifnot(!anyDuplicated(d$year), all(is.finite(d$contrast)))
  d <- d[order(d$year), ]; d <- d[d$year != start, ]
  if (endpoint == "before_after") {
    if (sum(d$year < start) < 3L || sum(d$year > start) < 5L) return(NULL)
    d$post <- as.numeric(d$year > start)
    fit <- lm(contrast ~ post, d); term <- "post"
  } else {
    d <- d[d$year > start, ]
    if (nrow(d) < 8L || diff(range(d$year)) < 8L) return(NULL)
    fit <- lm(contrast ~ year, d); term <- "year"
  }
  v <- uncertainty_hac(fit, d$year, lag)
  se <- sqrt(max(0, v[term, term])); est <- unname(coef(fit)[term])
  critical <- qt(.975, df.residual(fit))
  r <- residuals(fit); adjacent <- which(diff(d$year) == 1L)
  lag1 <- if (length(adjacent) >= 3L && sd(r[adjacent]) > 1e-12 &&
    sd(r[adjacent + 1L]) > 1e-12) cor(r[adjacent], r[adjacent + 1L]) else NA_real_
  list(result = data.frame(endpoint, hac_lag = lag, n_years = nrow(d),
    n_before = sum(d$year < start), n_after = sum(d$year > start),
    estimate = est, se, ci_low = est - critical * se, ci_high = est + critical * se,
    residual_lag1 = lag1, years = paste(d$year, collapse = ";")),
    predictions = data.frame(year = d$year, observed = d$contrast, fitted = fitted(fit)))
}

fish_siting_tracking <- function(year, inside, reference) {
  if (length(year) < 6L) return(NA_real_)
  a <- residuals(lm(inside ~ year)); b <- residuals(lm(reference ~ year))
  if (sd(a) < 1e-12 || sd(b) < 1e-12) return(NA_real_)
  cor(a, b)
}

fish_siting_complete_years <- function(site_years, roster) {
  if (!length(roster) || !nrow(site_years)) return(integer())
  years <- split(site_years$site, site_years$year)
  as.integer(names(years)[vapply(years, function(x) all(roster %in% x), logical(1))])
}
