# Fixed-cohort descriptive contrasts; annual observations never supply MPA df.
nonfished_metrics <- function(d, first_after) {
  stopifnot(!anyDuplicated(d$year), all(is.finite(d$contrast)))
  pre <- d[d$year < first_after, ]; post <- d[d$year >= first_after, ]
  if (nrow(pre) < 6L || nrow(post) < 6L) return(NULL)
  model <- lm(contrast ~ year, data = pre)
  track <- function(x, y) {
    if (length(x) < 6L || sd(x) < 1e-12 || sd(y) < 1e-12) return(NA_real_)
    cor(x, y)
  }
  tracking <- function(z, detrend = FALSE) {
    if (!all(c("log_inside", "log_reference") %in% names(z))) return(NA_real_)
    x <- z$log_inside; y <- z$log_reference
    if (detrend) {
      x <- residuals(lm(log_inside ~ year, data = z))
      y <- residuals(lm(log_reference ~ year, data = z))
    }
    track(x, y)
  }
  residual <- residuals(model)
  consecutive <- which(diff(pre$year) == 1L)
  lag1 <- if (length(consecutive) >= 3L)
    track(residual[consecutive], residual[consecutive + 1L]) else NA_real_
  data.frame(n_before = nrow(pre), n_after = nrow(post),
             before = mean(pre$contrast), after = mean(post$contrast),
             change = mean(post$contrast) - mean(pre$contrast),
             pretrend = unname(coef(model)["year"]),
             pretrend_residual_lag1 = lag1,
             tracking_before = tracking(pre), tracking_before_detrended = tracking(pre, TRUE),
             tracking_all = tracking(d), tracking_all_detrended = tracking(d, TRUE))
}

nonfished_pair <- function(site, roster, taxon, years) {
  z <- site[site$proj_taxon_id == taxon & site$year %in% years &
              site$site_id %in% roster$site_id, , drop = FALSE]
  # Never let a partial reference roster masquerade as the fixed reference.
  complete <- vapply(split(z$site_id, z$year), function(s) setequal(s, roster$site_id), logical(1))
  z <- z[z$year %in% as.integer(names(complete)[complete]), ]
  z <- merge(z, roster[c("site_id", "status")], by = "site_id", all.x = TRUE)
  if (!nrow(z)) return(data.frame())
  status <- aggregate(density ~ year + status, z, mean)
  wide <- reshape(status, idvar = "year", timevar = "status", direction = "wide")
  names(wide) <- sub("density.mpa", "inside", names(wide), fixed = TRUE)
  names(wide) <- sub("density.reference", "reference", names(wide), fixed = TRUE)
  wide[order(wide$year), ]
}
