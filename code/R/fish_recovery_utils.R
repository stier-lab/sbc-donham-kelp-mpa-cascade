fish_recovery_select <- function(footprint, roster, years, start, endpoint, threshold) {
  years <- setdiff(years, start)
  periods <- if (endpoint == "post_trend") list(after = years[years > start]) else
    list(before = years[years < start], after = years[years > start])
  scores <- lapply(roster$site, function(s) {
    available <- unique(footprint$year[footprint$site == s])
    fractions <- vapply(periods, function(y) if (length(y)) mean(y %in% available) else 0, numeric(1))
    data.frame(site = s, min_coverage = min(fractions),
      selected = all(fractions >= threshold))
  })
  do.call(rbind, scores)
}

fish_recovery_eligible <- function(years, start, endpoint) {
  if (endpoint == "before_after") return(sum(years < start) >= 3 && sum(years > start) >= 5)
  y <- years[years > start]
  length(y) >= 8 && diff(range(y)) >= 8
}
