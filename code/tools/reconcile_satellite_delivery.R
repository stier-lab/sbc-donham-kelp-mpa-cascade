# Diagnostic reconstruction only; never replaces the delivered primary series.
# Rscript code/tools/reconcile_satellite_delivery.R NETCDF_PATH
local({
  library(ncdf4)
  library(readxl)
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) == 1L)
  out <- "data/external/source_reconciliation_2026-10-17"
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  d <- nc_open(args[1]); on.exit(nc_close(d))
  lat <- ncvar_get(d, "lat"); lon <- ncvar_get(d, "lon")
  year <- ncvar_get(d, "year")
  key <- function(a, b) paste(round(a, 6), round(b, 6), sep = "|")
  keys <- key(lat, lon); stopifnot(!anyDuplicated(keys))
  wb <- read_excel("local_archive/metadata_audit_2026-10-16/MPA_Ref_Locations.xlsx",
                   sheet = "Pixel Locations", col_names = FALSE, .name_repair = "minimal")
  sites <- read.csv("data/external/source_metadata_2026-10-16/satellite_site_crosswalk.csv")
  raw <- read.csv("data/LANDSAT/MPA_Runs_new.csv", check.names = FALSE)
  sets <- lapply(seq(2L, ncol(wb), 2L), function(j) {
    a <- suppressWarnings(as.numeric(wb[[j]][-(1:2)]))
    b <- suppressWarnings(as.numeric(wb[[j+1L]][-(1:2)]))
    ok <- is.finite(a) & is.finite(b)
    idx <- match(key(a[ok], b[ok]), keys)
    stopifnot(!anyNA(idx))
    list(site = as.integer(wb[[j]][1]), group = as.character(wb[[j]][2]), idx = idx)
  })
  sums <- matrix(NA_real_, length(sets), length(year))
  for (t in seq_along(year)) {
    v <- ncvar_get(d, "biomass", start = c(1L, t), count = c(-1L, 1L))
    for (j in seq_along(sets)) {
      x <- v[sets[[j]]$idx]
      if (any(is.finite(x))) sums[j,t] <- sum(x, na.rm = TRUE)
    }
  }
  years <- sort(unique(year))
  annual <- sapply(years, function(y) apply(sums[, year == y, drop = FALSE], 1,
    function(x) if (all(!is.finite(x))) NA_real_ else max(x, na.rm = TRUE)))
  rows <- list(); geometry <- list(); missing_years <- list()
  km <- function(a,b,c,e) {
    a <- a*pi/180; b <- b*pi/180; c <- c*pi/180; e <- e*pi/180
    6371*2*asin(pmin(1,sqrt(sin((a-c)/2)^2+cos(a)*cos(c)*sin((b-e)/2)^2)))
  }
  for (j in seq_along(sets)) {
    s <- sets[[j]]
    if (!length(s$idx)) next
    mpa <- sites$MPA[sites[[1L]] == s$site]
    if (length(mpa) != 1L || is.na(mpa)) next
    inside <- Filter(function(x) x$site == s$site && x$group == "MPA", sets)
    stopifnot(length(inside) == 1L)
    if (length(inside[[1]]$idx) && s$group != "MPA") {
      a <- mean(lat[inside[[1]]$idx]); b <- mean(lon[inside[[1]]$idx])
      distances <- km(lat[s$idx], lon[s$idx], a, b)
      geometry[[length(geometry)+1L]] <- data.frame(MPA = mpa, group = s$group,
        n_pixels = length(s$idx), centroid_to_centroid_km = km(mean(lat[s$idx]),mean(lon[s$idx]),a,b),
        min_pixel_km = min(distances), max_pixel_km = max(distances),
        shared_with_own_mpa = length(intersect(s$idx, inside[[1]]$idx)))
    }
    candidates <- which(raw$MPA == mpa & (raw$status == "mpa") == (s$group == "MPA"))
    for (i in candidates) {
      supplied <- as.numeric(raw[i,as.character(years)])
      reconstructed <- annual[j,]
      paired <- is.finite(supplied) & is.finite(reconstructed)
      missing_reconstruction <- !is.finite(reconstructed) & is.finite(supplied)
      missing_delivery <- is.finite(reconstructed) & !is.finite(supplied)
      if (any(paired) && max(abs(supplied[paired]-reconstructed[paired])) < 1e-6 &&
          any(missing_delivery)) {
        missing_years[[length(missing_years)+1L]] <- data.frame(MPA=mpa,group=s$group,
          delivered_rep=raw$rep[i],year=years[missing_delivery],
          reconstructed_biomass=reconstructed[missing_delivery],
          n_finite_quarters=vapply(years[missing_delivery],function(y)
            sum(is.finite(sums[j,year==y])),integer(1)))
      }
      zero_filled <- reconstructed; zero_filled[!is.finite(zero_filled)] <- 0
      rows[[length(rows)+1L]] <- data.frame(MPA = mpa, group = s$group, delivered_rep = raw$rep[i],
        paired_years = sum(paired), max_abs_finite_difference = if(any(paired)) max(abs(supplied[paired]-reconstructed[paired])) else NA_real_,
        reconstructed_missing_delivered_finite = sum(missing_reconstruction),
        of_those_delivered_zero = sum(supplied[missing_reconstruction] == 0),
        delivered_missing_reconstructed_finite = sum(missing_delivery),
        same_missingness = identical(is.finite(supplied),is.finite(reconstructed)),
        exact_after_all_missing_year_to_zero = all(is.finite(supplied)) && max(abs(supplied-zero_filled)) < 1e-6)
    }
  }
  result <- do.call(rbind, rows)
  write.csv(result, file.path(out,"satellite_missingness_diagnostics.csv"), row.names=FALSE)
  write.csv(do.call(rbind,missing_years),file.path(out,"satellite_missing_years.csv"),row.names=FALSE)
  geo <- do.call(rbind,geometry)
  write.csv(geo,file.path(out,"satellite_reference_centroid_geometry.csv"),row.names=FALSE)
  exact <- result[result$max_abs_finite_difference < 1e-6 & !is.na(result$max_abs_finite_difference),]
  cat("Finite-value exact candidates:",nrow(exact),"; complete matches:",sum(exact$same_missingness),
      "; exact after zero fill:",sum(exact$exact_after_all_missing_year_to_zero),"\n")
  cat("Reference centroid pairs outside 1-3 km:",sum(geo$centroid_to_centroid_km < 1 | geo$centroid_to_centroid_km > 3),"of",nrow(geo),"\n")
})
