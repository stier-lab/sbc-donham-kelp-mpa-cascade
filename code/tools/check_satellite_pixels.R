# Read-only check of the recovered 2022 pixel workbook and linked NetCDF.
# Rscript code/tools/check_satellite_pixels.R /tmp/CAkelpCanopyEnv_2021.nc
local({
  library(ncdf4)
  library(readxl)
  library(jsonlite)
  args <- commandArgs(trailingOnly = TRUE)
  stopifnot(length(args) == 1L)
  out <- "data/external/source_metadata_2026-10-16"
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  d <- nc_open(args[1])
  on.exit(nc_close(d))
  lat <- ncvar_get(d, "lat")
  lon <- ncvar_get(d, "lon")
  years <- ncvar_get(d, "year")
  quarters <- ncvar_get(d, "quarter")
  key <- function(a, b) paste(round(a, 6), round(b, 6), sep = "|")
  nc_keys <- key(lat, lon)
  stopifnot(!anyDuplicated(nc_keys))
  wb <- read_excel("local_archive/metadata_audit_2026-10-16/MPA_Ref_Locations.xlsx",
                   sheet = "Pixel Locations", col_names = FALSE, .name_repair = "minimal")
  sets <- lapply(seq(2L, ncol(wb), 2L), function(j) {
    a <- suppressWarnings(as.numeric(wb[[j]][-(1:2)]))
    b <- suppressWarnings(as.numeric(wb[[j + 1L]][-(1:2)]))
    finite <- is.finite(a) & is.finite(b)
    idx <- match(key(a[finite], b[finite]), nc_keys)
    list(site = as.integer(wb[[j]][1]), group = as.character(wb[[j]][2]),
         n_pixels = sum(finite), n_unmatched = sum(is.na(idx)), idx = idx[!is.na(idx)])
  })
  sites <- as.data.frame(read_excel("local_archive/metadata_audit_2026-10-16/MPA_Ref_Locations.xlsx",
                                    sheet = "SItes"))
  raw <- read.csv("data/LANDSAT/MPA_Runs_new.csv", check.names = FALSE)
  sites$MPA <- vapply(seq_len(nrow(sites)), function(i) {
    hit <- which(raw$status == "mpa" & abs(raw$lat - sites$Lat[i]) < 1e-8 &
                   abs(raw$lon - sites$Lon[i]) < 1e-8)
    if (length(hit) == 1L) raw$MPA[hit] else NA_character_
  }, character(1))
  write.csv(sites, file.path(out, "satellite_site_crosswalk.csv"), row.names = FALSE)
  # Read each quarter once. Compute spatial sums and means before annual maxima.
  summed <- averaged <- matrix(NA_real_, length(sets), length(years))
  for (t in seq_along(years)) {
    v <- ncvar_get(d, "biomass", start = c(1L, t), count = c(-1L, 1L))
    for (j in seq_along(sets)) {
      x <- v[sets[[j]]$idx]
      if (length(x) && any(is.finite(x))) {
        summed[j, t] <- sum(x, na.rm = TRUE)
        averaged[j, t] <- mean(x, na.rm = TRUE)
      }
    }
  }
  annual <- function(m) sapply(sort(unique(years)), function(y) {
    apply(m[, years == y, drop = FALSE], 1, function(x) {
      if (all(!is.finite(x))) NA_real_ else max(x, na.rm = TRUE)
    })
  })
  annual_sum <- annual(summed)
  annual_mean <- annual(averaged)
  raw_years <- names(raw)[grepl("^[0-9]{4}$", names(raw))]
  x <- as.matrix(raw[, raw_years])
  matches <- list()
  for (j in seq_along(sets)) {
    for (method in c("sum", "mean")) {
      v <- if (method == "sum") annual_sum[j, ] else annual_mean[j, ]
      v <- v[match(as.integer(raw_years), sort(unique(years)))]
      for (i in seq_len(nrow(raw))) {
        mpa <- sites$MPA[sites[[1L]] == sets[[j]]$site]
        if (length(mpa) != 1L || is.na(mpa) || raw$MPA[i] != mpa) next
        if ((sets[[j]]$group == "MPA") != (raw$status[i] == "mpa")) next
        ok <- is.finite(v) & is.finite(x[i, ])
        if (sum(ok) < 20L) next
        matches[[length(matches) + 1L]] <- data.frame(
          site_number = sets[[j]]$site, pixel_group = sets[[j]]$group, method = method,
          MPA = raw$MPA[i], status = raw$status[i], rep = raw$rep[i],
          paired_finite_years = sum(ok), same_missing = all(is.na(v) == is.na(x[i, ])),
          max_abs_difference = max(abs(v[ok] - x[i, ok])),
          mean_abs_difference = mean(abs(v[ok] - x[i, ok])))
      }
    }
  }
  comparisons <- do.call(rbind, matches)
  write.csv(comparisons, file.path(out, "satellite_pixel_comparison.csv"), row.names = FALSE)
  coverage <- do.call(rbind, lapply(sets, function(s) {
    i <- which(sites[[1L]] == s$site)
    a <- lat[s$idx] * pi / 180
    b <- lon[s$idx] * pi / 180
    a0 <- sites$Lat[i] * pi / 180
    b0 <- sites$Lon[i] * pi / 180
    km <- 6371 * 2 * asin(pmin(1, sqrt(sin((a - a0)/2)^2 + cos(a)*cos(a0)*sin((b - b0)/2)^2)))
    data.frame(site_number = s$site, group = s$group, n_pixels = s$n_pixels,
               n_unmatched = s$n_unmatched,
               min_km_from_supplied_mpa_centroid = if(length(km)) min(km) else NA_real_,
               median_km_from_supplied_mpa_centroid = if(length(km)) median(km) else NA_real_,
               max_km_from_supplied_mpa_centroid = if(length(km)) max(km) else NA_real_)
  }))
  write.csv(coverage, file.path(out, "satellite_pixel_coverage.csv"), row.names = FALSE)
  empty <- do.call(rbind, lapply(seq_along(sets), function(j) {
    s <- sets[[j]]
    if (s$n_pixels != 0L) return(NULL)
    mpa <- sites$MPA[sites[[1L]] == s$site]
    i <- which(raw$MPA == mpa & raw$status == "reference" & raw$rep == 2L)
    if (s$group != "Ref_2" || length(i) != 1L) return(NULL)
    data.frame(MPA = mpa, status = "reference", rep = 2L, site_number = s$site,
               pixel_group = s$group, n_pixels = 0L, all_delivered_years_zero = all(x[i, ] == 0))
  }))
  write.csv(empty, file.path(out, "satellite_empty_reference_rows.csv"), row.names = FALSE)
  overlap <- list()
  refs <- which(vapply(sets, function(s) s$group != "MPA" && length(s$idx) > 0L, logical(1)))
  for (a in seq_along(refs)) for (b in seq_along(refs)) {
    if (a >= b) next
    u <- sets[[refs[a]]]; v <- sets[[refs[b]]]
    n <- length(intersect(u$idx, v$idx))
    if (n > 0L) overlap[[length(overlap) + 1L]] <- data.frame(
      site_a = u$site, group_a = u$group, site_b = v$site, group_b = v$group,
      n_shared_pixels = n, fraction_a = n / length(u$idx), fraction_b = n / length(v$idx))
  }
  overlap <- if (length(overlap)) do.call(rbind, overlap) else data.frame(
    site_a=integer(), group_a=character(), site_b=integer(), group_b=character(),
    n_shared_pixels=integer(), fraction_a=double(), fraction_b=double())
  write.csv(overlap, file.path(out, "satellite_reference_overlap.csv"), row.names = FALSE)
  write_json(list(file = basename(args[1]), bytes = file.info(args[1])$size,
                  md5 = unname(tools::md5sum(args[1])),
                  global_attributes = ncatt_get(d, 0),
                  attributes = setNames(lapply(names(d$var), function(v) ncatt_get(d, v)), names(d$var)),
                  n_stations = length(lat), n_quarters = length(years),
                  year_range = range(years), quarter_values = sort(unique(quarters)),
                  coordinate_matching_precision = "6 decimal degrees",
                  exact_annual_sum_matches = sum(comparisons$method == "sum" &
                    comparisons$max_abs_difference < 1e-6 & comparisons$same_missing)),
             file.path(out, "satellite_netcdf_metadata.json"), pretty = TRUE, auto_unbox = TRUE)
  cat("Pixel sets:", nrow(coverage), "; unmatched pixels:", sum(coverage$n_unmatched), "\n")
  cat("Exact complete annual-vector matches:", sum(comparisons$method == "sum" &
       comparisons$max_abs_difference < 1e-6 & comparisons$same_missing), "\n")
  print(empty)
})
