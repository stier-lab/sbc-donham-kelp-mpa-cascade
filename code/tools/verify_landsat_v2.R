# Reconstruct the corrected masks independently from the cached NetCDF.
# Run from the analysis root; the external delivery remains locally held.
local({
  library(ncdf4)
  library(readxl)
  base <- Sys.getenv("LANDSAT_V2_DIR", "data/external/tom_bell_2026-10-21")
  ncpath <- Sys.getenv("LANDSAT_NETCDF", "local_archive/metadata_audit_2026-10-16/CAkelpCanopyEnv_2021.nc")
  out <- "outputs/landsat_v2_verification"
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  nc <- nc_open(ncpath)
  on.exit(nc_close(nc))
  key <- function(lat, lon) paste(round(lat, 6), round(lon, 6), sep = "|")
  keys <- key(ncvar_get(nc, "lat"), ncvar_get(nc, "lon"))
  stopifnot(!anyDuplicated(keys))
  qyear <- ncvar_get(nc, "year")
  years <- 1984:2021
  stopifnot(identical(as.integer(sort(unique(qyear))), years), all(table(qyear) == 4))
  masks <- read_excel(file.path(base, "MPA_Ref_Locations_v2.xlsx"),
                      sheet = "Pixel Locations", col_names = FALSE, .name_repair = "minimal")
  groups <- lapply(seq(2, ncol(masks), 2), function(j) {
    lat <- suppressWarnings(as.numeric(masks[[j]][-(1:2)]))
    lon <- suppressWarnings(as.numeric(masks[[j+1]][-(1:2)]))
    stopifnot(identical(is.finite(lat), is.finite(lon)))
    idx <- match(key(lat[is.finite(lat)], lon[is.finite(lon)]), keys)
    stopifnot(!anyNA(idx), !anyDuplicated(idx))
    list(site = as.integer(masks[[j]][1]), group = masks[[j]][2], idx = idx)
  })
  stopifnot(length(groups) == 72)
  q <- matrix(NA_real_, 72, length(qyear))
  for (t in seq_along(qyear)) {
    values <- ncvar_get(nc, "biomass", start = c(1,t), count = c(-1,1))
    for (j in seq_along(groups)) {
      x <- values[groups[[j]]$idx]
      if (length(x) && mean(is.na(x)) <= .25) q[j,t] <- sum(x, na.rm = TRUE)
    }
  }
  annual <- sapply(years, function(y) apply(q[, qyear == y, drop = FALSE], 1,
    function(x) if (sum(is.finite(x)) >= 3) max(x, na.rm = TRUE) else NA_real_))
  sheets <- c(MPA = "MPAs", Ref_1 = "MPAC1", Ref_2 = "MPAC2")
  delivered <- lapply(sheets, function(s) as.matrix(read_excel(
    file.path(base, "MPA_Runs_new.xlsx"), sheet = s, col_names = FALSE,
    col_types = "numeric", na = c("", "NaN"), .name_repair = "minimal")))
  stopifnot(all(vapply(delivered, function(x) identical(dim(x), c(24L,38L)), logical(1))))
  records <- lapply(seq_along(groups), function(j) {
    g <- groups[[j]]; a <- annual[j,]; b <- delivered[[g$group]][g$site,]
    finite <- is.finite(a) & is.finite(b)
    empty <- !length(g$idx)
    same_missing <- all(is.na(a) == is.na(b))
    delta <- if (any(finite)) max(abs(a[finite] - b[finite])) else NA_real_
    pass <- if (empty) all(is.na(b) | b == 0) else same_missing && is.finite(delta) && delta < 1e-6
    data.frame(site = g$site, group = g$group, n_pixels = length(g$idx),
      empty = empty, matching_missingness = same_missing, finite_years = sum(finite),
      max_abs_difference_kg = delta, pass = pass)
  })
  checks <- do.call(rbind, records)
  write.csv(checks, file.path(out, "mask_array_checks.csv"), row.names = FALSE)
  # Compare the exact historical site-centroid and replicate keys, not values alone.
  sites <- read.csv("data/external/tom_bell_2026-10-20/03_Site_List.csv", check.names = FALSE)
  raw <- read.csv("data/LANDSAT/MPA_Runs_new.csv", check.names = FALSE)
  input <- lapply(seq_len(nrow(raw)), function(i) {
    r <- raw[i,]
    s <- which(abs(sites$Lat - r$lat) < 1e-8 & abs(sites$Lon - r$lon) < 1e-8)
    stopifnot(length(s) == 1)
    group <- if (r$status == "mpa") "MPA" else paste0("Ref_", r$rep)
    a <- as.numeric(r[,as.character(years)]); b <- delivered[[group]][s,]
    ok <- all(is.na(a) == is.na(b)) && all(abs(a[is.finite(a)]-b[is.finite(a)]) < 1e-6)
    data.frame(site = s, MPA = r$MPA, status = r$status, rep = r$rep, pass = ok)
  })
  input <- do.call(rbind,input)
  write.csv(input, file.path(out,"input_checks.csv"), row.names = FALSE)
  registry <- read.csv("../Donham-Stier-CA-MPA-Data-2026/data/landsat_empty_reference_corrections.csv")
  slots <- merge(input, registry[c("MPA","status","rep")])
  slot_checks <- merge(slots, checks[checks$group == "Ref_2",], by = "site")
  stopifnot(nrow(slot_checks) == 6, all(slot_checks$empty))
  cat("Nonempty groups:",sum(!checks$empty),"; empty groups:",sum(checks$empty),
      "; maximum absolute difference (kg):",max(checks$max_abs_difference_kg,na.rm=TRUE),
      "; matched input rows:",sum(input$pass),"of",nrow(input),"\n")
  print(checks[!checks$pass,])
  stopifnot(all(checks$pass), nrow(input) == 60, all(input$pass))
  fingerprint <- function(path) list(path = path, bytes = file.info(path)$size,
    md5 = unname(tools::md5sum(path)), sha256 = digest::digest(file = path, algo = "sha256"))
  required <- c("data/harmonized/harmonized_landsat_rr.csv",
    "data/harmonized/harmonized_landsat_annual_index.csv",
    file.path(out,c("mask_array_checks.csv","input_checks.csv")))
  sources <- c(ncpath, "data/LANDSAT/MPA_Runs_new.csv",
    "data/external/tom_bell_2026-10-20/03_Site_List.csv",
    "../Donham-Stier-CA-MPA-Data-2026/data/landsat_empty_reference_corrections.csv",
    file.path(base,c("MPA_Final_2022.mat","MPA_Ref_Locations_v2.xlsx","MPA_Runs_new.xlsx",
      "reproduce_MPA_Runs_new.py","MPA_Analysis_2022.m","p_poly_dist.m","vline.m",
      "README.md","Mask_changes_v1_to_v2.csv")))
  manifest <- list(release = "landsat_v2_masks_2026_10_21", verified = TRUE,
    years = c(1984L,2021L), mpa_arrays = 24L, reference_arrays = 34L,
    empty_groups = 14L, input_rows = 60L, confirmed_removed_slots = 6L,
    max_abs_difference_kg = max(checks$max_abs_difference_kg,na.rm=TRUE),
    required_files = lapply(required, fingerprint), local_source_files = lapply(sources, fingerprint),
    reconstruction = "Independent R reconstruction from cached NetCDF; corrected masks reproduce delivered annual arrays including missingness. Empty masks are absent, not observed zeros.",
    limitation = "Some reference splits were recovered by the supplier against the delivered arrays. Reproduction does not establish unique historical clustering or resolve site eligibility.")
  jsonlite::write_json(manifest, "data/harmonized/landsat_v2_verification.json",
    pretty = TRUE, auto_unbox = TRUE, digits = NA)
  cat("PASS: NetCDF reconstruction, missingness, 60 delivered input rows and six empty-slot removals.\n")
})
