source("code/R/40_landsat_reference_check.R")
z <- check_landsat_reference_inputs(normalizePath("."))
stopifnot(z$new_rr_rows == 718L, z$old_rr_rows == 719L)
rr <- read.csv("data/harmonized/harmonized_landsat_rr.csv")
annual <- read.csv("data/harmonized/harmonized_landsat_annual_index.csv")
wide <- tidyr::pivot_wider(dplyr::mutate(annual, value = index + log_offset),
  id_cols = c(CA_MPA_Name_Short, year), names_from = status, values_from = value)
key <- function(x) paste(x$CA_MPA_Name_Short, x$year)
i <- match(key(rr), key(wide))
stopifnot(!anyNA(i), max(abs(rr$mpa - wide$mpa[i])) < 1e-10,
          max(abs(rr$reference - wide$reference[i])) < 1e-10)
helper <- "../Donham-Stier-CA-MPA-Data-2026/code/R/landsat_reference_utils.R"
if (file.exists(helper)) {
  source(helper)
  registry <- read.csv("../Donham-Stier-CA-MPA-Data-2026/data/landsat_empty_reference_corrections.csv")
  d <- registry[c("MPA", "status", "rep")]
  d$X2000 <- 0
  occupied <- d[1, ]; occupied$rep <- 1; occupied$X2000 <- 0
  d <- rbind(d, occupied)
  out <- remove_empty_landsat_references(d, registry)
  stopifnot(nrow(out) == 1L, out$rep == 1L, out$X2000 == 0)
  d$X2000[1] <- 1
  stopifnot(inherits(try(remove_empty_landsat_references(d, registry), silent=TRUE), "try-error"))
}
cat("PASS: exact six-slot removal, observed zeros retained, changed slots rejected, annual/ratio concordance\n")
