source("code/R/primary_uncertainty_utils.R")
source("code/R/fish_siting_utils.R")
source("code/R/review_concern_utils.R")

sy <- data.frame(site = c("a", "b", "a"), year = c(2000, 2000, 2001))
stopifnot(identical(fish_siting_complete_years(sy, c("a", "b")), 2000L),
  !length(fish_siting_complete_years(sy, c("a", "b", "c"))))
d <- data.frame(year = 1999:2015, contrast = ifelse(1999:2015 > 2003, 2, 0))
f <- fish_siting_fit(d, "before_after", 2003)
stopifnot(abs(f$result$estimate - 2) < 1e-10, f$result$n_before == 4,
  !2003 %in% f$predictions$year,
  is.null(fish_siting_fit(d[d$year >= 2002, ], "before_after", 2003)))
d$contrast <- .2 * d$year
stopifnot(abs(fish_siting_fit(d, "post_trend", 2003)$result$estimate - .2) < 1e-10,
  is.na(fish_siting_tracking(1:5, 1:5, 1:5)),
  is.na(fish_siting_tracking(1:8, rep(0, 8), 1:8)))

# Verify covariance of the direct difference, rather than comparing two P values.
set.seed(47)
d$contrast <- rnorm(nrow(d))
e <- d; e$contrast <- e$contrast + .3 * e$year
delta <- d; delta$contrast <- e$contrast - d$contrast
stopifnot(abs(fish_siting_fit(delta, "post_trend", 2003)$result$estimate - .3) < 1e-10,
  fish_siting_fit(delta, "post_trend", 2003)$result$se < 1e-8)
stopifnot(length(unique(review_reference_components(c("a", "b", "c"), c("x", "x;y", "z")))) == 2)

if (file.exists("tables/table_s_fish_siting_status.csv")) {
  read <- function(x) read.csv(paste0("tables/table_s_fish_siting_", x, ".csv"))
  inputs <- read("inputs")
  stopifnot(identical(inputs$md5, unname(tools::md5sum(inputs$path))))
  guard <- read("protected")
  stopifnot(identical(guard$before_md5, guard$after_md5),
    identical(guard$after_md5, unname(tools::md5sum(guard$path))))
  coverage <- read("coverage"); annual <- read("annual"); fits <- read("fits")
  density <- read("paired_density")
  stopifnot(!anyDuplicated(annual[c("MPA", "year", "frame", "size", "offset_count", "series")]),
    all(is.finite(annual$contrast)), all(density$mpa >= 0), all(density$reference >= 0))
  for (d in split(annual, interaction(annual[c("MPA", "year", "frame", "size", "offset_count")], drop = TRUE))) {
    stopifnot(nrow(d) == 14L)
    v <- setNames(d$contrast, d$series)
    stopifnot(abs(v["Difference_warm"] - (v["Targeted"] - v["Nontargeted_warm"])) < 1e-10,
      abs(v["Targeted"] - mean(v[c("SPUL", "PCLA")])) < 1e-10)
  }
  a <- subset(density, size == "all_sizes"); b <- subset(density, size == "size_filtered")
  key <- function(x) do.call(paste, c(x[c("MPA", "frame", "species", "year")], sep = "|"))
  b <- b[match(key(a), key(b)), ]
  stopifnot(!anyNA(b$mpa), all(b$mpa <= a$mpa + 1e-10), all(b$reference <= a$reference + 1e-10))
  stopifnot(all(fits$n_before[fits$endpoint == "before_after"] >= 3),
    all(fits$n_after[fits$endpoint == "before_after"] >= 5),
    all(fits$n_after[fits$endpoint == "post_trend"] >= 8))
  for (m in coverage$MPA) {
    d <- subset(annual, MPA == m)
    expected <- if (nzchar(coverage$years[coverage$MPA == m]))
      as.integer(strsplit(coverage$years[coverage$MPA == m], ";", fixed = TRUE)[[1]]) else integer()
    stopifnot(setequal(d$year, expected))
  }
}
cat("PASS: fish siting synthetic, matched-cohort, direct-difference and checksum checks\n")
