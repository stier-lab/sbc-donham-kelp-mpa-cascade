suppressPackageStartupMessages({library(dplyr); library(tidyr)})
source(here::here("code/R/primary_uncertainty_utils.R"))
source(here::here("code/R/programme_pathway_utils.R"))
read_table <- function(stem) read.csv(here::here("tables", paste0("table_s_programme_pathway_", stem, ".csv")))

# Mismatched years and endpoints cannot enter the common programme cohort.
toy <- expand.grid(MPA = c("A", "B"), Source = c("PISCO", "KFM"), year = 2001:2006,
  Taxa = c("S. pulcher", "S. purpuratus"), stringsAsFactors = FALSE) %>%
  mutate(Resp = "Den", BA = "After", value = sin(year) + as.numeric(factor(MPA)), primary_history = TRUE)
toy <- toy[!(toy$MPA == "A" & toy$Source == "KFM" & toy$year == 2001 & toy$Taxa == "S. purpuratus"), ]
match_toy <- programme_pathway_match(toy, programme_pathway_links()[4, ], c("KFM", "PISCO"), "retained_primary_histories")
stopifnot(!any(match_toy$annual$MPA == "A" & match_toy$annual$year == 2001),
  all((match_toy$annual %>% count(MPA, year))$n == 3L))
one <- programme_pathway_effects(match_toy$annual, 3)
for (mpa in c("A", "B")) {
  z <- one[one$MPA == mpa, ]
  # Identical programme measurements are not independent replicates: variance must not halve.
  stopifnot(abs(z$vi[z$Source == "equal_programme_mean"] - z$vi[z$Source == "KFM"]) < 1e-12,
    abs(z$x[z$Source == "equal_programme_mean"] - mean(z$x[z$Source != "equal_programme_mean"])) < 1e-12)
}
toy$primary_history[toy$MPA == "B" & toy$Source == "KFM"] <- FALSE
strict <- programme_pathway_match(toy, programme_pathway_links()[4, ], c("KFM", "PISCO"), "retained_primary_histories")
broad <- programme_pathway_match(toy, programme_pathway_links()[4, ], c("KFM", "PISCO"), "available_diver_histories")
stopifnot(!"B" %in% strict$annual$MPA, "B" %in% broad$annual$MPA)
stopifnot(programme_pathway_fit(data.frame(x = 1:2, y = 2:3, vi = c(1, 1), n_years = c(4, 4)))$result$status == "fewer_than_four_MPAs")
stopifnot(programme_pathway_fit(data.frame(x = 1:4, y = 2:5, vi = c(0, 1, 1, 1), n_years = rep(4, 4)))$result$status == "nonpositive_outcome_variance")

status <- read_table("status"); stopifnot(all(status$status == "completed"), !any(status$primary_changed))
hashes <- read_table("inputs")
stopifnot(identical(unname(tools::md5sum(here::here(hashes$path))), as.character(hashes$md5)))
fits <- read_table("fits"); e <- read_table("effects"); annual <- read_table("annual")
stopifnot(nrow(fits) == 216L, !anyDuplicated(fits[c("scope", "programme_pair", "link", "lag_years", "Source")]),
  !anyDuplicated(annual[c("scope", "programme_pair", "link", "MPA", "Source", "year")]),
  all((annual %>% count(scope, programme_pair, link, MPA, year))$n == 3L))
check_effects <- e %>% group_by(scope, programme_pair, link, lag_years, MPA) %>%
  summarise(n_sources = n(), n_year_sets = n_distinct(years), .groups = "drop")
stopifnot(all(check_effects$n_sources == 3), all(check_effects$n_year_sets == 1))
for (g in split(fits, interaction(fits$scope, fits$programme_pair, fits$lag_years, fits$Source)))
  stopifnot(isTRUE(all.equal(g$q, p.adjust(g$p, "BH", n = 6), check.attributes = FALSE)))
fitted <- fits[fits$status == "fitted", ]
stopifnot(all(fitted$n_mpa >= 4), all(fitted$df == fitted$n_mpa - 2), all(fitted$ci_low < fitted$ci_high))
for (j in seq_len(nrow(fitted))) {
  f <- fitted[j, ]
  d <- e %>% filter(scope == f$scope, programme_pair == f$programme_pair, link == f$link,
    lag_years == f$lag_years, Source == f$Source)
  m <- programme_pathway_fit(d)$result
  stopifnot(abs(f$slope - m$slope) < 1e-10, abs(f$ci_low - m$ci_low) < 1e-10)
}
lag3 <- fits %>% filter(scope == "available_diver_histories", programme_pair == "KFM_PISCO", lag_years == 3, link <= 4)
stopifnot(nrow(lag3) == 12, all(lag3$n_mpa == 5), all(lag3$status == "fitted"),
  all(lag3$ci_low < 0 & lag3$ci_high > 0))
stopifnot(all(fits$status[fits$programme_pair %in% c("KFM_LTER", "LTER_PISCO")] == "fewer_than_four_MPAs"))
omissions <- read_table("omissions")
stopifnot(all(omissions$n_mpa < 5), all(omissions$n_mpa[omissions$status == "fitted"] >= 4))
for (g in split(omissions, interaction(omissions$scope, omissions$programme_pair,
    omissions$lag_years, omissions$Source, omissions$omitted, drop = TRUE)))
  stopifnot(isTRUE(all.equal(g$q, p.adjust(g$p, "BH", n = 6), check.attributes = FALSE)))
registry <- data.frame(MPA = c("A", "B", "B", "C", "D"), Source = "PISCO",
  reference_id = c("one", "one", "two", "two", "other"), included = TRUE, identity_verified = TRUE)
groups <- programme_pathway_reference_groups(c("A", "B", "C", "D"), c("PISCO", "KFM"), registry)
stopifnot(all(groups$reference_group[groups$MPA != "D"] == "A;B;C"), groups$reference_group[groups$MPA == "D"] == "D")
cat("PASS: endpoint/year matching, scope separation, pooling covariance, sparse statuses, input hashes, BH families and exact model reproduction\n")
