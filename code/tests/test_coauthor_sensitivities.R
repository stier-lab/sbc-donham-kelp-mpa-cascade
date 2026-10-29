library(dplyr)
library(tidyr)
source(here::here("code/R/coauthor_sensitivity_utils.R"))
source(here::here("code/R/taxon_display.R"))
stopifnot(identical(taxon_display(c("S. pulcher", "Semicossyphus pulcher", "S. purpuratus")),
                    c("B. pulcher", "Bodianus pulcher", "S. purpuratus")))
stopifnot(max(abs(coauthor_offset_ratio(c(0, 2, 3), c(1, 4, 5)) -
                 coauthor_offset_ratio(c(0, 2, 3)*1000, c(1, 4, 5)*1000))) < 1e-12)
stopifnot(all(is.na(coauthor_offset_ratio(c(0, 0), c(0, 0)))))
raw <- expand.grid(CA_MPA_Name_Short = "A", source = "PISCO", year = 2001:2006,
                   status = c("mpa", "reference"),
                   taxon_name = c("Panulirus interruptus", "Semicossyphus pulcher", "Strongylocentrotus purpuratus"))
raw$resp <- "Bio"; raw$BA <- "After"
raw$value <- ifelse(raw$status == "mpa", 10, 5)
primary <- data.frame(MPA = "A", Source = "PISCO", Resp = "Bio",
                      Taxa = c("P. interruptus", "S. pulcher", "S. purpuratus"))
raw$value[raw$year == 2001 & raw$taxon_name == "Semicossyphus pulcher" & raw$status == "reference"] <- NA
d <- coauthor_predator_cohort(raw, primary, "Bio")
stopifnot(nrow(d) == 5L, !2001L %in% d$year, all(d$total_predator_mpa == 20),
          all(d$total_predator_reference == 10), all(d$urchin_mpa == 10))
stopifnot(nrow(coauthor_predator_cohort(raw, primary[-1, ], "Bio")) == 0L)
stopifnot(inherits(try(coauthor_predator_cohort(rbind(raw, raw[1, ]), primary, "Bio"), silent=TRUE), "try-error"))
t <- coauthor_spearman(1:4, 4:1)
stopifnot(t$permutations == factorial(4), t$rho == -1, abs(t$p - 2/24) < 1e-12)
stopifnot(is.na(coauthor_spearman(rep(1, 4), 1:4)$p))
cat("Coauthor common-cohort and permutation tests passed.\n")
