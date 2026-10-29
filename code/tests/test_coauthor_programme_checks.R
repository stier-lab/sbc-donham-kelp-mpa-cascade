library(dplyr)
library(tidyr)
read_table <- function(name) read.csv(here::here("tables", paste0(name, ".csv")))
p <- read.csv(here::here("data/sumstats_final.csv")) %>% filter(Primary == "Y")
z <- read_table("table_s_programme_omissions")
stopifnot(nrow(z) == 45L, !anyDuplicated(z[c("scenario", "Taxa", "Resp")]))
baseline <- z %>% filter(scenario == "retained_primary") %>%
  select(Taxa, Resp, primary_estimate = estimate, primary_q = q)
comparisons <- z %>% inner_join(baseline, by = c("Taxa", "Resp"))
unchanged <- comparisons %>% filter(scenario %in% c("omit_kfm_sheephead", "omit_skunk_point", "omit_gull_kfm"))
stopifnot(all(sign(unchanged$estimate) == sign(unchanged$primary_estimate)),
  all((unchanged$q < 0.05) == (unchanged$primary_q < 0.05)))
no_fish <- z %>% filter(scenario == "omit_kfm_sheephead", Taxa == "S. pulcher", Resp == "Den")
stopifnot(no_fish$k == sum(p$Taxa == "S. pulcher" & p$Resp == "Den" & p$Source != "KFM"))
links <- read_table("table_s_extended_pathways")
points <- read_table("table_s_extended_pathway_pairs")
stopifnot(nrow(links) == 8L, all(links$n_mpa >= 3),
  all(links$q >= .05),
  max(abs(links$q - p.adjust(links$p, "BH"))) < 1e-12,
  !anyDuplicated(points[c("link", "MPA")]),
  identical(as.integer(table(points$link)), as.integer(links$n_mpa)))
overlap <- links %>% filter(predator=="P. interruptus",prey=="S. purpuratus",
  predictor_response=="Bio",outcome_response=="Bio")
primary_link <- read_table("table_03_cross_taxa_meta_regression") %>%
  filter(Predictor=="P. interruptus biomass",Response_Var=="S. purpuratus biomass")
stopifnot(nrow(overlap)==1,nrow(primary_link)==1,overlap$ci_high<0,
  abs(overlap$slope-primary_link$Slope_Est)<5e-5,
  abs(overlap$p-primary_link$Slope_pval)<1e-10,
  sum(links$ci_high<0 | links$ci_low>0)==1)
fish <- read_table("table_s_sheephead_programme_common_years")
stopifnot(!anyDuplicated(fish[c("CA_MPA_Name_Short", "source", "year")]),
  all((fish %>% count(CA_MPA_Name_Short, year))$n == 2L))
if (file.exists(here::here("data/PISCO/MLPA_kelpforest_fish_2024.csv"))) {
  x <- read_table("table_s_sheephead_size_annual") %>%
    select(CA_MPA_Name_Short, year, measure, threshold_cm, mpa, reference) %>%
    pivot_wider(names_from = threshold_cm, values_from = c(mpa, reference))
  stopifnot(nrow(x) > 0, all(is.finite(as.matrix(x[c("mpa_0", "mpa_30", "reference_0", "reference_30")]))),
    all(x$mpa_30 <= x$mpa_0 + 1e-10), all(x$reference_30 <= x$reference_0 + 1e-10))
  l <- read_table("table_s_lobster_size_coverage")
  stopifnot(all(l$large_density <= l$all_density + 1e-10),
    all(l$large_biomass <= l$all_biomass + 1e-10),
    all(l$large_density[l$all_density == 0] == 0),
    all(l$large_fraction[is.finite(l$large_fraction)] >= 0 & l$large_fraction[is.finite(l$large_fraction)] <= 1))
}
cat("PASS: primary refits, omissions, BH families, one-MPA regressions, matched years and size bounds\n")
