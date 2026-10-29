suppressPackageStartupMessages({library(dplyr); library(here)})
source(here("code/R/review_concern_utils.R"))
read <- function(s) read.csv(here("tables", paste0("table_s_coverage_cohort_", s, ".csv")))
s <- read("series"); r <- read("roster"); p <- read("pairs"); d <- read("difference")
stopifnot(nrow(s) == 145L, sum(s$n_before >= 3 & s$n_after >= 3) == 60L,
  !anyDuplicated(r[c("MPA", "Taxa", "Resp", "source_subset", "min_before")]))
for (m in c(1, 3)) {
  a <- filter(r, source_subset == "All programmes", min_before == m)
  stopifnot(nrow(a) == nrow(distinct(s, MPA, Taxa, Resp)),
    all((a$cohort == "With before-after coverage") == (a$max_before_years >= m)))
}
old <- read.csv(here("tables/table_s_design_matched_pairs.csv")) %>%
  filter(correlation == "Annual means", source_subset == "All programmes", design == "BACIPS_step")
new <- filter(p, min_before == 3, source_subset == "All programmes", metric == "Before-after change")
z <- inner_join(new, old, by = c("MPA", "Taxa", "Resp"))
stopifnot(nrow(z) == nrow(old), nrow(z) == nrow(new), max(abs(z$estimate.x - z$estimate.y)) < 1e-12)
expected <- s %>% group_by(MPA, Taxa, Resp) %>% summarise(value = mean(after_mean), .groups = "drop")
z <- inner_join(filter(p, min_before == 3, source_subset == "All programmes", metric == "Post-designation difference"), expected,
  by = c("MPA", "Taxa", "Resp"))
stopifnot(nrow(z) == nrow(expected), max(abs(z$estimate - z$value)) < 1e-12,
  all(is.na(d$ci_low[d$groups_with < 3 | d$groups_without < 3])))
z <- filter(p, min_before == 3, source_subset == "All programmes", Taxa == "M. franciscanus", Resp == "Bio", metric == "Post-designation difference")
g <- review_reference_components(z$MPA, z$reference_ids); yes <- z$cohort == "With before-after coverage"
j <- vapply(unique(g), function(id) mean(z$estimate[g != id & yes]) - mean(z$estimate[g != id & !yes]), numeric(1))
se <- sqrt((length(j) - 1) / length(j) * sum((j - mean(j))^2))
check <- filter(d, min_before == 3, source_subset == "All programmes", region_scope == "All", Taxa == "M. franciscanus", Resp == "Bio")
stopifnot(check$cross_cohort_groups == 1, abs(check$se - se) < 1e-12)
inputs <- read("inputs")
stopifnot(identical(unname(tools::md5sum(here(inputs$path))), inputs$md5))
cat("PASS: grouping, source means, matched-history agreement, sparse intervals, shared-reference deletion and input integrity\n")
