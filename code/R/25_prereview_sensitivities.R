# =============================================================================
# 25_prereview_sensitivities.R
# -----------------------------------------------------------------------------
# Sensitivity analyses added 2026-09-21 in response to the Craig pre-review
# (Donham-Stier-CA-MPA-MS-2026/reviews/craig_prereview_2026-09-21.md).
#
#   A. Residual autocorrelation: the per-combination Durbin-Watson distribution,
#      which the single pooled statistic in table_s_ar1_sensitivity.csv hides.
#   B. Cross-taxa regressions: permutation test, because each regression puts one
#      estimated effect size on another and the weights use only the response
#      variance, so the slope is attenuated and its parametric p is optimistic.
#   C. Point Conception SMR: the kelp meta-analysis refitted with the MPA that is
#      currently excluded on the response variable.
#   D. Zero-offset constant: how far the annual lnRR series move when the adaptive
#      offset (half the minimum non-zero proportion) is replaced by alternatives.
#
# Outputs: tables/table_s_prereview_*.csv
# =============================================================================

suppressMessages({
  library(dplyr); library(tidyr); library(metafor); library(nlme)
})
set.seed(20260921)

out <- function(df, name) {
  path <- file.path("tables", name)
  write.csv(df, path, row.names = FALSE)
  cat(sprintf("  wrote %s (%d rows)\n", path, nrow(df)))
}

# --- A. Residual autocorrelation -------------------------------------------
cat("A. Durbin-Watson distribution\n")
ms <- read.csv("tables/table_model_selection.csv", stringsAsFactors = FALSE)
dw <- ms %>% filter(!is.na(DW_pval))
dw_summary <- dw %>%
  group_by(Model) %>%
  summarise(n_tested = n(),
            n_p_lt_05 = sum(DW_pval < 0.05),
            pct_p_lt_05 = round(100 * mean(DW_pval < 0.05), 1),
            median_DW = round(median(DW_stat), 3), .groups = "drop") %>%
  bind_rows(tibble(Model = "ALL", n_tested = nrow(dw),
                   n_p_lt_05 = sum(dw$DW_pval < 0.05),
                   pct_p_lt_05 = round(100 * mean(dw$DW_pval < 0.05), 1),
                   median_DW = round(median(dw$DW_stat), 3)))
out(dw_summary, "table_s_prereview_dw_distribution.csv")

# --- B. Cross-taxa permutation test -----------------------------------------
# Each cross-taxa regression puts one estimated effect size on another, weights
# only by the response variance, and (for the focal relationship) is fitted to
# four points with two residual degrees of freedom. A permutation test drops the
# parametric assumptions. With k = 4 the permutation distribution is enumerated
# exhaustively, so the test is exact.
cat("B. Cross-taxa permutation test\n")
t3 <- read.csv("tables/table_03_cross_taxa_meta_regression.csv", stringsAsFactors = FALSE)
ms_es <- read.csv("tables/table_model_selection.csv", stringsAsFactors = FALSE)
agg <- ms_es %>% filter(!is.na(Mean), !is.na(SE), SE > 0) %>%
  group_by(Taxa, MPA, Resp) %>%
  summarise(es = sum(Mean / SE^2) / sum(1 / SE^2),
            se = sqrt(1 / sum(1 / SE^2)), .groups = "drop")
defs <- list(c("S. purpuratus","Den","M. pyrifera","Bio"),
             c("M. franciscanus","Den","M. pyrifera","Bio"),
             c("P. interruptus","Den","S. purpuratus","Den"),
             c("S. pulcher","Den","S. purpuratus","Den"),
             c("P. interruptus","Bio","S. purpuratus","Bio"),
             c("S. pulcher","Bio","S. purpuratus","Bio"))
perm_rows <- list()
for (i in seq_along(defs)) {
  d <- defs[[i]]
  xx <- agg %>% filter(Taxa == d[1], Resp == d[2]) %>% select(MPA, x = es)
  yy <- agg %>% filter(Taxa == d[3], Resp == d[4]) %>% select(MPA, y = es, yse = se)
  j <- inner_join(xx, yy, by = "MPA")
  if (nrow(j) < 4) next
  fitslope <- function(yv) {
    m <- try(rma(yi = yv, sei = j$yse, mods = ~ j$x, method = "REML", test = "knha"),
             silent = TRUE)
    if (inherits(m, "try-error")) NA_real_ else unname(m$b[2])
  }
  obs <- fitslope(j$y)
  k <- nrow(j)
  if (k <= 7) {                      # exhaustive: k! <= 5040
    grid <- unique(as.matrix(expand.grid(rep(list(seq_len(k)), k))))
    grid <- grid[apply(grid, 1, function(r) length(unique(r)) == k), , drop = FALSE]
    sl <- apply(grid, 1, function(r) fitslope(j$y[r]))
    mode <- "exact"
  } else {
    sl <- replicate(10000, fitslope(sample(j$y)))
    mode <- "monte-carlo(1e4)"
  }
  sl <- sl[!is.na(sl)]
  perm_rows[[length(perm_rows) + 1]] <- tibble(
    Relationship = t3$Relationship[i], k = k,
    slope_reproduced = round(obs, 4),
    slope_published = round(t3$Slope_Est[i], 4),
    reproduces_published = abs(obs - t3$Slope_Est[i]) < 0.005,
    p_parametric_published = round(t3$Slope_pval[i], 4),
    permutation_mode = mode, n_permutations = length(sl),
    p_permutation = round(mean(abs(sl) >= abs(obs) - 1e-12), 4),
    p_permutation_floor = round(1 / length(sl), 4))
}
out(bind_rows(perm_rows), "table_s_prereview_crosstaxa_permutation.csv")

# --- C. Point Conception SMR -------------------------------------------------
cat("C. Point Conception SMR inclusion\n")
land <- read.csv("data/harmonized/harmonized_landsat_rr.csv", stringsAsFactors = FALSE)
pc <- land %>% filter(CA_MPA_Name_Short == "Point Conception SMR")
if (nrow(pc) > 0) {
  after <- pc %>% filter(BA == "After")
  pc_es <- tibble(
    MPA = "Point Conception SMR", source = "Landsat", taxon = "Macrocystis pyrifera",
    n_before = sum(pc$BA == "Before"), n_after = nrow(after),
    mean_lnRR_after = round(mean(after$lnDiff, na.rm = TRUE), 4),
    se_lnRR_after = round(sd(after$lnDiff, na.rm = TRUE) / sqrt(nrow(after)), 4))
  t2 <- read.csv("tables/table_02_meta_analysis.csv", stringsAsFactors = FALSE)
  kelp <- t2 %>% filter(Taxa == "M. pyrifera", Response == "Biomass")
  pc_es$published_kelp_estimate <- round(kelp$Estimate[1], 4)
  pc_es$published_kelp_k <- kelp$k[1]
  out(pc_es, "table_s_prereview_point_conception.csv")
} else {
  cat("   SKIPPED: no Point Conception rows in harmonized_landsat_rr.csv\n")
}

# --- D. Zero-offset sensitivity ----------------------------------------------
cat("D. Zero-offset sensitivity\n")
raw <- read.csv("data/harmonized/harmonized_raw_responses.csv", stringsAsFactors = FALSE)
rr  <- read.csv("data/harmonized/harmonized_response_ratios.csv", stringsAsFactors = FALSE)

lnrr_under <- function(rule) {
  raw %>%
    group_by(CA_MPA_Name_Short, taxon_name, source, resp) %>%
    mutate(maxv = max(value, na.rm = TRUE),
           Prop = ifelse(maxv > 0, value / maxv, NA_real_)) %>%
    mutate(a = rule(Prop), PropCorr = Prop + a) %>%
    ungroup() %>%
    select(CA_MPA_Name_Short, year, taxon_name, source, resp, status, PropCorr) %>%
    pivot_wider(names_from = status, values_from = PropCorr) %>%
    filter(!is.na(mpa), !is.na(reference)) %>%
    mutate(lnRR = log(mpa) - log(reference)) %>%
    select(CA_MPA_Name_Short, year, taxon_name, source, resp, lnRR)
}
adaptive <- lnrr_under(function(p) min(p[p > 0], na.rm = TRUE) / 2)
fixed01  <- lnrr_under(function(p) 0.01)
tenth    <- lnrr_under(function(p) min(p[p > 0], na.rm = TRUE) / 10)

# Fidelity: which reconstructed rows reproduce the pipeline exactly?
fid <- rr %>% select(CA_MPA_Name_Short, year, taxon_name = y, source, resp, lnDiff) %>%
  inner_join(adaptive, by = c("CA_MPA_Name_Short","year","taxon_name","source","resp")) %>%
  mutate(exact = abs(lnDiff - lnRR) < 1e-8)
cat(sprintf("   reconstruction reproduces %d of %d pipeline lnRR values exactly (%.1f%%)\n",
            sum(fid$exact), nrow(fid), 100 * mean(fid$exact)))

keep <- fid %>% filter(exact) %>% select(CA_MPA_Name_Short, year, taxon_name, source, resp)
cmp <- keep %>%
  inner_join(adaptive, by = names(keep)) %>% rename(adaptive = lnRR) %>%
  inner_join(fixed01,  by = names(keep)) %>% rename(fixed_0.01 = lnRR) %>%
  inner_join(tenth,    by = names(keep)) %>% rename(adaptive_over_10 = lnRR)

offset_summary <- cmp %>%
  group_by(taxon_name, resp) %>%
  summarise(n = n(),
            median_abs_shift_fixed = round(median(abs(fixed_0.01 - adaptive)), 4),
            max_abs_shift_fixed = round(max(abs(fixed_0.01 - adaptive)), 4),
            median_abs_shift_tenth = round(median(abs(adaptive_over_10 - adaptive)), 4),
            max_abs_shift_tenth = round(max(abs(adaptive_over_10 - adaptive)), 4),
            .groups = "drop") %>%
  arrange(desc(median_abs_shift_fixed))
attr(offset_summary, "coverage") <- mean(fid$exact)
out(offset_summary, "table_s_prereview_offset_sensitivity.csv")
writeLines(sprintf("reconstruction_fidelity_exact_fraction,%.4f\nrows_compared,%d",
                   mean(fid$exact), nrow(cmp)),
           "tables/table_s_prereview_offset_coverage.csv")
cat("  wrote tables/table_s_prereview_offset_coverage.csv\n")
cat("done\n")
