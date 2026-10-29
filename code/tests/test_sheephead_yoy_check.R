# Unit checks by default; --integration repeats all raw/baseline guards without writes.
root <- here::here()
source(file.path(root, "code/R/44_sheephead_yoy_check.R"))
fish <- data.frame(classcode = c(rep("SPUL", 8), "OTHER"),
  fish_tl = c(1, 9, 10, 11, NA, 0, Inf, -1, 5), count = c(2, 3, 7, 5, 4, 0, 2, 1, 8),
  site = letters[1:9], transect = seq_len(9))
stopifnot(identical(which(yoy_mask(fish, "remove_lt10")), 1:2),
          identical(which(yoy_mask(fish, "remove_le10")), 1:3),
          identical(yoy_filter_counts(fish, "baseline"), fish))
lt <- yoy_filter_counts(fish, "remove_lt10")
le <- yoy_filter_counts(fish, "remove_le10")
stopifnot(identical(lt[setdiff(names(fish), "count")], fish[setdiff(names(fish), "count")]),
          identical(le[setdiff(names(fish), "count")], fish[setdiff(names(fish), "count")]),
          all(le$count <= lt$count), all(lt$count <= fish$count),
          identical(is.na(le), is.na(fish)), nrow(le) == nrow(fish),
          identical(which(lt$count == 0), c(1L, 2L, 6L)))
stopifnot(inherits(try(yoy_equal(1, 2, "intentional mismatch"), silent = TRUE), "try-error"))
# The reused-allometry branch must equal evaluating the original row loop.
toy <- fish
toy$campus <- "VRG"
block <- expression(Fish <- tryCatch(stop("must use supplied raw")),
  Fish <- subset(Fish, campus == "VRG"), Fish.merged <- Fish,
  Fish.merged$biomass <- NA_real_,
  for (i in seq_len(nrow(Fish.merged))) {
    Fish.merged$biomass[i] <- Fish.merged$count[i] * Fish.merged$fish_tl[i]^3
  }, Fish.SPUL.All <- Fish.merged, Site.yr.fish <- Fish[c("site", "transect")])
toy_base <- yoy_aggregate_fish(globalenv(), block, root, "baseline", raw = toy)
for (scenario in c("remove_lt10", "remove_le10")) {
  original <- yoy_aggregate_fish(globalenv(), block, root, scenario, raw = toy)
  reused <- yoy_aggregate_fish(globalenv(), block, root, scenario, raw = toy,
                              prepared = toy_base$prepared)
  stopifnot(identical(original$site, reused$site),
            identical(original$roster, reused$roster))
}
skip <- run_sheephead_yoy_check(root = tempdir(), source_root = tempdir(), write = FALSE)
stopifnot(all(skip$status$status == "skipped_raw_unavailable"))

d <- data.frame(lnDiff = c(.1, .4, .5, .9), time.model = c(0, 1, 3, 7))
fit <- lm(lnDiff ~ time.model, d)
rebuilt <- list(data = d, gradient = c("(Intercept)" = 0, time.model = 11))
check <- yoy_fixed_fit(data.frame(Model = "Linear"), rebuilt, d$lnDiff)
stopifnot(isTRUE(all.equal(unname(check["estimate"]), 11 * unname(coef(fit)[2]))),
  isTRUE(all.equal(unname(check["se"]), 11 * sqrt(vcov(fit)[2, 2]))))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) <= 1L, !length(args) || args == "--integration")
status_path <- file.path(root, "tables/table_s_sheephead_yoy_status.csv")
if (file.exists(status_path) &&
    read.csv(status_path)$status[1] == "complete_baseline_verified") {
  read_output <- function(name) read.csv(file.path(root, "tables",
    paste0("table_s_sheephead_yoy_", name, ".csv")))
  contrasts <- read_output("contrasts")
  pooled <- read_output("pooled")
  protected <- read_output("protected")
  primary <- read.csv(file.path(root, "data/sumstats_final.csv"))
  baseline <- subset(contrasts, scenario == "baseline")
  j <- match(yoy_key(baseline, c("MPA", "Taxa", "Source", "Resp")),
             yoy_key(primary, c("MPA", "Taxa", "Source", "Resp")))
  stopifnot(!anyNA(j), all(table(pooled$scenario) == 9L),
            all(protected$before_md5 == protected$after_md5))
  yoy_equal(baseline$estimate, primary$Mean[j], "exported baseline contrasts")
  yoy_equal(baseline$se, primary$SE[j], "exported baseline SEs")
  annual <- read_output("annual")
  a <- subset(annual, scenario == "remove_lt10")
  b <- subset(annual, scenario == "remove_le10")
  keys <- c("CA_MPA_Name_Short", "year", "status", "resp")
  stopifnot(identical(yoy_key(a, keys), yoy_key(b, keys)),
    all(b$filtered <= a$filtered + 1e-8), all(a$filtered <= a$baseline + 1e-8),
    identical(is.na(a$filtered), is.na(a$baseline)))
}
if (length(args)) {
  result <- run_sheephead_yoy_check(root, write = FALSE)
  stopifnot(result$status$status == "complete_baseline_verified",
    all(result$protected$before_md5 == result$protected$after_md5),
    max(abs(subset(result$contrasts, scenario == "baseline")$delta)) < 1e-7,
    all(table(result$pooled$scenario) == 9L))
}
cat("PASS: cutoff indices/equality, zero roster, missingness, monotonicity, fixed-route clock",
    if (length(args)) "; exact baseline annual/contrast/SE/joint and primary preservation" else "", "\n")
