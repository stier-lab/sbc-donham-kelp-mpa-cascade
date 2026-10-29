# Read-only helper tests. Never source module 08 or write production outputs.
local({
  load_helpers <- function(path) {
    env <- new.env(parent = globalenv())
    for (name in c("add_step_effect_size", "run_dw_test", "add_time_columns")) {
      definitions <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
        identical(x[[2]], as.name(name)), as.list(parse(path)))
      stopifnot(length(definitions) == 1L)
      eval(definitions[[1]], env)
    }
    env
  }
  current <- load_helpers("code/R/08_effect_sizes.R")
  original <- load_helpers("local_archive/lter_count_refit_2026-10-10/before/R/08_effect_sizes.R")
  check <- function(data, expected) {
    model <- lm(lnDiff ~ BA, data)
    before_n <- sum(data$BA == "Before")
    after_n <- sum(data$BA == "After")
    manual <- mean(data$lnDiff[data$BA == "After"]) - mean(data$lnDiff[data$BA == "Before"])
    expected_se <- summary(model)$sigma * sqrt(1 / before_n + 1 / after_n)
    result <- current$add_step_effect_size(data, "test", "test", "test", "Den", run_diagnostics = FALSE)
    stopifnot(nrow(result) == 1L, abs(result$Mean - expected) < 1e-10,
      abs(result$Mean - manual) < 1e-10, abs(result$SE - expected_se) < 1e-10,
      abs(result$CI - expected_se * qt(.975, df.residual(model))) < 1e-10,
      result$N == nrow(data), result$Model == "Step", result$Type == "pBACIPS")
    result
  }

  # Change row order, BA storage class and factor-level order independently.
  toy <- data.frame(BA = rep(c("Before", "After"), c(3L, 4L)),
                    lnDiff = c(0, 1, 2, 2, 3, 4, 5))
  for (storage in c("character", "before_first", "after_first")) {
    for (reversed in c(FALSE, TRUE)) {
      data <- if (reversed) toy[nrow(toy):1L, ] else toy
      if (storage != "character") data$BA <- factor(data$BA,
        levels = if (storage == "before_first") c("Before", "After") else c("After", "Before"))
      result <- check(data, 2.5)
      old <- original$add_step_effect_size(data, "test", "test", "test", "Den", run_diagnostics = FALSE)
      stopifnot(abs(result$SE - old$SE) < 1e-10, abs(result$CI - old$CI) < 1e-10)
    }
  }

  baseline <- read.csv("local_archive/lter_count_refit_2026-10-10/count_only/sumstats_final.csv",
                       stringsAsFactors = FALSE)
  selected <- baseline[baseline$Model == "Step", ]
  stopifnot(nrow(selected) == 6L)
  annual_path <- "data/harmonized/harmonized_response_ratios.csv"
  protected <- c(annual_path, "data/sumstats_final.csv")
  fingerprint <- tools::md5sum(protected)
  annual <- read.csv(annual_path, stringsAsFactors = FALSE)
  taxa <- c("M. franciscanus" = "Mesocentrotus franciscanus",
            "S. purpuratus" = "Strongylocentrotus purpuratus",
            "S. pulcher" = "Semicossyphus pulcher")
  for (i in seq_len(nrow(selected))) {
    row <- selected[i, ]
    data <- annual[annual$CA_MPA_Name_Short == row$MPA & annual$source == row$Source &
                   annual$resp == row$Resp & annual$y == taxa[[row$Taxa]], ]
    data <- current$add_time_columns(data)
    stopifnot(nrow(data) == row$N, is.character(data$BA))
    em <- emmeans::emmeans(lm(lnDiff ~ BA, data), ~ BA)
    old_contrast <- as.data.frame(pairs(em, reverse = TRUE))
    stopifnot(identical(old_contrast$contrast, "Before - After"),
              abs(old_contrast$estimate - row$Mean) < 1e-10)
    result <- check(data, -row$Mean)
    stopifnot(abs(result$SE - row$SE) < 1e-10, abs(result$CI - row$CI) < 1e-10,
              abs(result$SD - row$SD) < 1e-10)
    cat(row$MPA, row$Taxa, row$Source, row$Resp,
        sprintf("original=%.12f corrected=%.12f SE=%.12f CI_half=%.12f\n",
                row$Mean, result$Mean, result$SE, result$CI))
  }
  stopifnot(identical(fingerprint, tools::md5sum(protected)))
  cat("PASS: label-matched After-minus-Before across BA classes/orders; all six stored Step signs reversed; SE/CI/N unchanged; no production writes.\n")
})
