local({
  if (!requireNamespace("here", quietly = TRUE)) stop("Package here is required.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package dplyr is required.")
  `%>%` <- dplyr::`%>%`
  source(here::here("code", "R", "review_concern_utils.R"), local = TRUE)
  set.seed(281997L)
  d <- data.frame(year = c(1984:1992, 1995:2001, 2004:2014))
  d$BA <- factor(ifelse(d$year <= 2003, "Before", "After"), c("Before", "After"))
  d$lnDiff <- 2 + 0.7 * (d$BA == "After") + rnorm(nrow(d), sd = 0.2)
  fit <- review_fit_gls(d, lnDiff ~ BA, "IID")
  stopifnot(!is.null(fit$model))
  contrast <- review_gls_contrast(fit$model, c("(Intercept)" = 0, BAAfter = 1))
  stopifnot(abs(contrast$estimate - (mean(d$lnDiff[d$BA == "After"]) - mean(d$lnDiff[d$BA == "Before"]))) < 1e-10,
            abs(contrast$se - summary(lm(lnDiff ~ BA, d))$coefficients[2, 2]) < 1e-10,
            abs(contrast$estimate - mean(d$lnDiff[d$BA == "After"])) > 1)
  ar <- review_fit_gls(d, lnDiff ~ BA, "AR1")
  stopifnot(!is.null(ar$model), all(nlme::getCovariate(ar$model$modelStruct$corStruct) == d$year))
  phi <- as.numeric(stats::coef(ar$model$modelStruct$corStruct, unconstrained = FALSE))
  stopifnot(max(abs(nlme::corMatrix(ar$model$modelStruct$corStruct) - phi^abs(outer(d$year, d$year, "-")))) < 1e-10)
  for (result in list(fit, ar)) {
    a <- result$attempts
    stopifnot(nrow(a) == if (identical(result, fit)) 2L else 10L,
      setequal(a$optimizer, c("nlminb", "optim_BFGS")), sum(a$selected) == 1L,
      a$log_lik[a$selected] == max(a$log_lik[a$success]),
      result$selection$best_log_lik == max(a$log_lik[a$success]))
  }
  stopifnot(is.na(review_gls_fit_audit(fit$model)$apvar_ok),
    grepl("not applicable", review_gls_fit_audit(fit$model)$apvar_status))
  diagnostics <- review_gls_diagnostics(ar$model, d, "BACIPS_step")
  stopifnot(diagnostics$residual_lag1_n_pairs == sum(diff(d$year) == 1L))
  failed <- review_fit_gls(d, lnDiff ~ absent_covariate, "AR1")
  stopifnot(is.null(failed$model), nrow(failed$attempts) == 10L, !any(failed$attempts$success),
    !any(failed$attempts$selected), all(nzchar(failed$attempts$error)), nzchar(failed$error))
  duplicate_failed <- tryCatch({review_fit_gls(rbind(d, d[1, ]), lnDiff ~ BA); FALSE}, error = function(e) TRUE)
  stopifnot(duplicate_failed)
  stopifnot(abs(review_signflip(c(1, 1, 1, 1)) - 2/16) < 1e-10,
            abs(review_signflip(c(-1, -1, -1, -1)) - 2/16) < 1e-10)
  groups <- review_reference_components(c("A", "B", "C", "D", "E"),
                                       c("ref1", "ref1", "ref2", "ref3", "ref2;ref3"))
  stopifnot(length(unique(groups)) == 2L, groups[1] == groups[2],
            groups[3] == groups[4], groups[4] == groups[5])
  x <- c(-0.4, 0.1, 0.5, 0.7, 1.2)
  independent <- review_group_summary(x, paste0("G", seq_along(x)))
  stopifnot(abs(independent$se - review_summary(x)$se) < 1e-10,
            !is.finite(review_group_summary(x, rep("shared", length(x)))$ci_low))
  groups <- c("shared", "shared", "B", "C", "C")
  summary <- review_group_summary(x, groups, mapping_complete = FALSE)
  deletions <- vapply(unique(groups), function(group) mean(x[groups != group]), numeric(1))
  expected_se <- sqrt(2/3 * sum((deletions - mean(deletions))^2))
  stopifnot(summary$estimate == mean(x), abs(summary$se - expected_se) < 1e-10,
    abs(summary$ci_low - (mean(x) - qt(.975, 2) * expected_se)) < 1e-10,
    summary$n_reference_groups == 3L, summary$few_reference_groups, !summary$reference_mapping_complete,
    grepl("incomplete", summary$uncertainty_status), summary$reference_group_size_ratio == 2)
  two_groups <- review_group_summary(x, c("A", "A", "B", "B", "B"))
  stopifnot(is.na(two_groups$se), is.na(two_groups$p), is.na(two_groups$ci_low), is.finite(two_groups$iid_se))
  # Native-unit and common-max scaling must cancel in the historical contrast.
  historical <- c(mpa_base = 10, ref_base = 20, mpa_event = 10, ref_event = 10)
  modern <- c(mpa_base = 30, ref_base = 60, mpa_event = 60, ref_event = 30)
  adv <- function(x) log((x[3] / x[1]) / (x[4] / x[2]))
  delta <- adv(modern) - adv(historical)
  stopifnot(abs(delta - log(4/2)) < 1e-10,
            abs(adv(modern/100) - adv(historical/100) - delta) < 1e-10)

  # Tamper only with temporary copies, never the tracked inputs or checksum manifest.
  source(here::here("code", "R", "00c_analysis_constants.R"), local = TRUE)
  root <- here::here()
  input <- review_design_input(dplyr::bind_rows(
    utils::read.csv(file.path(root, "data/harmonized/harmonized_response_ratios.csv")),
    utils::read.csv(file.path(root, "data/harmonized/harmonized_landsat_rr.csv"))),
    utils::read.csv(file.path(root, "data/sumstats_final.csv")), RESILIENCE_TAXA_SHORT)
  registry <- review_read_reference_registry(input)
  stopifnot(setequal(registry$primary_keys$series_id, unique(input$series_id)))
  fixture <- tempfile("review_registry_tests_")
  dir.create(file.path(fixture, "data/harmonized"), recursive = TRUE)
  files <- c(review_reference_analysis_files(), "data/harmonized/reference_sites.csv",
    "data/harmonized/reference_primary_keys.csv", "data/harmonized/mapping_provenance.csv")
  stopifnot(all(file.copy(file.path(root, files), file.path(fixture, files))))
  expect_error <- function(expression, pattern) {
    error <- tryCatch({ force(expression); NULL }, error = identity)
    stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
  }
  provenance_path <- file.path(fixture, "data/harmonized/mapping_provenance.csv")
  provenance <- utils::read.csv(provenance_path)
  provenance$path[provenance$kind == "source_audit"] <- paste0("/missing/source_audit/", provenance$role[provenance$kind == "source_audit"])
  utils::write.csv(provenance, provenance_path, row.names = FALSE)
  stopifnot(nrow(review_read_reference_registry(input, fixture)$reference_sites) == nrow(registry$reference_sites))
  changed <- input; changed$year[1] <- changed$year[1] + 1L
  expect_error(review_read_reference_registry(changed, fixture), "calendar coverage changed")
  input_path <- file.path(fixture, files[1])
  changed <- utils::read.csv(input_path); changed$lnDiff[1] <- changed$lnDiff[1] + 0.1
  utils::write.csv(changed, input_path, row.names = FALSE)
  expect_error(review_read_reference_registry(input, fixture), "MD5 mismatch")
  stopifnot(file.copy(file.path(root, files[1]), input_path, overwrite = TRUE))
  site_path <- file.path(fixture, "data/harmonized/reference_sites.csv")
  sites <- utils::read.csv(site_path); sites$n_reference_years[1] <- sites$n_reference_years[1] + 1L
  utils::write.csv(sites, site_path, row.names = FALSE)
  expect_error(review_read_reference_registry(input, fixture), "MD5 mismatch")
  provenance$md5[provenance$role == "reference_sites"] <- unname(tools::md5sum(site_path))
  utils::write.csv(provenance, provenance_path, row.names = FALSE)
  expect_error(review_read_reference_registry(input, fixture), "survey footprint is inconsistent")
  sites <- registry$reference_sites
  sites <- sites[sites$series_id != registry$primary_keys$series_id[1], ]
  utils::write.csv(sites, site_path, row.names = FALSE)
  provenance$md5[provenance$role == "reference_sites"] <- unname(tools::md5sum(site_path))
  utils::write.csv(provenance, provenance_path, row.names = FALSE)
  expect_error(review_read_reference_registry(input, fixture), "cover every primary series")
  cat("Review-concern estimator, optimizer, calendar-AR1, group uncertainty and stale-registry tests passed.\n")
})
