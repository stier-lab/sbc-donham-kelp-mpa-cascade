# Matched PISCO size refits. These never replace the primary biomass reconstruction.
size_refit_annual <- function(root) {
  sheep <- read.csv(file.path(root, "tables/table_s_sheephead_size_annual.csv"))
  lobster <- read.csv(file.path(root, "tables/table_s_lobster_size_annual.csv"))
  fish <- sheep %>% transmute(MPA = CA_MPA_Name_Short, year, Taxa = "S. pulcher",
    Resp = ifelse(measure == "density", "Den", "Bio"),
    size = ifelse(threshold_cm == 0, "all", "large"), mpa, reference, lnDiff = lnRR)
  lob <- lobster %>% transmute(MPA = CA_MPA_Name_Short, year, Taxa = "P. interruptus",
    Resp = ifelse(response == "density", "Den", "Bio"),
    size = ifelse(grepl("^all_", measure), "all", "large"), mpa, reference, lnDiff = lnRR)
  d <- bind_rows(fish, lob)
  keys <- c("MPA", "year", "Taxa", "Resp", "size")
  if (anyDuplicated(d[keys]) || any(!is.finite(as.matrix(d[c("mpa", "reference", "lnDiff")]))))
    stop("Invalid size-refit annual keys/values")
  if (any(d$mpa < 0 | d$reference < 0)) stop("Negative survey values")
  paired <- d %>% select(all_of(keys), mpa, reference) %>%
    pivot_wider(names_from = size, values_from = c(mpa, reference))
  if (anyNA(paired) || any(paired$mpa_large > paired$mpa_all + 1e-10) ||
      any(paired$reference_large > paired$reference_all + 1e-10)) stop("Unmatched size cohorts")
  # The offset is shared across thresholds within MPA/response, as in module 38.
  rebuilt <- d %>% group_by(MPA, Taxa, Resp) %>%
    mutate(check = coauthor_offset_ratio(mpa, reference)) %>% ungroup()
  if (max(abs(rebuilt$check - rebuilt$lnDiff)) > 1e-10) stop("Stored size ratios disagree")
  starts <- read.csv(file.path(root, "data/harmonized/harmonized_site_metadata.csv")) %>%
    distinct(MPA = CA_MPA_Name_Short, MPA_Start)
  if (anyDuplicated(starts$MPA)) stop("Conflicting implementation dates")
  d <- left_join(d, starts, by = "MPA") %>% mutate(time = year - MPA_Start, source = "PISCO")
  if (anyNA(d$time) || any(d$time <= 0)) stop("Expected post-designation size observations")
  d
}

size_refit_temporal_data <- function(predators, resource, primary, scope, threshold, response) {
  p <- filter(predators, size == threshold, Resp == response)
  if (scope == "retained_membership") {
    keys <- primary %>% filter(Source == "PISCO") %>% distinct(MPA, Taxa, Resp)
    p <- semi_join(p, keys, by = c("MPA", "Taxa", "Resp"))
  }
  r <- resource %>% filter(resp == response, !Species %in% c("P. interruptus", "S. pulcher"),
    CA_MPA_Name_Short %in% p$MPA) %>%
    transmute(MPA = CA_MPA_Name_Short, year, Taxa = as.character(Species), Resp = resp,
      size = threshold, lnDiff, time, source)
  d <- bind_rows(select(p, MPA, year, Taxa, Resp, size, lnDiff, time, source), r)
  d <- d %>% mutate(CA_MPA_Name_Short = MPA, Species = factor(Taxa),
    history = paste(MPA, Taxa, source, Resp, sep = "|")) %>% arrange(history, year)
  if (anyDuplicated(d[c("history", "year")]) || any(!is.finite(d$lnDiff))) stop("Invalid model cohort")
  d
}

size_refit_temporal <- function(d) {
  structures <- c(correlated_slopes = "(1 + time | CA_MPA_Name_Short)",
    independent_slopes = "(1 + time || CA_MPA_Name_Short)", intercepts = "(1 | CA_MPA_Name_Short)")
  fits <- list(); audits <- list()
  for (i in seq_along(structures)) {
    warnings <- character()
    formula <- as.formula(paste("lnDiff ~ Species * time +", structures[i]))
    m <- tryCatch(withCallingHandlers({
      m <- lme4::lmer(formula, data = d, REML = TRUE,
        control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
      if (!lme4::isSingular(m)) m <- lmerTest::as_lmerModLmerTest(m)
      m
    }, warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")},
       message = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleMessage")}),
      error = function(e) e)
    ok <- !inherits(m, "error")
    converged <- ok && identical(m@optinfo$conv$opt, 0L) &&
      !any(grepl("failed|negative eigenvalue|degenerate|gradient|Hessian", c(warnings, unlist(m@optinfo$conv$lme4$messages)), ignore.case = TRUE))
    singular <- if (ok) lme4::isSingular(m) else NA
    valid <- converged && !isTRUE(singular)
    audits[[i]] <- data.frame(candidate = names(structures)[i], formula = paste(deparse(formula), collapse = " "),
      converged, singular, accepted = valid, n_obs = nrow(d), n_mpa = n_distinct(d$MPA),
      n_history = n_distinct(d$history), warnings = if (ok) paste(unique(warnings), collapse = "; ") else conditionMessage(m), selected = FALSE)
    fits[[i]] <- m
    if (valid) break
  }
  audit <- bind_rows(audits); chosen <- which(audit$accepted)[1]
  if (is.na(chosen)) stop("No valid size temporal fit; inspect covariance/Hessian")
  audit$selected[chosen] <- TRUE
  list(model = fits[[chosen]], candidate = names(structures)[chosen], audit = audit)
}

size_refit_ar1 <- function(d) {
  d$year_f <- factor(d$year, levels = seq.int(min(d$year), max(d$year)))
  d$history <- factor(d$history)
  structures <- c(independent_slopes = "(1 + time || CA_MPA_Name_Short)",
    intercepts = "(1 | CA_MPA_Name_Short)")
  fits <- list(); audits <- list()
  for (i in seq_along(structures)) {
    warnings <- character()
    m <- tryCatch(withCallingHandlers(glmmTMB::glmmTMB(
      as.formula(paste("lnDiff ~ Species * time +", structures[i], "+ ar1(year_f + 0 | history)")),
      data = d, family = gaussian(), REML = TRUE,
      control = glmmTMB::glmmTMBControl(optCtrl = list(iter.max = 1500, eval.max = 2000))),
      warning = function(w) {warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")}),
      error = function(e) e)
    valid <- !inherits(m, "error") && m$fit$convergence == 0 && isTRUE(m$sdr$pdHess) &&
      all(is.finite(vcov(m)$cond))
    audits[[i]] <- data.frame(candidate = names(structures)[i], valid,
      convergence = if (inherits(m, "error")) NA_integer_ else m$fit$convergence,
      pdHess = if (inherits(m, "error")) FALSE else isTRUE(m$sdr$pdHess),
      max_gradient = if (valid) max(abs(m$sdr$gradient.fixed)) else NA_real_,
      warnings = if (inherits(m, "error")) conditionMessage(m) else paste(unique(warnings), collapse = "; "),
      selected = FALSE)
    fits[[i]] <- m
    if (valid) break
  }
  audit <- bind_rows(audits)
  chosen <- which(audit$valid)[1]
  if (is.na(chosen)) return(list(model = NULL, audit = audit))
  audit$selected[chosen] <- TRUE
  list(model = fits[[chosen]], audit = audit)
}

size_refit_ar1_slope <- function(m, species) {
  beta <- glmmTMB::fixef(m)$cond; v <- vcov(m)$cond
  nd <- data.frame(Species = factor(species, levels = levels(m$frame$Species)), time = c(0, 1))
  x <- model.matrix(~ Species * time, nd)[, names(beta), drop = FALSE]
  g <- x[2, ] - x[1, ]; estimate <- sum(g * beta)
  se <- sqrt(drop(t(g) %*% v %*% g))
  data.frame(estimate, se, ci_low = estimate - qnorm(.975) * se,
    ci_high = estimate + qnorm(.975) * se, p = 2 * pnorm(-abs(estimate / se)))
}

size_refit_cluster_slope <- function(m, species) {
  v <- clubSandwich::vcovCR(m, cluster = m@frame$CA_MPA_Name_Short, type = "CR2")
  nd <- data.frame(Species = factor(species, levels = levels(m@frame$Species)), time = c(0, 1))
  x <- model.matrix(~ Species * time, nd)[, names(lme4::fixef(m)), drop = FALSE]
  g <- x[2, ] - x[1, ]
  test <- clubSandwich::Wald_test(m, constraints = matrix(g, nrow = 1), vcov = v, test = "HTZ")
  estimate <- sum(g * lme4::fixef(m)); se <- sqrt(drop(t(g) %*% v %*% g))
  df <- test$df_denom
  data.frame(estimate, se, df, t_value = estimate / se, p = test$p_val,
    ci_low = estimate - qt(.975, df) * se, ci_high = estimate + qt(.975, df) * se)
}

size_refit_ar1_predictions <- function(m, response) {
  beta <- glmmTMB::fixef(m)$cond; v <- vcov(m)$cond
  bind_rows(lapply(levels(m$frame$Species), function(sp) {
    nd <- data.frame(Species = factor(sp, levels = levels(m$frame$Species)), time = seq(0, 15, .5))
    x <- model.matrix(~ Species * time, nd)[, names(beta), drop = FALSE]
    estimate <- drop(x %*% beta); se <- sqrt(rowSums((x %*% v) * x))
    data.frame(resp = response, Species = sp, time = nd$time, estimate, se,
      ci_low = estimate - qnorm(.975) * se, ci_high = estimate + qnorm(.975) * se)
  }))
}

size_refit_history <- function(value, year, model = "Mean") {
  if (length(value) < 4 || anyDuplicated(year)) stop("Insufficient/nonunique paired years")
  time <- year - min(year)
  fit <- if (model == "Linear") lm(value ~ time) else lm(value ~ 1)
  g <- if (model == "Linear") c(0, 11) else 1
  estimate <- sum(g * coef(fit))
  iid <- drop(t(g) %*% vcov(fit) %*% g)
  hac <- drop(t(g) %*% uncertainty_hac(fit, year, lag = 1) %*% g)
  variance <- max(iid, hac)
  data.frame(estimate, variance, n_years = length(year), first_year = min(year),
    last_year = max(year), model, extrapolated_11 = model == "Linear" && max(time) < 11)
}

size_refit_pathways <- function(predators, resources, primary) {
  rows <- annual <- effects <- omissions <- list(); counter <- 0L
  for (scope in c("retained_membership", "expanded_size_cohort")) {
    for (response in c("Den", "Bio")) for (taxon in c("S. pulcher", "P. interruptus")) {
      p <- filter(predators, Taxa == taxon, Resp == response)
      if (scope == "retained_membership") p <- semi_join(p,
        primary %>% filter(Source == "PISCO") %>% distinct(MPA, Taxa, Resp), by = c("MPA", "Taxa", "Resp"))
      p <- p %>% select(MPA, year, size, lnDiff) %>% pivot_wider(names_from = size, values_from = lnDiff)
      r <- resources %>% filter(resp == response, Species == "S. purpuratus") %>%
        transmute(MPA = CA_MPA_Name_Short, year, resource = lnDiff)
      d <- inner_join(p, r, by = c("MPA", "year")) %>% group_by(MPA) %>% filter(n() >= 4) %>% ungroup()
      if (anyDuplicated(d[c("MPA", "year")])) stop("Duplicated pathway annual observations")
      annual[[length(annual) + 1L]] <- cbind(scope, Taxa = taxon, Resp = response, d)
      for (estimand in c("post_only_mean", "common_CI_route")) {
        e <- bind_rows(lapply(split(d, d$MPA), function(z) {
          select_route <- function(value) {
            f <- lm(value ~ year, data = z)
            prob <- coef(summary(f))[2, 4]
            if (is.finite(prob) && prob <= .05) "Linear" else "Mean"
          }
          xm <- if (estimand == "common_CI_route") select_route(z$all) else "Mean"
          ym <- if (estimand == "common_CI_route") select_route(z$resource) else "Mean"
          y <- size_refit_history(z$resource, z$year, ym)
          bind_rows(lapply(c("all", "large"), function(threshold) {
            x <- size_refit_history(z[[threshold]], z$year, xm)
            data.frame(MPA = z$MPA[1], size = threshold, x = x$estimate, x_vi = x$variance,
              y = y$estimate, vi = y$variance, n_years = nrow(z), years = paste(z$year, collapse = ";"),
              x_model = xm, y_model = ym, x_extrapolated_11 = x$extrapolated_11,
              y_extrapolated_11 = y$extrapolated_11)
          }))
        }))
        effects[[length(effects) + 1L]] <- cbind(scope, Taxa = taxon, Resp = response, estimand, e)
        for (threshold in c("all", "large")) {
          z <- filter(e, size == threshold)
          counter <- counter + 1L
          status <- if (nrow(z) < 4) "fewer_than_four_MPAs" else if (any(z$vi <= 0))
            "zero_outcome_variance" else if (sd(z$x) == 0) "constant_predictor" else "fitted"
          fit <- if (status == "fitted") tryCatch(metafor::rma(y, vi = vi, mods = ~ x,
            data = z, method = "REML", test = "knha"), error = function(e) e) else NULL
          if (inherits(fit, "error")) {status <- conditionMessage(fit); fit <- NULL}
          rows[[counter]] <- data.frame(scope, Taxa = taxon, Resp = response, estimand, size = threshold,
            n_mpa = nrow(z), n_paired_years = sum(z$n_years), status,
            slope = if (is.null(fit)) NA_real_ else as.numeric(coef(fit)[2]),
            se = if (is.null(fit)) NA_real_ else fit$se[2],
            ci_low = if (is.null(fit)) NA_real_ else fit$ci.lb[2],
            ci_high = if (is.null(fit)) NA_real_ else fit$ci.ub[2],
            p = if (is.null(fit)) NA_real_ else fit$pval[2],
            n_extrapolated = sum(z$x_extrapolated_11 | z$y_extrapolated_11))
          omit_sets <- c(setNames(as.list(z$MPA), z$MPA), list(shared_Naples_Campus = c("Naples SMCA", "Campus Point SMCA")))
          if (!is.null(fit) && nrow(z) >= 5) for (omit_name in names(omit_sets)) {
            omit <- omit_sets[[omit_name]]
            sub <- z[!z$MPA %in% omit, ]
            m <- tryCatch(metafor::rma(y, vi = vi, mods = ~ x, data = sub,
              method = "REML", test = "knha"), error = function(e) NULL)
            omissions[[length(omissions) + 1L]] <- data.frame(scope, Taxa = taxon, Resp = response,
              estimand, size = threshold, omitted = omit_name, n_mpa = nrow(sub),
              slope = if (is.null(m)) NA_real_ else as.numeric(coef(m)[2]),
              ci_low = if (is.null(m)) NA_real_ else m$ci.lb[2], ci_high = if (is.null(m)) NA_real_ else m$ci.ub[2])
          }
        }
      }
    }
  }
  summary <- bind_rows(rows) %>% group_by(scope, estimand, size) %>%
    mutate(q = p.adjust(p, "BH", n = 4)) %>% ungroup() %>%
    mutate(inference = "Conditional exploratory: ignores moderator/size-pool uncertainty and between-MPA shared references",
      quantity = ifelse(estimand == "post_only_mean", "Mean post-designation log ratio",
        "Legacy route: mixture of mean levels and 11-year changes; conditional on IID selection; route locked across sizes"))
  list(summary = summary, annual = bind_rows(annual), effects = bind_rows(effects), omissions = bind_rows(omissions))
}
