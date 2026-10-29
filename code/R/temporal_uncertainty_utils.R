# Shared by model fitting and plotting so annual cohorts cannot drift.
temporal_primary_data <- function(input, primary) {
  taxa <- c("Panulirus interruptus" = "P. interruptus", "Semicossyphus pulcher" = "S. pulcher",
    "Strongylocentrotus purpuratus" = "S. purpuratus", "Mesocentrotus franciscanus" = "M. franciscanus",
    "Macrocystis pyrifera" = "M. pyrifera")
  key <- function(mpa, taxa, source, resp) paste(mpa, taxa, source, resp, sep = "|")
  k <- key(input$CA_MPA_Name_Short, unname(taxa[input$y]), input$source, input$resp)
  pk <- key(primary$MPA, primary$Taxa, primary$Source, primary$Resp)
  if (anyDuplicated(pk) || anyDuplicated(data.frame(k, year = input$year)))
    stop("Duplicate temporal primary keys or annual observations.")
  keep <- k %in% pk & input$source != "Landsat" & input$BA == "After" &
    input$time >= 0 & input$resp %in% c("Bio", "Den")
  d <- input[keep, ]
  d$history <- k[keep]
  if (any(!is.finite(d$lnDiff)) || any(!is.finite(d$time))) stop("Invalid temporal observations.")
  d[order(d$history, d$year), ]
}

fit_temporal_candidates <- function(d) {
  structures <- c(correlated_slopes_source = "(1 + time | CA_MPA_Name_Short) + (1 | source)",
    independent_slopes_source = "(1 + time || CA_MPA_Name_Short) + (1 | source)",
    independent_slopes = "(1 + time || CA_MPA_Name_Short)",
    intercepts_source = "(1 | CA_MPA_Name_Short) + (1 | source)",
    intercepts = "(1 | CA_MPA_Name_Short)")
  fits <- vector("list", length(structures)); audits <- vector("list", length(structures))
  for (i in seq_along(structures)) {
    warnings <- messages <- character()
    formula <- as.formula(paste("lnDiff ~ Species * time +", structures[i]))
    m <- tryCatch(withCallingHandlers(lme4::lmer(formula, data = d, REML = TRUE,
      control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") },
      message = function(w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleMessage") }),
      error = function(e) e)
    ok <- !inherits(m, "error")
    convergence <- if (ok) unlist(m@optinfo$conv$lme4$messages) else conditionMessage(m)
    converged <- ok && identical(m@optinfo$conv$opt, 0L) &&
      !any(grepl("failed|unable|degenerate|gradient|Hessian", convergence, ignore.case = TRUE))
    singular <- if (ok) lme4::isSingular(m) else NA
    covariance <- if (ok) as.data.frame(lme4::VarCorr(m)) else NULL
    audits[[i]] <- data.frame(candidate = names(structures)[i], formula = paste(deparse(formula), collapse = " "),
      converged, singular, accepted = converged && !isTRUE(singular),
      n_obs = nrow(d), n_mpa = length(unique(d$CA_MPA_Name_Short)), n_history = length(unique(d$history)),
      AIC = if (ok) AIC(m) else NA_real_, logLik = if (ok) as.numeric(logLik(m)) else NA_real_,
      variance_components = if (ok) paste(paste(covariance$grp, covariance$var1,
        signif(covariance$vcov, 7), sep = ":"), collapse = "; ") else "",
      warnings = paste(unique(c(warnings, messages, convergence)), collapse = "; "), selected = FALSE)
    if (ok) fits[[i]] <- m
  }
  audit <- do.call(rbind, audits)
  selected <- which(audit$accepted)[1]
  if (is.na(selected)) stop("No converged nonsingular temporal structure; inspect candidates.")
  audit$selected[selected] <- TRUE
  # The order is a diagnostic hierarchy, not a search for small P values or minimum AIC.
  list(model = fits[[selected]], audit = audit, candidate = names(structures)[selected],
    is_fallback = grepl("intercepts", names(structures)[selected]))
}

temporal_contrast <- function(model, species, time = 0, slope = FALSE) {
  beta <- lme4::fixef(model)
  nd <- data.frame(Species = factor(species, levels = levels(model@frame$Species)), time = time)
  formula <- reformulas::nobars(formula(model))
  x <- model.matrix(delete.response(terms(formula)), nd)
  x <- x[, names(beta), drop = FALSE]
  if (slope) {
    nd$time <- time + 1
    x <- model.matrix(delete.response(terms(formula)), nd)[, names(beta), drop = FALSE] - x
  }
  m <- if (inherits(model, "lmerModLmerTest")) model else lmerTest::as_lmerModLmerTest(model)
  test <- lmerTest::contest1D(m, as.numeric(x), ddf = "Satterthwaite", confint = TRUE)
  data.frame(estimate = test$Estimate, se = test$`Std. Error`, df = test$df,
    t_value = test$`t value`, p = test$`Pr(>|t|)`, ci_low = test$lower, ci_high = test$upper)
}

temporal_prediction_grid <- function(model, response, times = seq(0, 15, by = .1)) {
  model <- if (inherits(model, "lmerModLmerTest")) model else lmerTest::as_lmerModLmerTest(model)
  do.call(rbind, lapply(levels(model@frame$Species), function(sp) {
    do.call(rbind, lapply(times, function(t) cbind(resp = response, Species = sp, time = t,
      temporal_contrast(model, sp, t))))
  }))
}
