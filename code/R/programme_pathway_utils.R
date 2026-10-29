# Matched programme comparisons use one post-designation level per MPA, not mixed routes.
programme_pathway_links <- function() data.frame(
  link = 1:6,
  predictor = c("S. purpuratus", "M. franciscanus", "P. interruptus", "S. pulcher", "P. interruptus", "S. pulcher"),
  predictor_response = c("Den", "Den", "Den", "Den", "Bio", "Bio"),
  outcome = c("M. pyrifera", "M. pyrifera", rep("S. purpuratus", 4)),
  outcome_response = c("Bio", "Bio", "Den", "Den", "Bio", "Bio"))

programme_pathway_annual <- function(ratios, primary) {
  required <- c("CA_MPA_Name_Short", "source", "year", "y", "resp", "BA", "lnDiff")
  if (!all(required %in% names(ratios))) stop("Missing annual ratio columns.")
  taxa <- c("Panulirus interruptus" = "P. interruptus", "Semicossyphus pulcher" = "S. pulcher",
    "Strongylocentrotus purpuratus" = "S. purpuratus", "Mesocentrotus franciscanus" = "M. franciscanus",
    "Macrocystis pyrifera" = "M. pyrifera")
  d <- data.frame(MPA = ratios$CA_MPA_Name_Short, Source = ratios$source, year = ratios$year,
    Taxa = unname(taxa[ratios$y]), Resp = ratios$resp, BA = ratios$BA, value = ratios$lnDiff)
  d <- d[d$Source %in% c("PISCO", "KFM", "LTER") & d$Resp %in% c("Bio", "Den"), ]
  if (anyNA(d$Taxa) || anyDuplicated(d[c("MPA", "Source", "Taxa", "Resp", "year")]))
    stop("Unknown taxon or duplicated annual key.")
  if (any(!is.finite(d$year)) || any(d$year != as.integer(d$year))) stop("Invalid calendar year.")
  d <- d[d$BA == "After" & is.finite(d$value) & d$MPA %in% primary$MPA, ]
  d$primary_history <- uncertainty_key(d) %in% uncertainty_key(primary)
  d
}

programme_pathway_match <- function(annual, link, sources, scope, min_years = 4L) {
  d <- annual
  if (scope == "retained_primary_histories") d <- d[d$primary_history, ]
  if (!scope %in% c("retained_primary_histories", "available_diver_histories")) stop("Unknown scope.")
  x <- d[d$Taxa == link$predictor & d$Resp == link$predictor_response & d$Source %in% sources,
         c("MPA", "Source", "year", "value", "primary_history")]
  y <- d[d$Taxa == link$outcome & d$Resp == link$outcome_response & d$Source %in% sources,
         c("MPA", "Source", "year", "value", "primary_history")]
  names(x)[4:5] <- c("x", "x_primary"); names(y)[4:5] <- c("y", "y_primary")
  paired <- dplyr::inner_join(x, y, by = c("MPA", "Source", "year"))
  common <- paired %>% dplyr::group_by(MPA, year) %>%
    dplyr::filter(dplyr::n_distinct(Source) == length(sources)) %>% dplyr::ungroup()
  coverage <- paired %>% dplyr::group_by(MPA, Source) %>%
    dplyr::summarise(n_available_years = dplyr::n(), .groups = "drop") %>%
    dplyr::left_join(common %>% dplyr::distinct(MPA, year) %>% dplyr::count(MPA, name = "n_common_years"), by = "MPA") %>%
    dplyr::mutate(n_common_years = tidyr::replace_na(n_common_years, 0L),
      eligible = n_common_years >= min_years)
  common <- common %>% dplyr::group_by(MPA) %>%
    dplyr::filter(dplyr::n_distinct(year) >= min_years) %>% dplyr::ungroup() %>% dplyr::arrange(MPA, year, Source)
  # Pool annual ratios before estimating variance, preserving their observed cofluctuation.
  pooled <- common %>% dplyr::group_by(MPA, year) %>%
    dplyr::summarise(x = mean(x), y = mean(y), x_primary = all(x_primary),
      y_primary = all(y_primary), .groups = "drop") %>% dplyr::mutate(Source = "equal_programme_mean")
  list(annual = dplyr::bind_rows(common, pooled), coverage = coverage)
}

programme_pathway_effects <- function(annual, lag) {
  if (!nrow(annual)) return(data.frame())
  dplyr::bind_rows(lapply(split(annual, interaction(annual$MPA, annual$Source, drop = TRUE)), function(d) {
    d <- d[order(d$year), ]
    xm <- stats::lm(x ~ 1, data = d); ym <- stats::lm(y ~ 1, data = d)
    data.frame(MPA = d$MPA[1], Source = d$Source[1], n_years = nrow(d),
      year_min = min(d$year), year_max = max(d$year), years = paste(d$year, collapse = ";"),
      x = mean(d$x), y = mean(d$y), x_vi = uncertainty_hac(xm, d$year, lag)[1, 1],
      vi = uncertainty_hac(ym, d$year, lag)[1, 1],
      x_primary = all(d$x_primary), y_primary = all(d$y_primary))
  }))
}

programme_pathway_fit <- function(d, min_mpas = 4L) {
  warnings <- character()
  status <- if (nrow(d) < min_mpas) "fewer_than_four_MPAs" else if (
    any(!is.finite(d$vi) | d$vi <= 0)) "nonpositive_outcome_variance" else if (
    !is.finite(stats::sd(d$x)) || stats::sd(d$x) < 1e-12) "constant_predictor" else "fitted"
  fit <- NULL
  if (status == "fitted") {
    fit <- tryCatch(withCallingHandlers(metafor::rma(yi = y, vi = vi, mods = ~ x,
      data = d, method = "REML", test = "knha"), warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
      }), error = function(e) e)
    if (inherits(fit, "error")) {status <- paste0("fit_failed: ", conditionMessage(fit)); fit <- NULL}
    if (!is.null(fit) && any(!is.finite(c(coef(fit), fit$se, fit$ci.lb, fit$ci.ub, fit$pval)))) {
      status <- "nonfinite_fit"; fit <- NULL
    }
  }
  get <- function(name) if (is.null(fit)) NA_real_ else as.numeric(fit[[name]][2])
  result <- data.frame(n_mpa = nrow(d), n_paired_mpa_years = if (nrow(d)) sum(d$n_years) else 0L,
    status, slope = if (is.null(fit)) NA_real_ else as.numeric(coef(fit)[2]),
    se = get("se"), ci_low = get("ci.lb"), ci_high = get("ci.ub"), p = get("pval"),
    df = if (is.null(fit)) NA_real_ else fit$k - fit$p,
    tau2 = if (is.null(fit)) NA_real_ else fit$tau2,
    residual_Q = if (is.null(fit)) NA_real_ else fit$QE,
    residual_Q_p = if (is.null(fit)) NA_real_ else fit$QEp,
    max_cooks_distance = if (is.null(fit)) NA_real_ else max(stats::cooks.distance(fit)),
    warnings = paste(unique(warnings), collapse = "; "))
  list(result = result, model = fit)
}

programme_pathway_reference_groups <- function(mpas, sources, registry) {
  r <- registry[registry$MPA %in% mpas & registry$Source %in% sources &
    registry$included & registry$identity_verified, ]
  groups <- setNames(as.list(mpas), mpas)
  ids <- split(r$MPA, paste(r$Source, r$reference_id, sep = "|"))
  for (shared in ids) {
    joined <- unique(unlist(groups[unique(shared)]))
    for (mpa in joined) groups[[mpa]] <- joined
  }
  data.frame(MPA = mpas, reference_group = vapply(groups, function(x) paste(sort(x), collapse = ";"), ""))
}
