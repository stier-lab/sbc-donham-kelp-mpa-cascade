# Exact retained estimands with calendar-gap HAC uncertainty.
uncertainty_key <- function(d) do.call(paste, c(d[c("MPA", "Taxa", "Source", "Resp")], sep = "|"))

uncertainty_annual <- function(root, primary) {
  taxa <- c("Panulirus interruptus" = "P. interruptus", "Semicossyphus pulcher" = "S. pulcher",
    "Strongylocentrotus purpuratus" = "S. purpuratus", "Mesocentrotus franciscanus" = "M. franciscanus",
    "Macrocystis pyrifera" = "M. pyrifera")
  read <- function(name) {
    x <- read.csv(file.path(root, "data/harmonized", name))
    data.frame(MPA = x$CA_MPA_Name_Short, Taxa = unname(taxa[x$y]), Source = x$source,
      Resp = x$resp, year = x$year, BA = x$BA, lnDiff = x$lnDiff, input_time = x$time)
  }
  d <- rbind(read("harmonized_response_ratios.csv"), read("harmonized_landsat_rr.csv"))
  d$key <- uncertainty_key(d)
  if (anyDuplicated(d[c("key", "year")])) stop("Duplicate annual history/year.")
  d <- d[d$key %in% uncertainty_key(primary), ]
  if (any(!is.finite(d$lnDiff)) || anyNA(d$Taxa)) stop("Invalid retained annual observations.")
  d
}

uncertainty_definitions <- function(root) {
  e <- new.env(parent = globalenv())
  for (file in c("08_effect_sizes.R", "02_pBACIPS_function.R")) {
    for (x in as.list(parse(file.path(root, "code/R", file)))) {
      if (is.call(x) && identical(x[[1]], as.name("<-")) && is.call(x[[3]]) &&
          identical(x[[3]][[1]], as.name("function"))) eval(x, e)
    }
  }
  e$MODEL_FIT_LOG <- list()
  e$nls.lm <- minpack.lm::nls.lm
  e$nls.lm.control <- minpack.lm::nls.lm.control
  e$nls2 <- nls2::nls2
  e
}

uncertainty_reconstruct <- function(row, annual, definitions, horizon = 11) {
  d <- annual[annual$key == uncertainty_key(row), ]
  d <- d[order(d$year), ]
  if (nrow(d) != row$N) stop("Retained N mismatch: ", uncertainty_key(row))
  d$CA_MPA_Name_Short <- d$MPA
  d <- definitions$add_time_columns(d)
  supplied <- row$AnalysisType == "CI" || row$Source == "Landsat" ||
    uncertainty_key(row) == "Campus Point SMCA|P. interruptus|LTER|Bio"
  if (supplied) d$time.model <- d$input_time
  d$BA <- factor(d$BA, c("Before", "After"))
  fit <- switch(row$Model,
    Mean = lm(lnDiff ~ 1, data = d),
    Linear = lm(lnDiff ~ time.model, data = d),
    Step = lm(lnDiff ~ BA, data = d),
    Sigmoid = definitions$mySIGfun_standalone(d$lnDiff, d$time.model,
      max(which(d$time.model == 0)), d$time.true),
    stop("Unsupported retained model: ", row$Model))
  if (is.null(fit) || (row$Model == "Sigmoid" && !inherits(fit, "nls")))
    stop("Retained model reconstruction failed.")
  if (inherits(fit, "nls") && !isTRUE(fit$convInfo$isConv)) stop("Nonconverged retained nonlinear fit.")
  theta <- coef(fit)
  if (row$Model == "Mean") g <- setNames(1, names(theta))
  if (row$Model == "Linear") g <- setNames(c(0, horizon), names(theta))
  if (row$Model == "Step") g <- setNames(c(0, 1), names(theta))
  if (row$Model == "Sigmoid") {
    nd <- data.frame(time.model = c(0, horizon), time_offset = c(.01, horizon + .01),
      time.true = max(which(d$time.model == 0)) + c(0, horizon))
    f <- function(par, data) {
      e <- list2env(c(as.list(data), as.list(par)), parent = fit$m$getEnv())
      as.numeric(eval(fit$m$formula()[[3]], e))
    }
    g <- drop(numDeriv::jacobian(function(par) diff(f(setNames(par, names(theta)), nd)), theta))
    names(g) <- names(theta)
  }
  estimate <- if (row$Model == "Sigmoid") diff(as.numeric(predict(fit, nd))) else sum(g * theta)
  iid_se <- sqrt(drop(t(g) %*% vcov(fit) %*% g))
  if (abs(estimate - row$Mean) > 1e-7 || abs(iid_se - row$SE) > max(1e-7, row$SE * 1e-5))
    stop("Retained contrast/SE reconstruction mismatch: ", uncertainty_key(row))
  list(fit = fit, data = d, gradient = g, estimate = estimate, iid_se = iid_se,
    clock = if (row$Model == "Mean") "not applicable" else if (supplied) "supplied calendar time" else "after observation index")
}

uncertainty_hac <- function(fit, year, lag) {
  if (anyDuplicated(year) || any(!is.finite(year)) || any(year != as.integer(year)))
    stop("HAC needs unique integer calendar years.")
  scores <- sandwich::estfun(fit)
  n <- nrow(scores); p <- ncol(scores)
  if (length(year) != n || n <= p) stop("Invalid HAC sample/parameter count.")
  span <- seq.int(min(year), max(year))
  padded <- matrix(0, length(span), p)
  padded[match(year, span), ] <- scores
  # Zero scores at unobserved years preserve actual lag distances, not imputed data.
  registerS3method("estfun", "mpa_calendar_scores", function(x, ...) x$scores,
                  envir = asNamespace("sandwich"))
  object <- structure(list(scores = padded), class = "mpa_calendar_scores")
  lag <- min(as.integer(lag), length(span) - 1L)
  weights <- 1 - (0:lag) / (lag + 1)
  meat <- sandwich::meatHAC(object, weights = weights, prewhite = FALSE, adjust = FALSE)
  inverse <- sandwich::bread(fit) / n
  v <- inverse %*% (meat * length(span)) %*% inverse * n / (n - p)
  if (any(!is.finite(v)) || min(eigen(v, symmetric = TRUE, only.values = TRUE)$values) < -1e-8)
    stop("Invalid HAC covariance.")
  v
}

uncertainty_joint <- function(primary, se, scenario) {
  d <- primary; d$vi <- se^2
  result <- do.call(rbind, lapply(c("Bio", "Den"), function(resp) {
    x <- d[d$Resp == resp, ]; x$Taxa <- factor(x$Taxa); x$Source <- factor(x$Source)
    m <- metafor::rma.mv(Mean, V = vi, mods = ~ Taxa - 1,
      random = list(~1 | MPA, ~1 | Source), method = "REML", test = "t", data = x)
    data.frame(scenario, Taxa = sub("^Taxa", "", names(coef(m))), Resp = resp,
      estimate = as.numeric(coef(m)), se = m$se, ci_low = m$ci.lb, ci_high = m$ci.ub, p = m$pval)
  }))
  result$q <- p.adjust(result$p, "BH")
  result
}
