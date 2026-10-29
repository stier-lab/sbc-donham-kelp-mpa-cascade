# Common-cohort exploratory checks; primary contrasts remain unchanged.
coauthor_offset_ratio <- function(inside, outside) {
  if (length(inside) != length(outside) || any(!is.finite(c(inside, outside))) ||
      any(c(inside, outside) < 0)) stop("Paired responses must be finite and nonnegative.")
  positive <- c(inside, outside)[c(inside, outside) > 0]
  if (!length(positive)) return(rep(NA_real_, length(inside)))
  offset <- min(positive) / 2
  log((inside + offset) / (outside + offset))
}

coauthor_predator_cohort <- function(raw, primary, response) {
  required <- c("CA_MPA_Name_Short", "source", "year", "status", "taxon_name", "resp", "value", "BA")
  if (!all(required %in% names(raw))) stop("Missing raw-response columns.")
  if (anyDuplicated(raw[c("CA_MPA_Name_Short", "source", "year", "status", "taxon_name", "resp")]))
    stop("Duplicate raw response key.")
  taxa <- c("Panulirus interruptus", "Semicossyphus pulcher", "Strongylocentrotus purpuratus")
  short <- c("P. interruptus", "S. pulcher", "S. purpuratus")
  key <- function(mpa, source) paste(mpa, source, sep = "|")
  eligible <- Reduce(intersect, lapply(short, function(taxon) {
    d <- primary[primary$Taxa == taxon & primary$Resp == response, ]
    key(d$MPA, d$Source)
  }))
  d <- raw[raw$resp == response & raw$BA == "After" & raw$taxon_name %in% taxa &
             key(raw$CA_MPA_Name_Short, raw$source) %in% eligible, ]
  d <- d[is.finite(d$value) & d$value >= 0 & d$status %in% c("mpa", "reference"), ]
  names(d)[names(d) == "CA_MPA_Name_Short"] <- "MPA"
  wide <- tidyr::pivot_wider(d[c("MPA", "source", "year", "status", "taxon_name", "value")],
                            names_from = taxon_name, values_from = value)
  for (taxon in taxa) if (!taxon %in% names(wide)) wide[[taxon]] <- NA_real_
  wide <- wide[stats::complete.cases(wide[taxa]), ]
  wide$total_predator <- wide[[taxa[1]]] + wide[[taxa[2]]]
  wide$urchin <- wide[[taxa[3]]]
  paired <- tidyr::pivot_wider(wide[c("MPA", "source", "year", "status", "total_predator", "urchin")],
                              names_from = status, values_from = c(total_predator, urchin))
  cols <- c("total_predator_mpa", "total_predator_reference", "urchin_mpa", "urchin_reference")
  for (col in cols) if (!col %in% names(paired)) paired[[col]] <- NA_real_
  paired <- paired[stats::complete.cases(paired[cols]), ]
  paired <- paired %>% dplyr::group_by(MPA, source) %>%
    dplyr::filter(dplyr::n() >= 4L) %>%
    dplyr::mutate(predator_lnRR = coauthor_offset_ratio(total_predator_mpa, total_predator_reference),
                  urchin_lnRR = coauthor_offset_ratio(urchin_mpa, urchin_reference)) %>%
    dplyr::ungroup()
  paired$response <- response
  paired
}

coauthor_spearman <- function(x, y, seed = 341012L) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]; n <- length(x)
  if (n < 3L || stats::sd(x) == 0 || stats::sd(y) == 0)
    return(data.frame(n_mpa = n, rho = NA_real_, p = NA_real_, permutations = 0L, method = "not estimable"))
  rx <- rank(x); ry <- rank(y)
  statistic <- function(indices) stats::cor(rx, ry[indices])
  observed <- statistic(seq_len(n))
  if (n <= 8L) {
    permute <- function(x) {
      if (length(x) == 1L) return(matrix(x, nrow = 1L))
      do.call(rbind, lapply(seq_along(x), function(i) cbind(x[i], permute(x[-i]))))
    }
    permutations <- permute(seq_len(n))
    values <- apply(permutations, 1L, statistic)
    p <- mean(abs(values) >= abs(observed) - 1e-12)
    method <- "all permutations; descriptive, assumes exchangeable MPA values"
  } else {
    set.seed(seed)
    values <- replicate(19999L, statistic(sample.int(n)))
    p <- (1 + sum(abs(values) >= abs(observed) - 1e-12)) / (length(values) + 1)
    method <- "19999 sampled permutations plus observed; descriptive, assumes exchangeable MPA values"
  }
  data.frame(n_mpa = n, rho = observed, p = p, permutations = length(values), method = method)
}
