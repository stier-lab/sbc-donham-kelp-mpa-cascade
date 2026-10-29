# Descriptive summaries for manuscript Table 1, using tracked annual responses.
summarise_raw_responses <- function(data) {
  required <- c("taxon_name", "source", "resp", "value")
  if (!all(required %in% names(data)) || anyNA(data[required]) ||
      any(!is.finite(data$value))) stop("Invalid raw response summary inputs.")
  dplyr::summarise(dplyr::group_by(data, taxon_name, source, resp),
    N = dplyr::n(), mean_value = mean(value), sd = stats::sd(value), .groups = "drop") %>%
    dplyr::mutate(value = mean_value, se = sd / sqrt(N),
      ci = stats::qt(.975, N - 1L) * se) %>%
    dplyr::select(taxon_name, source, resp, N, value, sd, se, ci) %>%
    dplyr::arrange(taxon_name, source, resp)
}
