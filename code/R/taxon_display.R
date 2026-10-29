# Only display names change; historical source names and statistical keys remain stable.
taxon_display <- function(x) {
  if (is.expression(x)) return(as.expression(lapply(x, taxon_display)))
  if (is.call(x)) {
    x[] <- lapply(as.list(x), function(node) {
      if (is.character(node) || is.call(node) || is.expression(node)) taxon_display(node) else node
    })
    return(x)
  }
  x <- gsub("Semicossyphus pulcher", "Bodianus pulcher", as.character(x), fixed = TRUE)
  gsub("S. pulcher", "B. pulcher", x, fixed = TRUE)
}
