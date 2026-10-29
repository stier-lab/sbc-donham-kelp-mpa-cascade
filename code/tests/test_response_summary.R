suppressPackageStartupMessages(library(dplyr))
source("code/R/response_summary_utils.R")
d <- data.frame(taxon_name = c("A", "A", "B", "B"), source = "LTER",
                resp = "Bio", value = c(2, 4, 10, 14))
r <- summarise_raw_responses(d)
stopifnot(identical(names(r), c("taxon_name", "source", "resp", "N", "value", "sd", "se", "ci")),
          all(r$N == 2L), all(r$value == c(3, 12)), all(r$se == c(1, 2)),
          isTRUE(all.equal(r$ci, stats::qt(.975, 1) * c(1, 2))))
d$value[1] <- NA_real_
stopifnot(inherits(try(summarise_raw_responses(d), silent = TRUE), "try-error"))
cat("Raw response summary tests passed.\n")
