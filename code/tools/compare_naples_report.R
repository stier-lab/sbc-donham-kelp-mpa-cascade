# A report-style calendar trend on current inputs, not a report reproduction.
d <- read.csv("data/harmonized/harmonized_response_ratios.csv")
r <- read.csv("data/external/source_reconciliation_2026-10-17/report_naples_verified.csv")
rows <- lapply(seq_len(nrow(r)), function(i) {
  x <- d[d$CA_MPA_Name_Short == "Naples SMCA" & d$source == "PISCO" &
           d$y == r$taxon[i] & d$resp == r$response[i] & d$year <= 2020 &
           is.finite(d$lnDiff), ]
  stopifnot(nrow(x) > 2L, !anyDuplicated(x$year))
  fit <- summary(lm(lnDiff ~ year, data=x))$coefficients
  data.frame(r[i,], current_n_years=nrow(x), current_first_year=min(x$year),
    current_last_year=max(x$year), current_ln_slope=fit["year","Estimate"],
    current_ln_se=fit["year","Std. Error"], current_p=fit["year","Pr(>|t|)"],
    current_log10_slope=fit["year","Estimate"]/log(10),
    caveat="Current corrected inputs and zero offset; report log base/zero treatment and exact annual source version unverified")
})
write.csv(do.call(rbind,rows),
  "data/external/source_reconciliation_2026-10-17/naples_report_style_diagnostics.csv",row.names=FALSE)
print(do.call(rbind,rows))
