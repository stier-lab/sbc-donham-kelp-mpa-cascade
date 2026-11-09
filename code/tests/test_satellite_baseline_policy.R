source('code/R/satellite_baseline_policy.R')
x <- read.csv('data/harmonized/harmonized_landsat_rr.csv')
s <- satellite_baseline_screen(x)
stopifnot(nrow(s) == 20L,
  setequal(s$MPA[!s$baseline_eligible], c('Santa Barbara Island SMR','Swamis SMCA')),
  s$baseline_eligible[s$MPA == 'Naples SMCA'])
# Translation of calendar origin and row ordering cannot change the decision.
y <- x[nrow(x):1, ]; y$year <- y$year - 1900
s2 <- satellite_baseline_screen(y)
stopifnot(identical(s$baseline_eligible,s2$baseline_eligible),
  max(abs(s$baseline_p-s2$baseline_p)) < 1e-10)
# Missing and duplicate baseline years must not masquerade as replication.
z <- subset(x,CA_MPA_Name_Short == 'Naples SMCA')
b <- which(z$BA == 'Before'); z$lnDiff[b[-c(1,2)]] <- NA
stopifnot(!satellite_baseline_screen(z)$baseline_eligible)
stopifnot(inherits(try(satellite_baseline_screen(rbind(x,x[1,])),silent=TRUE),'try-error'))
cat('PASS: satellite baseline eligibility, calendar invariance and invalid-year handling.\n')
