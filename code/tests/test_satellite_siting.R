source("code/R/satellite_siting_utils.R")
d <- data.frame(year=1980:2020,inside=2,reference=1,contrast=log(2))
stopifnot(siting_window(d,2000,5)$change_log==0,
  siting_window(d,2000,5)$n_before==5,siting_window(d,2000,5)$n_after==5)
d$contrast[d$year>2000] <- log(4)
stopifnot(abs(siting_window(d,2000,5)$change_log-log(2))<1e-12)
stopifnot(!siting_window(d[!d$year %in% c(1995,1996),],2000,5)$eligible)
d$inside[d$year<2000] <- 0
stopifnot(!siting_window(d,2000,5)$positive_means,is.na(siting_window(d,2000,5)$change_period_log))
stopifnot(is.na(siting_tracking(rep(0,10),1:10)))
if(file.exists("tables/table_s_siting_windows.csv")) {
  x <- read.csv("tables/table_s_siting_windows.csv")
  stopifnot(all(x$anchor[x$scenario!="actual"]+x$width[x$scenario!="actual"]<x$start[x$scenario!="actual"]))
  for(w in unique(x$width)) {
    z <- x[x$width==w & x$fixed_cohort,]
    stopifnot(all(table(z$MPA)==3L))
  }
}
cat("Satellite siting tests passed\n")
