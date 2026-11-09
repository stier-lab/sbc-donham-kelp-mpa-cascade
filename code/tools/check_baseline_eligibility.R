suppressPackageStartupMessages(library(sandwich))
source('code/R/primary_uncertainty_utils.R')
d<-read.csv('outputs/researcher_scope_20261101/annual_input.csv');ledger<-read.csv('outputs/model_eligibility_20261101/model_eligibility.csv');err<-numeric()
for(id in ledger$series_id[is.finite(ledger$slope)]) {
 x<-subset(d,series_id==id&t<0);x<-x[order(x$year),];f<-lm(lnDiff~year,x);r<-ledger[ledger$series_id==id,]
 v<-list(vcov(f),vcovHC(f,type='HC3'),uncertainty_hac(f,x$year,1),uncertainty_hac(f,x$year,3),uncertainty_hac(f,x$year,5));se<-max(vapply(v,function(z)sqrt(z[2,2]),numeric(1)))
 err<-c(err,abs(coef(f)[2]-r$slope),abs(se-r$slope_se));stopifnot(tail(err,2)<1e-7)
}
write.csv(data.frame(linear_baselines=sum(is.finite(ledger$slope)),max_error=max(err)),'outputs/model_eligibility_20261101/independent_baseline_check.csv',row.names=FALSE)
print(max(err))
