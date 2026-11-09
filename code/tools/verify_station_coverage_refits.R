d<-read.csv('outputs/station_coverage_refits_20261101/station_model_rows.csv');e<-read.csv('outputs/station_coverage_refits_20261101/effects.csv');errors<-c();n<-0
for(id in unique(d$series_id)){
 z<-d[d$series_id==id,];z$protected<-as.integer(z$status=='mpa');z$station<-interaction(z$status,z$site,drop=TRUE);z$yr<-factor(z$year);z$pt<-z$protected*z$t
 for(sc in c('Log_station_adjusted','Log_station_unadjusted')){
  m<-if(sc=='Log_station_adjusted')lm(log_value~station+yr+pt,z,weights=weight) else lm(log_value~protected+yr+pt,z,weights=weight)
  target<-e[e$series_id==id&e$scenario==sc,];if(!nrow(target))next
  err<-abs(coef(m)['pt']-target$estimate);stopifnot(is.finite(err),err<1e-8);errors<-c(errors,err);n<-n+1
 }
}
write.csv(data.frame(models_verified=n,max_coefficient_error=max(errors)), 'outputs/station_coverage_refits_20261101/independent_R_check.csv',row.names=FALSE);print(c(models=n,max_error=max(errors)))
