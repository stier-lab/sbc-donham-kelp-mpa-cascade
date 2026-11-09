source('code/R/primary_uncertainty_utils.R')
O<-'outputs/all_program_station_refits_20261102';v<-read.csv(file.path(O,'verified_station_values.csv'));a<-read.csv('outputs/pisco_fish_rules_20261101/approved_annual_input.csv');o<-read.csv(file.path(O,'offsets.csv'));d<-read.csv(file.path(O,'reference_omission.csv'));err<-c()
for(i in seq_len(nrow(d))){
 q<-d[i,];z<-v[v$series_id==q$series_id & v$t>0 & !is.na(v$value),];z<-z[!(z$status=='reference'&z$site==q$omitted_reference),];yy<-as.integer(strsplit(q$years,';',fixed=TRUE)[[1]]);z<-z[z$year %in% yy,];g<-aggregate(value~year+status,z,mean);g<-reshape(g,direction='wide',idvar='year',timevar='status');g<-g[order(g$year),];off<-o$offset[o$series_id==q$series_id];g$lnDiff<-log((g$value.mpa+off)/(g$value.reference+off));g$t<-g$year-z$MPA_Start[1];m<-lm(lnDiff~t,g);se<-sqrt(uncertainty_hac(m,g$year,3)[2,2]);lo<-coef(m)[2]-qt(.975,df.residual(m))*se;hi<-coef(m)[2]+qt(.975,df.residual(m))*se
 x<-max(abs(c(coef(m)[2]-q$after,lo-q$after_low,hi-q$after_high)));stopifnot(is.finite(x),x<1e-7);err<-c(err,x)
}
write.csv(data.frame(omission_fits_checked=nrow(d),max_coefficient_interval_error=max(err)),file.path(O,'independent_omission_check.csv'),row.names=FALSE);print(c(checked=nrow(d),max_error=max(err)))
