# Independent R refit of every displayed diagnostic coefficient and uncertainty interval.
source('code/R/primary_uncertainty_utils.R')
O<-'outputs/all_program_station_refits_20261102';v<-read.csv(file.path(O,'verified_station_values.csv'));e<-read.csv(file.path(O,'effects.csv'));a<-read.csv('outputs/pisco_fish_rules_20261101/approved_annual_input.csv');offset<-read.csv(file.path(O,'offsets.csv'));errors<-list()
for(ii in seq_len(nrow(e))){
 q<-e[ii,];allz<-v[v$series_id==q$series_id,];z<-allz[allz$t>0&!is.na(allz$value),];off<-offset$offset[offset$series_id==q$series_id]
 if(grepl('Log_station',q$scenario)){
  z$protected<-as.integer(z$status=='mpa');z$station<-interaction(z$status,z$site,drop=TRUE);z$yr<-factor(z$year);z$pt<-z$protected*z$t;z$log_value<-log(z$value+off);z$w<-1/ave(z$value,interaction(z$year,z$status),FUN=length)
  m<-if(q$scenario=='Log_station_adjusted')lm(log_value~station+yr+pt,z,weights=w) else lm(log_value~protected+yr+pt,z,weights=w)
  X<-model.matrix(m);bread<-solve(crossprod(X,z$w*X));score<-X*as.numeric(z$w*residuals(m));ny<-length(unique(z$year));ss<-matrix(0,nrow=max(z$year)-min(z$year)+1,ncol=ncol(X))
  for(y in unique(z$year))ss[y-min(z$year)+1,]<-colSums(score[z$year==y,,drop=FALSE])
  meat<-crossprod(ss);lag<-min(3,nrow(ss)-1)
  for(h in seq_len(lag)){cc<-crossprod(ss[(h+1):nrow(ss),,drop=FALSE],ss[1:(nrow(ss)-h),,drop=FALSE]);meat<-meat+(1-h/(lag+1))*(cc+t(cc))}
  V<-bread%*%meat%*%bread*ny/(ny-1);idx<-which(names(coef(m))=='pt');estimate<-coef(m)[idx];se<-sqrt(max(0,V[idx,idx]));df<-ny-2;n<-nrow(z)
 }else{
  current<-a[a$series_id==q$series_id&a$t>0,]
  if(q$scenario=='Approved'){g<-current}else{
   r<-if(q$scenario=='Persistent_post_reference')z[z$status=='reference',] else allz[allz$status=='reference',]
   sets<-lapply(split(r,r$year),function(x)x$site[!is.na(x$value)]);fixed<-Reduce(intersect,sets)
   zz<-z[z$status=='mpa'|z$site %in% fixed,];means<-aggregate(value~year+status,zz,mean);g<-reshape(means,direction='wide',idvar='year',timevar='status')
   stopifnot(all(c('value.mpa','value.reference') %in% names(g)));g<-g[complete.cases(g),];g$lnDiff<-log((g$value.mpa+off)/(g$value.reference+off));g$t<-g$year-z$MPA_Start[1]
   if(q$scenario=='Available_on_persistent_years')g<-current[current$year %in% g$year,]
  }
  g<-g[order(g$year),];m<-lm(lnDiff~t,g);V<-uncertainty_hac(m,g$year,3);estimate<-coef(m)[2];se<-sqrt(max(0,V[2,2]));df<-df.residual(m);ny<-nrow(g);n<-ny
 }
 lo<-estimate-qt(.975,df)*se;hi<-estimate+qt(.975,df)*se;err<-max(abs(c(estimate-q$estimate,se-q$se,lo-q$low,hi-q$high)));stopifnot(is.finite(err),err<1e-7,ny==q$n_years,n==q$n_rows,df==q$df)
 errors[[ii]]<-data.frame(series_id=q$series_id,scenario=q$scenario,max_error=err)
}
x<-do.call(rbind,errors);write.csv(x,'outputs/dataset_inventory_20261101/model_accuracy_checks.csv',row.names=FALSE);print(c(fits=nrow(x),max_error=max(x$max_error)))
