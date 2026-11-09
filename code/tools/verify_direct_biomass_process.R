# Independent conditional GLS/REML checks; does not run or modify primary analyses.
root <- 'outputs/biomass_next_steps_20261102/uncertainty'
a <- read.csv(file.path(root,'direct_annual_observations.csv'))
covtab <- read.csv(file.path(root,'direct_sampling_covariance.csv'))
e <- read.csv(file.path(root,'direct_trend_scenarios.csv'))
checks <- list()
for (i in seq_len(nrow(e))) {
 r<-e[i,];z<-a[a$series_id==r$series_id,];z<-z[order(z$year),];n<-nrow(z);X<-cbind(1,z$year-mean(z$year));y<-z$observed_log_ratio
 S<-matrix(0,n,n);cv<-covtab[covtab$series_id==r$series_id,]
 for (j in seq_len(nrow(cv))) S[match(cv$year[j],z$year),match(cv$other_year[j],z$year)]<-cv$sampling_covariance[j]
 K<-r$rho^abs(outer(z$year,z$year,'-'))
 calc<-function(q,ret=FALSE){V<-S+q*K+diag(1e-9,n);iv<-solve(V);bc<-solve(t(X)%*%iv%*%X);b<-bc%*%t(X)%*%iv%*%y;res<-y-X%*%b
   if(ret)return(list(b=b,se=sqrt(bc[2,2])))
   as.numeric(determinant(V,logarithm=TRUE)$modulus+determinant(t(X)%*%iv%*%X,logarithm=TRUE)$modulus+t(res)%*%iv%*%res)
 }
 q<-optimize(calc,c(0,max(1,10*mean((y-mean(y))^2))),tol=1e-10)$minimum;if(calc(0)<=calc(q))q<-0
 f<-calc(q,TRUE);stopifnot(abs(f$b[2]-r$estimate)<1e-6,abs(f$se-r$conditional_gls_se)<1e-6,abs(q-r$process_variance)<1e-5)
 checks[[i]]<-data.frame(series_id=r$series_id,rho=r$rho,slope_error=abs(f$b[2]-r$estimate),se_error=abs(f$se-r$conditional_gls_se),process_variance_error=abs(q-r$process_variance))
}
write.csv(do.call(rbind,checks),file.path(root,'independent_R_checks.csv'),row.names=FALSE)
cat('Independent R checks passed for',nrow(e),'direct-only GLS/REML fits.\n')
