source('code/R/primary_uncertainty_utils.R')
o<-'outputs/pisco_biomass_borrowing_20261102/'
s<-read.csv('outputs/pisco_timing_donors_20261102/size_sample_inventory.csv');s$date<-as.Date(s$date)
c<-read.csv(paste0(o,'outer_site_calibration.csv'));c$cal_date<-as.Date(c$cal_date)
blocks<-split(s,paste(s$series_id,s$year,sep='~'))
for(i in seq_len(nrow(c))){t<-c[i,];d<-blocks[[paste(t$series_id,t$cal_year,sep='~')]];d<-d[d$status==t$status&d$zone==t$cal_zone&d$site!=t$outer_site,]
 if(t$kind=='other_site')d<-d[d$site!=t$cal_site,] else d<-d[d$site==t$cal_site&d$date!=t$cal_date,]
 d$gap<-abs(as.numeric(d$date-t$cal_date));d<-d[d$gap==ave(d$gap,d$site,FUN=min),];dd<-aggregate(cbind(total_mass,n_sized)~site,d,sum);pred<-mean(dd$total_mass/dd$n_sized);stopifnot(abs(pred-t$predicted)<1e-8,abs(log(t$observed/pred)-t$residual)<1e-9)
}
cv<-read.csv(paste0(o,'heldout_interval_validation.csv'));cv<-cv[tolower(as.character(cv$available))=='true',];calblocks<-split(c,paste(c$series_id,c$outer_site,c$status,c$kind,sep='~'))
for(i in seq_len(nrow(cv))){t<-cv[i,];z<-calblocks[[paste(t$series_id,t$site,t$status,t$kind,sep='~')]];z$weight<-1/ave(rep(1,nrow(z)),z$cluster,FUN=sum);z<-z[order(z$residual),];qq<-(cumsum(z$weight)-.5*z$weight)/sum(z$weight);q<-approx(qq,z$residual,xout=c(.025,.5,.975),rule=2)$y;stopifnot(max(abs(t$predicted*exp(q)-c(t$low,t$point,t$high)))<1e-8)}
a<-read.csv(paste0(o,'biomass_annual.csv'));e<-read.csv(paste0(o,'biomass_effects.csv'))
for(i in seq_len(nrow(e))){z<-a[a$series_id==e$series_id[i]&a$t>0&is.finite(a$biomass_lnratio),];f<-lm(biomass_lnratio~t,z);stopifnot(abs(coef(f)[2]-e$estimate[i])<1e-9,abs(sqrt(uncertainty_hac(f,z$year,3)[2,2])-e$se[i])<1e-9)}
d<-read.csv(paste0(o,'private/verification_draws.csv'));g<-split(d,interaction(d$series_id,d$draw,drop=TRUE));for(z in g)stopifnot(abs(coef(lm(lnDiff~t,z))[2]-z$expected_slope[1])<1e-9)
cat('Verified',nrow(c),'outer calibration predictions,',nrow(cv),'held-out intervals,',nrow(e),'biomass models and',length(g),'draw regressions.\n')
