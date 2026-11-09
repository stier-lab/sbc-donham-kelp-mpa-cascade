source('code/R/primary_uncertainty_utils.R')
o<-'outputs/pisco_timing_donors_20261102/'
e<-read.csv(paste0(o,'timing_effects.csv'));a<-read.csv(paste0(o,'timing_annual.csv'));a$complete<-tolower(as.character(a$complete))=='true'
for(i in seq_len(nrow(e))){z<-a[a$series_id==e$series_id[i]&a$tolerance_days==e$tolerance_days[i]&a$t>0,];if(e$comparison[i]!='baseline_all')z<-z[z$complete,];z$y<-z[[if(e$comparison[i]=='date_paired')'paired' else 'baseline']];f<-lm(y~t,z);stopifnot(nrow(z)==e$n_years[i],abs(coef(f)[2]-e$estimate[i])<1e-9,abs(sqrt(uncertainty_hac(f,z$year,3)[2,2])-e$se[i])<1e-9)}
p<-read.csv(paste0(o,'donor_predictions.csv'));s<-read.csv(paste0(o,'size_sample_inventory.csv'));s$date<-as.Date(s$date);p$date<-as.Date(p$date);blocks<-split(s,paste(s$series_id,s$year,sep='~'))
for(i in seq_len(nrow(p))){t<-p[i,];d<-blocks[[paste(t$series_id,t$year,sep='~')]];k<-t$strategy
 if(k=='same_site_other_visit')d<-d[d$site==t$site&d$zone==t$zone&d$date!=t$date,]
 if(k=='same_site_other_zone')d<-d[d$site==t$site&d$zone!=t$zone&d$date==t$date,]
 if(k=='same_role_other_site')d<-d[d$site!=t$site&d$status==t$status&d$zone==t$zone,]
 if(k=='other_role_other_site')d<-d[d$site!=t$site&d$status!=t$status&d$zone==t$zone,]
 if(k=='same_role_other_site_any_zone')d<-d[d$site!=t$site&d$status==t$status,]
 d$gap<-abs(as.numeric(d$date-t$date));d<-d[d$gap<=t$tolerance_days,];d<-d[d$gap==ave(d$gap,d$site,FUN=min),];dd<-aggregate(cbind(total_mass,n_sized)~site,d,sum);pred<-mean(dd$total_mass/dd$n_sized)
 stopifnot(abs(pred-t$predicted)<1e-8,abs(log(pred/t$observed)-t$log_error)<1e-9)
}
cat('Verified',nrow(e),'timing models and',nrow(p),'held-out donor predictions independently in R.\n')
