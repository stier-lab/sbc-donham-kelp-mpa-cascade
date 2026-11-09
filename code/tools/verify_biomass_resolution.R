source('code/R/primary_uncertainty_utils.R')
a<-read.csv('outputs/pisco_fish_rules_20261101/approved_annual_input.csv')
e<-read.csv('outputs/biomass_resolution_20261102/matched_donor_effects.csv')
l<-read.csv('outputs/biomass_resolution_20261102/donor_coverage_losses.csv')
for(i in which(e$scenario=='current')){z<-a[a$series_id==e$series_id[i]&a$t>0,];f<-lm(lnDiff~t,z);stopifnot(abs(coef(f)[2]-e$estimate[i])<1e-9,abs(sqrt(uncertainty_hac(f,z$year,3)[2,2])-e$se[i])<1e-9)}
p<-read.csv('outputs/biomass_support_20261102/private/pool_members.csv')
p$mass<-ifelse(p$Taxa=='M. franciscanus',.00059*p$size^2.917,ifelse(p$Taxa=='S. purpuratus',.00059*p$size^2.870,.001352821*p$size^2.913963113))
p$key<-paste(p$series_id,p$year,p$status,sep='~')
means<-tapply(seq_len(nrow(p)),p$key,function(ii)weighted.mean(p$mass[ii],p$count[ii]))
i<-read.csv('outputs/biomass_support_20261102/animal_pool_inventory.csv');i$key<-paste(i$series_id,i$year,i$status,sep='~')
o<-read.csv('outputs/all_program_station_refits_20261102/offsets.csv');n<-0
for(flag in c('same_program','same_count_site')){
 q<-p[p[[flag]],];m<-tapply(seq_len(nrow(q)),q$key,function(ii)weighted.mean(q$mass[ii],q$count[ii]));ratio<-m/means[names(m)]
 for(k in which(e$scenario==flag)){
  z<-a[a$series_id==e$series_id[k]&a$t>0,];v<-list()
  for(st in c('mpa','reference')){key<-paste(z$series_id,z$year,st,sep='~');value<-i$value[match(key,i$key)];rr<-ratio[key];rr[value==0]<-1;v[[st]]<-value*rr}
  offset<-o$offset[match(e$series_id[k],o$series_id)];z$new<-log((v$mpa+offset)/(v$reference+offset));z<-z[is.finite(z$new),];f<-lm(new~t,z)
  stopifnot(nrow(z)==e$n_years[k],abs(coef(f)[2]-e$estimate[k])<1e-9,abs(sqrt(uncertainty_hac(f,z$year,3)[2,2])-e$se[k])<1e-9);n<-n+1
 }
}
cat('Independently verified 35 original slopes and',n,'matched donor refits.\n')
