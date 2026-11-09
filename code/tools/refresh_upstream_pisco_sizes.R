Sys.setenv(DONHAM_KELP_CALIBRATION_CSV=normalizePath('../../../data/external/sbc1275/Macrocystis_Biomass_Relationships_20261016.csv'))
suppressPackageStartupMessages({library(here);library(dplyr)})
here::i_am('run_review.R');source('code/R/00c_analysis_constants.R');source('code/R/01_utils.R')
SizeFreq<-readRDS('../kfm_urchin_size_pool.rds');SizeFreq<-subset(SizeFreq,size>=PISCO_URCHIN_MIN_SIZE_MM)
u<-readRDS('data/cache/pisco_urchin_bootstrap.rds');u<-subset(u,CA_MPA_Name_Short %in% c('Anacapa Island SMCA','Painted Cave SMCA','Point Dume SMR','Point Dume SMCA'))
countcheck<-u%>%group_by(CA_MPA_Name_Short,site,year,y)%>%summarise(rebuilt=mean(count/transect)/60,.groups='drop')
s<-read.csv('../pisco_site_values.csv');fresh<-s%>%filter(y%in%c('MESFRAAD','STRPURAD'))%>%select(CA_MPA_Name_Short,site,year,y,density)%>%distinct()
fresh<-fresh%>%filter(CA_MPA_Name_Short%in%unique(u$CA_MPA_Name_Short))
v<-inner_join(countcheck,fresh,by=c('CA_MPA_Name_Short','site','year','y'));stopifnot(nrow(v)==nrow(fresh),nrow(v)>0,max(abs(v$rebuilt-v$density))<1e-8)
write.csv(v,'../pisco_urchin_count_validation.csv',row.names=FALSE)
set.seed(SEED_PISCO_URCHIN);u$size_rows<-0L
for(i in seq_len(nrow(u))){
 t<-which(SizeFreq$CA_MPA_Name_Short==u$CA_MPA_Name_Short[i]&SizeFreq$site_status==u$site_status[i]&SizeFreq$year==u$year[i]&SizeFreq$classcode==u$y[i]);u$size_rows[i]<-length(t)
 b<-bootstrap_biomass(u$count[i],t,SizeFreq,if(u$y[i]=='MESFRAAD')bio_redurch else bio_purpurch);u$biomass[i]<-b$biomass
}
bad<-u%>%filter(count>0&!is.finite(biomass))%>%distinct(site,year,y)
u<-u%>%anti_join(bad,by=c('site','year','y'))
z<-u%>%group_by(CA_MPA_Name_Short,site,year,y)%>%summarise(biomass=mean(biomass/transect),size_rows=min(size_rows),.groups='drop')
write.csv(z,'../pisco_refreshed_urchin_biomass.csv',row.names=FALSE);write.csv(u,'../pisco_refreshed_urchin_draws.csv',row.names=FALSE);write.csv(bad,'../pisco_unsized_urchin_siteyears.csv',row.names=FALSE)
cat('Validated',nrow(v),'source count cells; rebuilt',nrow(z),'site/year biomasses\n')
