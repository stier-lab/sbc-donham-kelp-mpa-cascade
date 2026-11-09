source('code/R/primary_uncertainty_utils.R')
o<-'outputs/biomass_joint_checks_20261102/'
a<-read.csv('outputs/biomass_resolution_20261102/review_annual_input.csv')
b<-read.csv(paste0(o,'baseline_models.csv'))
for(i in seq_len(nrow(b))){z<-a[a$series_id==b$series_id[i]&a$t>0,];f<-lm(lnDiff~t,z);stopifnot(abs(coef(f)[2]-b$estimate[i])<1e-9,abs(sqrt(uncertainty_hac(f,z$year,3)[2,2])-b$se[i])<1e-9)}
d<-read.csv(paste0(o,'private/verification_draws.csv'));gg<-split(d,interaction(d$series_id,d$scheme,d$draw,drop=TRUE));for(z in gg)stopifnot(abs(coef(lm(lnDiff~t,z))[2]-z$expected_slope[1])<1e-9)
# Verify independent arithmetic for KFM expected mass using original plant records.
p<-read.csv('/Users/adrianstier/sbc-donham-kelp-mpa-data/data/MBON/KFM_Macrocystis_RawData_1984-2023.csv');p<-unique(p);p$date<-as.Date(p$SurveyDate,'%m/%d/%Y');p$year<-as.integer(format(p$date,'%Y'));p$month<-as.integer(format(p$date,'%m'));p$site<-sprintf('a-k-%02d',p$SiteNumber)
c<-read.csv('data/external/sbc1275/Macrocystis_Biomass_Relationships_20261016.csv',fileEncoding='UTF-8-BOM');c<-c[c$SP_CODE=='MAPY'&c$model==1,];p$mass<-p$Stipe_Count*c$slope[match(p$month,c$MONTH)]*1000
k<-read.csv(paste0(o,'kelp_expected_biomass_check.csv'));u<-read.csv(paste0(o,'private/kfm.sum.csv'));u<-u[u$taxon_name=='Macrocystis pyrifera'&u$area==10,]
for(i in seq_len(nrow(k))){pp<-p[p$site==k$site[i]&p$year==k$year[i],];uu<-u[u$site_id==k$site[i]&u$year==k$year[i],];density<-mean(uu$count/uu$area);expected<-if(density==0)0 else density*mean(pp$mass);stopifnot(abs(expected-k$expected_biomass[i])<1e-8)}
cat('Verified',nrow(b),'baseline models,',length(gg),'draw regressions and',nrow(k),'kelp station-year expectations.\n')
