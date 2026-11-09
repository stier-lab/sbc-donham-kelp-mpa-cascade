# Decision audit only: no primary input, classification or manuscript is rewritten.
suppressPackageStartupMessages({library(data.table);library(dplyr)})
root<-here::here();read<-function(n)read.csv(file.path(root,'tables',paste0('table_s_pisco_fixed_',n,'.csv')))
out<-function(x,n)write.csv(x,file.path(root,'tables',paste0('table_s_pisco_baseline_',n,'.csv')),row.names=FALSE)
raw<-file.path(Sys.getenv('MPA_DATA_REPO',file.path(dirname(root),'sbc-donham-kelp-mpa-data')),'data/PISCO/MLPA_kelpforest_fish_2024.csv')
h<-tools::md5sum(raw);c<-read('coverage');ss<-unique(unlist(strsplit(paste(c$protected[c$eligible],c$reference[c$eligible],sep=';'),';')))
f<-fread(raw);dates<-unique(f[site %in% ss & year==2012 & campus %in%c('UCSB','VRG') & level!='CAN',.(site,year,month,day)])
dates[,survey_date:=as.Date(sprintf('%04d-%02d-%02d',year,month,day))]
dates[,effective_date:=as.Date('2012-01-01')];dates[,after_effective_date:=survey_date>=effective_date]
stopifnot(nrow(dates)==19,!anyNA(dates$survey_date),all(dates$after_effective_date),identical(h,tools::md5sum(raw)))
out(as.data.frame(dates[order(site,month,day)]),'2012_dates')
p<-read.csv(file.path(root,'data/external/pisco_policy_20261031/primary_before.csv'));p$key<-do.call(paste,c(p[c('MPA','Taxa','Source','Resp')],sep='|'))
z<-read('fits')%>%filter(scenario %in%c('strict_min3','strict_min4'),membership=='fixed',method=='pBACIPS')
b<-p%>%filter(Source=='PISCO',Taxa=='S. pulcher')%>%select(MPA,Resp,old_N=N,old_model=Model,old_estimate=Mean,old_SE=SE)
z<-left_join(z,b,by=c('MPA','Resp'));out(z,'candidate_changes')
summary<-z%>%group_by(scenario)%>%summarise(histories_reassigned=n(),old_years=sum(old_N),new_years=sum(N),new_total_annual_rows=sum(p$N)-sum(old_N)+sum(N),.groups='drop')
summary$proposed_CI<-sum(p$AnalysisType=='CI')-summary$histories_reassigned
summary$proposed_pBACIPS<-nrow(p)-summary$proposed_CI
stopifnot(summary$histories_reassigned[summary$scenario=='strict_min3']==6,summary$histories_reassigned[summary$scenario=='strict_min4']==2)
out(summary,'options')
out(data.frame(path=raw,md5=unname(h),effective_date='2012-01-01',date_source='https://wildlife.ca.gov/Conservation/Marine/MPAs/Statistics',checked='2026-10-31'),'sources')
print(summary);cat('PASS: 19 observed 2012 station-dates are after implementation; prospective changes exported without primary mutation.\n')
