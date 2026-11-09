# Conditional treatments of five antenna-only LTER lobster observations.
# No primary promotion. Keep all primary histories/years and the retained model structures.
local({
 suppressPackageStartupMessages({library(data.table);library(dplyr);library(metafor)})
 root<-normalizePath('.');dr<-Sys.getenv('DONHAM_DATA_REPO',file.path(dirname(root),'sbc-donham-kelp-mpa-data'))
 source('code/R/primary_uncertainty_utils.R',local=TRUE);source('code/R/temporal_uncertainty_utils.R',local=TRUE)
 protected<-c('data/sumstats_final.csv','data/harmonized/harmonized_response_ratios.csv','data/harmonized/harmonized_raw_responses.csv','tables/table_02_meta_analysis.csv','tables/table_03_cross_taxa_meta_regression.csv','data/cache/temporal_primary_models.rds');hash<-tools::md5sum(protected)
 primary<-read.csv(protected[1]);primary$key<-uncertainty_key(primary);annual<-uncertainty_annual(root,primary);defs<-uncertainty_definitions(root)
 selected<-primary[primary$Taxa=='P. interruptus' & primary$Source=='LTER' & primary$Resp=='Bio',];stopifnot(nrow(selected)==2)
 primary_scenario<-if(file.exists('data/harmonized/lobster_missing_lengths.json'))'site_year_mean_mass' else 'current_zero'
 ids<-selected$key
 baselinefits<-setNames(lapply(seq_len(nrow(selected)),function(i)uncertainty_reconstruct(selected[i,],annual,defs)),ids)
 lobpath<-file.path(dr,'data/LTER/Lobster_Abundance_All_Years_20230831.csv');geopath<-file.path(dr,'data/MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv')
 d<-fread(lobpath);d[,source_row:=.I];geo<-fread(geopath)
 d<-merge(d,geo[,.(SITE=site_id,MPA=CA_MPA_Name_Short,status)],by='SITE',allow.cartesian=TRUE)[MPA %in% selected$MPA]
 d[,missing:=COUNT>0 & (!is.finite(SIZE_MM)|SIZE_MM<=0)]
 stopifnot(sum(d[missing==TRUE,COUNT])==5,all(d[missing==TRUE,NUM_AO]==d[missing==TRUE,COUNT]))
 mass<-function(x)0.001352821*x^2.913963113
 donors<-d[COUNT>0 & SIZE_MM>0,.(n_measured=sum(COUNT),mean_mass=sum(COUNT*mass(SIZE_MM))/sum(COUNT),min_length=min(SIZE_MM),max_length=max(SIZE_MM),min_mass=min(mass(SIZE_MM)),max_mass=max(mass(SIZE_MM))),by=.(SITE,YEAR)]
 td<-d[COUNT>0 & SIZE_MM>0,.(transect_mean=sum(COUNT*mass(SIZE_MM))/sum(COUNT)),by=.(SITE,YEAR,TRANSECT)]
 d<-merge(d,donors,by=c('SITE','YEAR'),all.x=TRUE);d<-merge(d,td,by=c('SITE','YEAR','TRANSECT'),all.x=TRUE)
 stopifnot(all(is.finite(d[missing==TRUE,mean_mass])),all(is.finite(d[missing==TRUE,transect_mean])))
 summarize<-function(z){z[,.(b=sum(biomass)),by=.(MPA,YEAR,status,SITE,TRANSECT,REPLICATE)][,.(b=mean(b)),by=.(MPA,YEAR,status,SITE,TRANSECT)][,.(b=mean(b)),by=.(MPA,YEAR,status,SITE)][,.(value=mean(b)/300),by=.(MPA,year=YEAR,status)]}
 zero_d<-copy(d);zero_d[,biomass:=0];zero_d[COUNT>0 & SIZE_MM>0,biomass:=COUNT*mass(SIZE_MM)]
 zero_annual<-summarize(zero_d);setnames(zero_annual,'value','zero_contribution_value')
 observed<-fread(protected[3])[source=='LTER' & resp=='Bio' & taxon_name=='Panulirus interruptus'];setnames(observed,'CA_MPA_Name_Short','MPA')
 scenarios<-c('current_zero','site_year_mean_mass','transect_year_mean_mass','site_year_minimum','site_year_maximum','exclude_affected_replicates')
 results<-list();pooled<-list();pathways<-list();temporal<-list();annuals<-list();series<-list();fitchecks<-list()
 tc<-readRDS(protected[6]);tdata<-tc$annual[tc$annual$resp=='Bio',]
 stopifnot(nrow(tdata)==nobs(tc$bio));full<-temporal_primary_data(read.csv(protected[2]),primary)
 stopifnot(identical(as.character(full$history),as.character(tc$annual$history)),identical(full$year,tc$annual$year),identical(full$lnDiff,tc$annual$lnDiff))
 for(sc in scenarios){
  z<-copy(d);z[,biomass:=0];z[COUNT>0 & SIZE_MM>0,biomass:=COUNT*mass(SIZE_MM)]
  if(sc=='site_year_mean_mass')z[missing==TRUE,biomass:=COUNT*mean_mass]
  if(sc=='transect_year_mean_mass')z[missing==TRUE,biomass:=COUNT*transect_mean]
  if(sc=='site_year_minimum')z[missing==TRUE,biomass:=COUNT*min_mass]
  if(sc=='site_year_maximum')z[missing==TRUE,biomass:=COUNT*max_mass]
  if(sc=='exclude_affected_replicates'){bad<-unique(z[missing==TRUE,.(SITE,YEAR,TRANSECT,REPLICATE)]);z<-z[!bad,on=.(SITE,YEAR,TRANSECT,REPLICATE)]}
  a<-summarize(z);baseline<-merge(a,observed[,.(MPA,year,status,baseline=value)],by=c('MPA','year','status'));stopifnot(nrow(baseline)==48)
  if(sc==primary_scenario)stopifnot(max(abs(baseline$value-baseline$baseline))<1e-10)
  baseline<-merge(baseline,zero_annual,by=c('MPA','year','status'));baseline[,`:=`(scenario=sc,pct_change=100*(value/zero_contribution_value-1),pct_change_from_current=100*(value/baseline-1))];annuals[[sc]]<-baseline
  # Reproduce adaptive correction on each complete MPA history, including zeros.
  a[,corrected:=value+min(value[value>0])/2,by=MPA];rr<-dcast(a,MPA+year~status,value.var='corrected');rr[,lnDiff:=log(mpa/reference)]
  x<-primary;ta<-tdata
  for(i in seq_len(nrow(selected))){row<-selected[i,];r<-baselinefits[[row$key]];nd<-r$data;j<-match(nd$year,rr[MPA==row$MPA,year]);stopifnot(!anyNA(j));nd$lnDiff<-rr[MPA==row$MPA,lnDiff][j]
   if(sc==primary_scenario)stopifnot(max(abs(nd$lnDiff-r$data$lnDiff))<1e-12)
   f<-if(row$Model=='Mean')lm(lnDiff~1,nd) else lm(lnDiff~time.model,nd)
   g<-r$gradient;est<-sum(g*coef(f));se<-sqrt(drop(t(g)%*%vcov(f)%*%g));ii<-match(row$key,x$key);x$Mean[ii]<-est;x$SE[ii]<-se
   series[[paste(sc,row$key)]]<-data.frame(scenario=sc,key=row$key,estimate=est,se,baseline=row$Mean,baseline_se=row$SE)
   ix<-which(ta$history==row$key);ta$lnDiff[ix]<-nd$lnDiff[match(ta$year[ix],nd$year)]
  }
  stopifnot(identical(x$key,primary$key),identical(x$N,primary$N))
  pooled[[sc]]<-uncertainty_joint(x,x$SE,sc)
  mpa<-x %>% group_by(MPA,Taxa,Resp)%>%summarise(Mean=mean(Mean),vi=sum(SE^2)/n()^2,.groups='drop')
  px<-filter(mpa,Taxa=='P. interruptus',Resp=='Bio');py<-filter(mpa,Taxa=='S. purpuratus',Resp=='Bio');pd<-merge(px,py,by='MPA',suffixes=c('_x','_y'));stopifnot(nrow(pd)==4)
  fit<-rma(yi=Mean_y,vi=vi_y,mods=~Mean_x,data=pd,method='REML',test='knha')
  perms<-as.matrix(expand.grid(rep(list(1:4),4)));perms<-perms[apply(perms,1,function(q)length(unique(q))==4),]
  perm_slopes<-apply(perms,1,function(q)as.numeric(coef(rma(yi=pd$Mean_y,vi=pd$vi_y,mods=~I(pd$Mean_x[q]),method='REML',test='knha'))[2]))
  pathways[[sc]]<-data.frame(scenario=sc,k=4,slope=coef(fit)[2],ci_low=fit$ci.lb[2],ci_high=fit$ci.ub[2],p=fit$pval[2],exact_p=mean(abs(perm_slopes)>=abs(coef(fit)[2])-1e-12))
  tm<-if(sc==primary_scenario)tc$bio else lme4::refit(tc$bio,ta$lnDiff)
  fitchecks[[sc]]<-data.frame(scenario=sc,n=nobs(tm),singular=lme4::isSingular(tm),convergence=paste(tm@optinfo$conv$lme4$messages,collapse=';'))
  stopifnot(!lme4::isSingular(tm),is.null(tm@optinfo$conv$lme4$messages))
  temporal[[sc]]<-do.call(rbind,lapply(levels(tm@frame$Species),function(sp)cbind(scenario=sc,Species=sp,temporal_contrast(tm,sp,slope=TRUE))))
  message('Completed ',sc)
 }
 out<-function(x,n)fwrite(rbindlist(x,fill=TRUE),paste0('tables/table_s_lobster_missing_',n,'.csv'))
 out(annuals,'annual');out(series,'contrasts');out(pooled,'pooled');out(pathways,'pathway');out(temporal,'temporal');out(fitchecks,'fit_checks')
 base<-pooled[[primary_scenario]];table2<-read.csv(protected[4]);table2$Resp<-ifelse(table2$Response=='Biomass','Bio','Den');check<-merge(base,table2,by=c('Taxa','Resp'));stopifnot(nrow(check)==9,max(abs(check$estimate-check$Estimate))<1e-7,max(abs(check$q-check$pval_fdr))<1e-7)
 pt<-read.csv(protected[5]);pt<-pt[pt$Predictor=='P. interruptus biomass' & pt$Response_Var=='S. purpuratus biomass',];stopifnot(abs(round(pathways[[primary_scenario]]$slope,4)-pt$Slope_Est)<1e-8,abs(pathways[[primary_scenario]]$p-pt$Slope_pval)<1e-8)
 private<-'local_archive/lobster_missing_lengths_2026-10-30';dir.create(private,recursive=TRUE,showWarnings=FALSE)
 fwrite(d[missing==TRUE,.(MPA,SITE,YEAR,DATE,TRANSECT,REPLICATE,COUNT,NUM_AO,SIZE_MM,n_measured,mean_mass,transect_mean,min_length,max_length)],file.path(private,'source_records_and_donors.csv'))
 stopifnot(identical(hash,tools::md5sum(protected)))
 inputs<-c(protected,lobpath,geopath,'code/R/56_lobster_missing_lengths.R','code/R/primary_uncertainty_utils.R','code/R/temporal_uncertainty_utils.R');fwrite(data.table(path=inputs,md5=unname(tools::md5sum(inputs))),'tables/table_s_lobster_missing_inputs.csv')
 cat('All scenarios verified against current primary; no primary files changed.\n')
})
