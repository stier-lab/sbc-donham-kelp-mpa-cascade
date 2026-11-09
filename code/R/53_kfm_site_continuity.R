# Fixed-site sensitivity for the two primary KFM sheephead pBACIPS histories.
# Input audit: python3 code/audit_kfm_site_continuity.py. No primary files changed.
local({
  suppressPackageStartupMessages({library(here);library(metafor)})
  root<-here::here();source(file.path(root,'code/R/primary_uncertainty_utils.R'),local=TRUE)
  status<-system2('python3',shQuote(file.path(root,'code/audit_kfm_site_continuity.py')),stdout=FALSE)
  if(status!=0L)stop('Count-site continuity audit failed; do not fit stale inputs')
  paths<-c('data/sumstats_final.csv','data/harmonized/harmonized_response_ratios.csv',
    'data/harmonized/harmonized_landsat_rr.csv','tables/table_s_kfm_continuity_annual.csv',
    'code/R/primary_uncertainty_utils.R','code/R/53_kfm_site_continuity.R','code/audit_kfm_site_continuity.py')
  hashes<-tools::md5sum(file.path(root,paths))
  primary<-read.csv(file.path(root,paths[1]));annual<-uncertainty_annual(root,primary)
  inputs<-read.csv(file.path(root,paths[4]));defs<-uncertainty_definitions(root)
  out<-list();scenarios<-list(primary=primary,matched_years=primary,persistent_sites=primary)
  for(mpa in unique(inputs$MPA)) {
    idx<-which(primary$MPA==mpa & primary$Source=='KFM' & primary$Taxa=='S. pulcher' & primary$Resp=='Den')
    stopifnot(length(idx)==1);row<-primary[idx,]
    rebuilt<-uncertainty_reconstruct(row,annual,defs)
    x<-inputs[inputs$MPA==mpa,];d<-rebuilt$data[rebuilt$data$year %in% x$year,]
    x<-x[match(d$year,x$year),];stopifnot(all(d$year==x$year),max(abs(d$lnDiff-x$lnRR_original))<1e-10)
    for(scenario in c('matched_years','persistent_sites')) {
      d$lnDiff<-if(scenario=='matched_years')x$lnRR_original else x$lnRR_persistent
      fit<-switch(row$Model,Linear=lm(lnDiff~time.model,d),Step=lm(lnDiff~BA,d),stop('Unexpected model'))
      gradient<-rebuilt$gradient
      est<-sum(gradient*coef(fit));se<-sqrt(drop(t(gradient)%*%vcov(fit)%*%gradient))
      lag<-max(1L,floor(4*(nrow(d)/100)^(2/9)))
      hac<-max(se,sqrt(drop(t(gradient)%*%uncertainty_hac(fit,d$year,lag)%*%gradient)))
      stopifnot(is.finite(est),is.finite(se),se>0)
      scenarios[[scenario]]$Mean[idx]<-est;scenarios[[scenario]]$SE[idx]<-se
      out[[length(out)+1]]<-data.frame(MPA=mpa,scenario,Model=row$Model,n_years=nrow(d),
        years=paste(d$year,collapse=';'),estimate=est,se,ci_low=est-qt(.975,df.residual(fit))*se,
        ci_high=est+qt(.975,df.residual(fit))*se,HAC_SE=hac,
        HAC_ci_low=est-qt(.975,df.residual(fit))*hac,HAC_ci_high=est+qt(.975,df.residual(fit))*hac,
        primary_estimate=row$Mean,primary_SE=row$SE,primary_n_years=row$N)
    }
  }
  contrasts<-do.call(rbind,out)
  pooled<-do.call(rbind,lapply(names(scenarios),function(s)uncertainty_joint(scenarios[[s]],scenarios[[s]]$SE,s)))
  stopifnot(identical(hashes,tools::md5sum(file.path(root,paths))))
  write.csv(contrasts,file.path(root,'tables/table_s_kfm_continuity_contrasts.csv'),row.names=FALSE)
  write.csv(pooled,file.path(root,'tables/table_s_kfm_continuity_pooled.csv'),row.names=FALSE)
  write.csv(data.frame(path=paths,md5=unname(hashes)),file.path(root,'tables/table_s_kfm_continuity_inputs.csv'),row.names=FALSE)
  print(contrasts[,c('MPA','scenario','n_years','estimate','se','ci_low','ci_high')])
  print(pooled[pooled$Taxa=='S. pulcher' & pooled$Resp=='Den',])
})
