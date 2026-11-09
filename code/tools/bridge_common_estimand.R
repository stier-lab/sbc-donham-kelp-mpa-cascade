#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(metafor);library(clubSandwich)})
source('code/R/primary_uncertainty_utils.R')
out<-'outputs/pooled_common_estimand_20261101'
primary<-read.csv('data/sumstats_final.csv')
z<-uncertainty_joint(primary,primary$SE,'Reproduced_current')
old<-read.csv('tables/table_02_meta_analysis.csv');old$Resp<-ifelse(old$Response=='Biomass','Bio','Den')
j<-merge(z,old,by=c('Taxa','Resp'))
stopifnot(nrow(j)==9,max(abs(j$estimate-j$Estimate))<1e-6,max(abs(j$q-j$pval_fdr))<1e-6)
write.csv(z,file.path(out,'reproduced_current_pooled.csv'),row.names=FALSE)
# Evaluate only the pool function, not the analysis script's top-level code.
for(expr in as.list(parse('code/tools/refit_common_estimand.R'))) if(is.call(expr)&&identical(expr[[1]],as.name('<-'))&&identical(expr[[2]],as.name('pool')))eval(expr)
cl<-read.csv(file.path(out,'reference_clusters.csv'));logs<-list();e<-read.csv(file.path(out,'effects.csv'));rows<-list()
for(cohort in c('All_current','Change_eligible','Post_eligible')) for(resp in c('Bio','Den')) for(taxon in unique(primary$Taxa)) {
 x<-primary[primary$Resp==resp & primary$Taxa==taxon,];x$series_id<-uncertainty_key(x)
 if(cohort!='All_current') {
  metric<-if(cohort=='Change_eligible')'Change' else 'Post_difference'
  keys<-e$series_id[e$scenario=='Primary'&e$metric==metric]
  x<-x[x$series_id%in%keys,]
 }
 if(!nrow(x))next
 x<-merge(x,cl,by='MPA');x$estimate<-x$Mean;x$vi<-x$SE^2;x$scenario<-cohort;x$metric<-'Legacy_mixture'
 rows[[length(rows)+1]]<-pool(x,'All')
}
p<-do.call(rbind,rows);p$q<-ave(p$p,p$scenario,FUN=function(x)p.adjust(x,'BH',n=9))
write.csv(p,file.path(out,'legacy_bridge.csv'),row.names=FALSE)
print(p[p$Taxa=='M. pyrifera',c('scenario','estimate','low','high','k','q')],row.names=FALSE)
cat('Current nine-response Table2 reproduced. Legacy effects refit with new pooling/CR2 without changing effect definitions.\n')
# Exact same histories: isolate redefining effects from satellite-cohort expansion.
new<-subset(e,scenario=='Primary'&metric=='Change'&Taxa=='M. pyrifera'&Resp=='Bio')
new<-new[new$series_id%in%uncertainty_key(primary),]
new$scenario<-'New_change_same_legacy_histories'
b<-pool(new,'All');write.csv(b,file.path(out,'kelp_same_cohort_bridge.csv'),row.names=FALSE)
print(b[c('scenario','k','estimate','low','high')],row.names=FALSE)
