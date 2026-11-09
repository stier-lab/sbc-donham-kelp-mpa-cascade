suppressPackageStartupMessages({library(metafor);library(clubSandwich);library(sandwich)})
source('code/R/primary_uncertainty_utils.R')
out<-'outputs/change_attribution_20261101';dir.create(out,showWarnings=FALSE)
old<-read.csv('data/sumstats_final.csv');old$series_id<-uncertainty_key(old)
cl<-read.csv('outputs/pooled_common_estimand_20261101/reference_clusters.csv')
raw<-uncertainty_annual(getwd(),old);defs<-uncertainty_definitions(getwd())
up<-read.csv('outputs/upstream_sites_20261101/annual_input.csv');up$baseline_design_verified<-as.logical(up$baseline_design_verified);stopifnot(!anyNA(up$baseline_design_verified))
proposed<-read.csv('outputs/researcher_scope_20261101/effects.csv');proposed<-subset(proposed,scenario=='Primary')
old$cluster<-cl$cluster[match(old$MPA,cl$MPA)];stopifnot(!anyNA(old$cluster))
res<-list();effectlog<-list();obslog<-list();warnlog<-list();checks<-list()
record<-function(x,stage,metric,mode='cr2') {
 if(!nrow(x))return()
 x$stage<-stage;x$metric<-metric;effectlog[[length(effectlog)+1]]<<-x
 for(resp in c('Bio','Den'))for(tx in unique(old$Taxa)) {
  z<-subset(x,Resp==resp&Taxa==tx);if(!nrow(z))next
  k<-nrow(z);ng<-length(unique(z$cluster));nm<-length(unique(z$MPA));est<-z$estimate[1];se<-df<-lo<-hi<-pv<-NA_real_;status<-'Insufficient independent information';sig<-NA_character_
  if(k>=2&&nm>=2){
   random<-if(anyDuplicated(z$MPA))list(~1|MPA,~1|series_id) else list(~1|MPA)
   m<-withCallingHandlers(rma.mv(estimate,V=vi,random=random,data=z,method='REML',test='t'),warning=function(w){warnlog[[length(warnlog)+1]]<<-data.frame(stage,metric,tx,resp,message=conditionMessage(w));invokeRestart('muffleWarning')})
   est<-unname(coef(m)[1]);sig<-paste(m$sigma2,collapse=';')
   if(mode=='model'){se<-m$se[1];df<-m$ddf[1];lo<-m$ci.lb[1];hi<-m$ci.ub[1];pv<-m$pval[1];status<-'Conventional model interval'} else {
    ct<-tryCatch(coef_test(m,vcov='CR2',cluster=z$cluster,test='Satterthwaite'),error=function(e)NULL)
    if(!is.null(ct)){se<-ct$SE[1];df<-ct$df_Satt[1];if(ng>=3&&df>=2){lo<-est-qt(.975,df)*se;hi<-est+qt(.975,df)*se;pv<-ct$p_Satt[1];status<-if(df<4)'CR2; df<4' else 'CR2'}}
   }
  }
  res[[length(res)+1]]<<-data.frame(stage,metric,Taxa=tx,Resp=resp,k,n_mpa=nm,n_clusters=ng,estimate=est,se,df,low=lo,high=hi,p=pv,status,sigma=sig)
 }
}
legacy<-data.frame(old[c('series_id','MPA','Taxa','Resp','Source','cluster')],estimate=old$Mean,vi=old$SE^2)
# Exact original synthesis; keep its output separate because joint heterogeneity is shared across taxa.
joint<-function(se,stage){z<-uncertainty_joint(old,se,stage);z$stage<-stage;z$metric<-'Legacy';write.csv(z,file.path(out,paste0(stage,'.csv')),row.names=FALSE);z}
j0<-joint(old$SE,'A0_current');tab<-read.csv('tables/table_02_meta_analysis.csv');tab$Resp<-ifelse(tab$Response=='Biomass','Bio','Den');jj<-merge(j0,tab,by=c('Taxa','Resp'));stopifnot(max(abs(jj$estimate-jj$Estimate))<1e-7,max(abs(jj$q-jj$pval_fdr))<1e-7)
record(legacy,'A1_pool_structure','Legacy','model');record(legacy,'A2_cluster_uncertainty','Legacy')
hac<-numeric(nrow(old));rebuilt<-list()
for(i in seq_len(nrow(old))){r<-uncertainty_reconstruct(old[i,],raw,defs);rebuilt[[old$series_id[i]]]<-r$data;hac[i]<-sqrt(drop(t(r$gradient)%*%uncertainty_hac(r$fit,r$data$year,3)%*%r$gradient))}
legacyh<-legacy;legacyh$vi<-hac^2;record(legacyh,'A3_annual_HAC','Legacy')
joint(pmax(hac,old$SE),'A0b_current_pool_HAC_floor')
# Fit common estimands, logging exact observations and both covariance estimates.
common<-function(d,id,metric,period='original',window=FALSE,postonly=FALSE,tag='') {
 r<-old[match(id,old$series_id),];start<-if('MPA_Start'%in%names(d))d$MPA_Start[1] else r$MPA_Start
 d<-d[order(d$year),];d$t<-d$year-start
 if(period=='strict')d<-d[d$t!=0,]
 d$post<-if(period=='original')as.integer(d$BA=='After') else as.integer(d$t>0)
 if(window)d<-d[d$post==0|d$t<=9,]
 if(postonly)d<-d[d$post==1,]
 nb<-sum(d$post==0);na<-sum(d$post==1)
 if(na<3||(metric=='Change'&&nb<3))return(NULL)
 d$pt<-d$post*d$t
 single_period<-postonly || nb==0
 f<-if(single_period)lm(lnDiff~t,d) else lm(lnDiff~post+pt,d)
 if(f$rank<length(coef(f)))return(NULL)
 X<-model.matrix(f)
 g<-if(single_period)colMeans(X) else if(metric=='Change')c(0,1,mean(d$t[d$post==1])) else c(1,1,mean(d$t[d$post==1]))
 est<-sum(g*coef(f));iid<-drop(g%*%vcov(f)%*%g);hv<-drop(g%*%uncertainty_hac(f,d$year,3)%*%g)
 direct<-mean(d$lnDiff[d$post==1])-if(metric=='Change')mean(d$lnDiff[d$post==0]) else 0
 stopifnot(abs(est-direct)<1e-8)
 if(!is.finite(hv)||hv<1e-24)return(NULL)
 obslog[[length(obslog)+1]]<<-data.frame(stage=tag,metric,series_id=id,year=d$year,lnDiff=d$lnDiff,post=d$post,t=d$t)
 meta<-if(id%in%old$series_id)old[match(id,old$series_id),c('series_id','MPA','Taxa','Resp','Source','cluster')] else up[match(id,up$series_id),c('series_id','MPA','Taxa','Resp','Source','cluster')]
 data.frame(meta,estimate=est,vi=hv,iid_vi=iid,n=nrow(d),nb,na)
}
for(metric in c('Change','Post_difference')) {
 # Common cohort eligible in both the original data and final standardized window.
 ids<-intersect(old$series_id,proposed$series_id[proposed$metric==metric])
 # Include Gull for explicit researcher-removal contrasts later; membership here is intersection with final roster.
 fitset<-function(period='original',window=FALSE,postonly=FALSE,tag='')do.call(rbind,lapply(ids,function(id)common(rebuilt[[id]],id,metric,period,window,postonly,tag)))
 record(legacyh[legacyh$series_id%in%ids,],'B0_matched_legacy',metric)
 b1<-fitset(tag='B1_common_definition_same_observations');record(b1,'B1_common_definition_same_observations',metric)
 b2<-fitset(period='strict',tag='B2_calendar_designation_rule');record(b2,'B2_calendar_designation_rule',metric)
 b3<-fitset(period='strict',window=TRUE,tag='B3_nine_year_window');record(b3,'B3_nine_year_window',metric)
 b4<-fitset(period='strict',window=TRUE,postonly=metric=='Post_difference',tag='B4_target_observations');record(b4,'B4_target_observations',metric)
 iid<-b4;iid$vi<-iid$iid_vi
 record(iid,'U1_same_data_IID_model',metric,'model');record(iid,'U2_same_data_IID_CR2',metric)
 record(b4,'U3_same_data_HAC_model',metric,'model');record(b4,'U4_same_data_HAC_CR2',metric)
 # Source support/values checked on common years before full reconstructed support.
 sharedold<-sharednew<-fullnew<-list()
 for(id in ids){o<-rebuilt[[id]];n<-up[up$series_id==id,];years<-intersect(o$year,n$year)
 sharedold[[id]]<-common(o[o$year%in%years,],id,metric,'strict',TRUE,metric=='Post_difference','C1_common_source_years_old_values')
 sharednew[[id]]<-common(n[n$year%in%years,],id,metric,'strict',TRUE,metric=='Post_difference','C2_common_source_years_new_values')
 fullnew[[id]]<-common(n,id,metric,'strict',TRUE,metric=='Post_difference','C3_reconstructed_support')
 }
 for(pair in list(list(sharedold,'C1_common_source_years_old_values'),list(sharednew,'C2_common_source_years_new_values'),list(fullnew,'C3_reconstructed_support')))record(do.call(rbind,pair[[1]]),pair[[2]],metric)
}
# Membership factorial under identical proposed models and source data.
allnew<-list()
for(metric in c('Change','Post_difference')){
 keys<-unique(c(old$series_id,proposed$series_id[proposed$metric==metric]));ps<-list()
 for(id in keys){d<-up[up$series_id==id,];if(!nrow(d))next
  if(metric=='Change'&&!all(d$baseline_design_verified))next
  ps[[id]]<-common(d,id,metric,'strict',TRUE,metric=='Post_difference','membership_fixed')
 }
 x<-do.call(rbind,ps);oldids<-old$series_id;newids<-proposed$series_id[proposed$metric==metric]
 for(s in c('M0_old_roster','M1_removals_only','M2_additions_only','M3_final_roster')){
  keep<-switch(s,M0_old_roster=x$series_id%in%oldids,M1_removals_only=x$series_id%in%intersect(oldids,newids),M2_additions_only=x$series_id%in%union(oldids,newids),M3_final_roster=x$series_id%in%newids)
  record(x[keep,],s,metric)
 }
 allnew[[metric]]<-x
}
p<-do.call(rbind,res);p$q<-ave(p$p,interaction(p$stage,p$metric),FUN=function(x)p.adjust(x,'BH',n=max(9,sum(is.finite(x)))))
write.csv(p,file.path(out,'pooled.csv'),row.names=FALSE)
# Different effect tables have different optional columns.
cols<-unique(unlist(lapply(effectlog,names)));ee<-do.call(rbind,lapply(effectlog,function(x){for(n in setdiff(cols,names(x)))x[[n]]<-NA;x[cols]}));write.csv(ee,file.path(out,'effects.csv'),row.names=FALSE)
write.csv(do.call(rbind,obslog),file.path(out,'observations.csv'),row.names=FALSE)
if(length(warnlog))write.csv(do.call(rbind,warnlog),file.path(out,'warnings.csv'),row.names=FALSE)
expected<-read.csv('outputs/researcher_scope_20261101/pooled.csv');expected<-subset(expected,scenario=='Primary'&view=='All');jj<-merge(subset(p,stage=='M3_final_roster'),expected,by=c('metric','Taxa','Resp'));err<-sapply(c('estimate','low','high','q'),function(v)max(abs(jj[[paste0(v,'.x')]]-jj[[paste0(v,'.y')]]),na.rm=TRUE));stopifnot(nrow(jj)==nrow(expected),max(err)<1e-7)
write.csv(data.frame(check=names(err),max_error=err),file.path(out,'endpoint_verification.csv'),row.names=FALSE)
cat('Current table and proposed endpoint reproduced. Saved',nrow(p),'pooled comparisons.\n');print(err)

capture.output(warnings(),file=file.path(out,'all_warnings.txt'))
