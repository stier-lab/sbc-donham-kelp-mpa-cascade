suppressPackageStartupMessages({library(metafor);library(clubSandwich)})
source('code/R/primary_uncertainty_utils.R')
out<-'outputs/seasonal_zone_audit_20261101';d<-read.csv(file.path(out,'annual_input.csv'));d$baseline_design_verified<-as.logical(d$baseline_design_verified);stopifnot(!anyNA(d$baseline_design_verified))
for(expr in as.list(parse('code/tools/refit_common_estimand.R')))if(is.call(expr)&&identical(expr[[1]],as.name('<-'))&&identical(expr[[2]],as.name('pool')))eval(expr)
rows<-list();used<-list();cov<-list();logs<-list();diagnostics<-list()
for(id in unique(d$series_id))for(scenario in c('Original','Within_station_zones','Fixed_station_zones','Fixed_all_program_units')){
 z<-d[d$series_id==id&d$scenario==scenario,];z<-z[order(z$year),];z<-z[z$t!=0,];if(scenario=='Nine_year')z<-z[z$t<0|z$t<=9,]
 pre<-z$t<0;post<-z$t>0;nb<-sum(pre);na<-sum(post);ba<-nb>=3&&all(z$baseline_design_verified);lag<-switch(scenario,Full_HAC1=1,Full_HAC5=5,3)
 for(metric in c('Change','Post_slope','BA_post_slope','No_baseline_post_slope','Slope_change','Trend_adjusted_change')){
  needba<-metric%in%c('Change','BA_post_slope','Slope_change','Trend_adjusted_change')
  eligible<-na>=3&&(!needba||ba)&&(metric!='No_baseline_post_slope'||!ba)
  reason<-if(na<3)'Fewer than three post observations' else if(needba&&!ba)'No qualifying verified baseline' else if(metric=='No_baseline_post_slope'&&ba)'Assigned to baseline-supported group' else 'Eligible'
  cov[[length(cov)+1]]<-data.frame(series_id=id,scenario,metric,nb,na,eligible,reason)
  if(!eligible)next
  x<-if(metric%in%c('Post_slope','BA_post_slope','No_baseline_post_slope'))z[post,] else z
  x$post<-as.integer(x$t>0);x$pt<-x$post*x$t
  form<-if(metric%in%c('Post_slope','BA_post_slope','No_baseline_post_slope'))lnDiff~t else if(metric=='Change')lnDiff~post+pt else lnDiff~t+post+pt
  fit<-lm(form,x);if(fit$rank!=length(coef(fit))||df.residual(fit)<1){cov[[length(cov)]]$eligible<-FALSE;cov[[length(cov)]]$reason<-'Rank/df insufficient';next}
  g<-if(metric%in%c('Post_slope','BA_post_slope','No_baseline_post_slope'))c(0,1) else if(metric=='Change')c(0,1,mean(x$t[x$post==1])) else if(metric=='Slope_change')c(0,0,0,1) else c(0,0,1,mean(x$t[x$post==1]))
  v<-uncertainty_hac(fit,x$year,lag);est<-sum(g*coef(fit));se<-sqrt(drop(g%*%v%*%g));
  if(!is.finite(se)||se<1e-12){cov[[length(cov)]]$eligible<-FALSE;cov[[length(cov)]]$reason<-'Constant/exact-fit response: no estimable sampling variance';next}
  if(metric=='Change')stopifnot(abs(est-(mean(x$lnDiff[x$post==1])-mean(x$lnDiff[x$post==0])))<1e-8)
  rows[[length(rows)+1]]<-data.frame(series_id=id,MPA=x$MPA[1],Taxa=x$Taxa[1],Resp=x$Resp[1],Source=x$Source[1],cluster=x$cluster[1],scenario,metric,n=nrow(x),nb=if(needba)nb else 0,na,estimate=est,se,vi=se^2,df=df.residual(fit),first_after=min(x$year[x$post==1]),last_after=max(x$year[x$post==1]),first_age=min(x$t[x$post==1]),last_age=max(x$t[x$post==1]),before_years=paste(x$year[x$post==0],collapse=';'),after_years=paste(x$year[x$post==1],collapse=';'))
  used[[length(used)+1]]<-data.frame(series_id=id,scenario,metric,year=x$year,lnDiff=x$lnDiff,t=x$t,post=x$post,lag)
  diagnostics[[length(diagnostics)+1]]<-data.frame(series_id=id,scenario,metric,max_cook=max(cooks.distance(fit)),n=nrow(x),df=df.residual(fit))
 }
}
e<-do.call(rbind,rows);u<-do.call(rbind,used)
# Hold cohort membership fixed using the original full-followup eligibility.
base<-read.csv('outputs/full_followup_recovery_20261101/effects.csv');base<-base[base$scenario=='Full',]
e<-merge(e,unique(base[,c('series_id','metric')]),by=c('series_id','metric'))
early<-base[base$metric=='Post_slope' & base$na>=5 & base$last_after-base$first_after>=5 & base$first_age<=3,'series_id']
extra<-e[e$metric=='Post_slope' & e$series_id %in% early,];extra$scenario<-paste0(extra$scenario,'_early140');e<-rbind(e,extra)
pools<-list()
for(s in unique(e$scenario))for(mt in unique(e$metric))for(tx in unique(e$Taxa))for(resp in c('Bio','Den')){
 x<-e[e$scenario==s&e$metric==mt&e$Taxa==tx&e$Resp==resp,];if(nrow(x))pools[[length(pools)+1]]<-pool(x,'All')
}
p<-do.call(rbind,pools);capture.output(str(p),file=file.path(out,'pool_structure.txt'));p$q<-ave(p$p,interaction(p$scenario,p$metric),FUN=function(x)p.adjust(x,'BH',n=max(9,sum(is.finite(x)))))
write.csv(e,file.path(out,'effects.csv'),row.names=FALSE);write.csv(u,file.path(out,'model_rows.csv'),row.names=FALSE);write.csv(do.call(rbind,cov),file.path(out,'coverage.csv'),row.names=FALSE);write.csv(p,file.path(out,'pooled.csv'),row.names=FALSE);write.csv(do.call(rbind,logs),file.path(out,'fit_log.csv'),row.names=FALSE);write.csv(do.call(rbind,diagnostics),file.path(out,'linear_diagnostics.csv'),row.names=FALSE)
old<-read.csv('outputs/full_followup_recovery_20261101/pooled.csv');j<-merge(p[p$scenario=='Original',],old[old$scenario=='Full',],by=c('metric','Taxa','Resp'));stopifnot(nrow(j)==nrow(old[old$scenario=='Full',]),max(abs(j$estimate.x-j$estimate.y))<1e-7,max(abs(j$q.x-j$q.y),na.rm=TRUE)<1e-7)
capture.output(sessionInfo(),file=file.path(out,'session_info.txt'));cat('Original full-followup reproduced;',nrow(e),'effects;',nrow(p),'pools.\n')
