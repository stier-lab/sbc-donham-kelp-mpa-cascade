suppressPackageStartupMessages({library(metafor);library(clubSandwich)})
source('code/R/primary_uncertainty_utils.R')
O<-'outputs/cross_site_seasons_20261101'
for(expr in as.list(parse('code/tools/refit_common_estimand.R')))if(is.call(expr)&&identical(expr[[1]],as.name('<-'))&&identical(expr[[2]],as.name('pool')))eval(expr)
x<-read.csv(file.path(O,'local_annual.csv'));coverage<-read.csv(file.path(O,'selected_coverage.csv'));coverage$full_history<-as.logical(coverage$full_history)
a<-read.csv('outputs/pisco_source_reconciliation_20261101/annual_input.csv');a<-a[a$scenario=='Supported_corrections',]
base<-read.csv('outputs/pisco_source_reconciliation_20261101/effects.csv');base<-base[base$scenario=='Supported_corrections',]
ids<-unique(a$series_id[a$Source=='PISCO' & a$MPA%in%coverage$MPA & ((a$Taxa=='M. pyrifera'&a$Resp=='Bio')|(a$Taxa%in%c('M. franciscanus','S. purpuratus','P. interruptus')&a$Resp=='Den'))]);stopifnot(length(ids)==20)
stopifnot(!any(base$series_id[base$metric=='Change']%in%ids));ba<-base[base$metric=='Change',];stopifnot(nrow(ba)==55);write.csv(ba,file.path(O,'unchanged_before_after.csv'),row.names=FALSE)
base<-base[base$metric=='Post_slope',];stopifnot(nrow(base)==169);base$variant<-'Previous_review';base$scenario<-'Previous_review'
fits<-list();used<-list();diag<-list();failed<-list();logs<-list()
fitone<-function(z,sc,variant){
 if(nrow(z)<3){failed[[length(failed)+1]]<<-data.frame(scenario=sc,variant,series_id=if(nrow(z))z$series_id[1] else '',reason='Fewer than three common years');return(NULL)}
 z<-z[order(z$year),];stopifnot(all(z$t>0),!anyDuplicated(z$year));f<-lm(lnDiff~t,z);v<-uncertainty_hac(f,z$year,3);se<-sqrt(v[2,2])
 if(!is.finite(se)||se<1e-12){failed[[length(failed)+1]]<<-data.frame(scenario=sc,variant,series_id=z$series_id[1],reason='Constant/exact fit: no sampling variance');return(NULL)}
 used[[length(used)+1]]<<-data.frame(scenario=sc,variant,series_id=z$series_id,year=z$year,t=z$t,lnDiff=z$lnDiff)
 diag[[length(diag)+1]]<<-data.frame(scenario=sc,variant,series_id=z$series_id[1],max_cook=max(cooks.distance(f)),n=nrow(z),df=df.residual(f))
 q<-qt(.975,df.residual(f));data.frame(series_id=z$series_id[1],MPA=z$MPA[1],Taxa=z$Taxa[1],Resp=z$Resp[1],Source=z$Source[1],cluster=z$cluster[1],scenario=sc,metric='Post_slope',n=nrow(z),nb=0,na=nrow(z),estimate=unname(coef(f)[2]),se,vi=se^2,df=df.residual(f),first_after=min(z$year),last_after=max(z$year),first_age=min(z$t),last_age=max(z$t),before_years='',after_years=paste(z$year,collapse=';'),variant,low=unname(coef(f)[2])-q*se,high=unname(coef(f)[2])+q*se)
}
for(sc in unique(coverage$scenario))for(id in ids){
 z<-x[x$scenario==sc & x$series_id==id,];if(!nrow(z))next
 old<-a[a$series_id==id&a$t>0,];ys<-intersect(old$year,z$year)
 for(v in c('Recovered','Old_common_years','Corrected_common_years')){
  zz<-if(v=='Recovered')z else if(v=='Old_common_years')old[old$year%in%ys,] else z[z$year%in%ys,];r<-fitone(zz,sc,v);if(!is.null(r))fits[[length(fits)+1]]<-r
 }
}
f<-do.call(rbind,fits);write.csv(f,file.path(O,'local_effects.csv'),row.names=FALSE);write.csv(do.call(rbind,used),file.path(O,'model_rows.csv'),row.names=FALSE);write.csv(do.call(rbind,diag),file.path(O,'linear_diagnostics.csv'),row.names=FALSE);write.csv(do.call(rbind,failed),file.path(O,'failed_local_fits.csv'),row.names=FALSE)
cols<-intersect(names(base),names(f));base<-base[,cols];e<-list(base)
for(sc in unique(coverage$scenario)){
 pp<-f[f$scenario==sc & f$variant=='Recovered',];other<-base[!base$series_id%in%ids,];fullmpa<-coverage$MPA[coverage$scenario==sc & coverage$full_history]
 for(v in c('Recovered','Retention_only','Old_common_years','Corrected_common_years','Full_history_only')){
  add<-if(v=='Retention_only')base[base$series_id%in%pp$series_id,] else if(v=='Full_history_only')pp[pp$MPA%in%fullmpa,cols] else f[f$scenario==sc&f$variant==v,cols]
  # Match the two common-year steps on successful estimable fits in both versions.
  if(v%in%c('Old_common_years','Corrected_common_years')){
   common<-intersect(f$series_id[f$scenario==sc&f$variant=='Old_common_years'],f$series_id[f$scenario==sc&f$variant=='Corrected_common_years']);add<-add[add$series_id%in%common,]
  }
  z<-rbind(other,add);z$scenario<-sc;z$variant<-v;e[[length(e)+1]]<-z
 }
}
e<-do.call(rbind,e);e$rule<-e$scenario;e$scenario<-paste(e$rule,e$variant,sep='__');write.csv(e,file.path(O,'effects.csv'),row.names=FALSE)
# Reuse numerically identical inputs across rules; never reuse across differing membership/estimates/variance.
cache<-new.env();pools<-list()
for(sc in unique(e$scenario))for(tx in unique(e$Taxa))for(resp in c('Bio','Den')){
 z<-e[e$scenario==sc&e$Taxa==tx&e$Resp==resp,];if(!nrow(z))next;z<-z[order(z$series_id),]
 key<-digest::digest(z[,c('series_id','MPA','cluster','estimate','vi')]);
 if(exists(key,cache,inherits=FALSE)){r<-get(key,cache);r$scenario<-sc}else{r<-pool(z,'All');assign(key,r,cache)}
 r$rule<-z$rule[1];r$variant<-z$variant[1];pools[[length(pools)+1]]<-r
}
p<-do.call(rbind,pools);p$q<-ave(p$p,p$scenario,FUN=function(x)p.adjust(x,'BH',n=9));write.csv(p,file.path(O,'pooled.csv'),row.names=FALSE);write.csv(do.call(rbind,logs),file.path(O,'fit_log.csv'),row.names=FALSE)
oldp<-read.csv('outputs/pisco_source_reconciliation_20261101/pooled.csv');j<-merge(p[p$variant=='Previous_review',],oldp[oldp$scenario=='Supported_corrections'&oldp$metric=='Post_slope',],by=c('Taxa','Resp'));stopifnot(nrow(j)==9,max(abs(j$estimate.x-j$estimate.y))<1e-7,max(abs(j$q.x-j$q.y),na.rm=TRUE)<1e-7)
capture.output(sessionInfo(),file=file.path(O,'session_info.txt'));cat(nrow(f),'local fits;',nrow(p),'pools;',length(ls(cache)),'distinct pools; previous review reproduced\n')
