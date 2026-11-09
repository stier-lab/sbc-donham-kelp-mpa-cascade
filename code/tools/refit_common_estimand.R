#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(metafor);library(clubSandwich);library(ggplot2)})
source('code/R/primary_uncertainty_utils.R')
args<-commandArgs(trailingOnly=TRUE)
out<-if(length(args))args[1] else 'outputs/pooled_common_estimand_20261101'
frozen<-jsonlite::fromJSON(file.path(out,'freeze_manifest.json'))
stopifnot(all(vapply(names(frozen),function(p) digest::digest(file=p,algo='sha256')==frozen[[p]],logical(1))))
d<-read.csv(file.path(out,'annual_input.csv'));d$baseline_design_verified<-as.logical(d$baseline_design_verified);stopifnot(!anyNA(d$baseline_design_verified)); rows<-list();coverage<-list();used<-list()
scenarios<-c('Primary','Full_after','Baseline10','Original_timing','HAC1','HAC5','Trend')
if(grepl('diver_exclusion_refit|upstream_sites',out))scenarios<-c(scenarios,'Omit_Harris_lobster','Omit_Gull_KFM')
if(grepl('upstream_sites',out))scenarios<-c(scenarios,'Omit_Anacapa_KFM')
for(id in unique(d$series_id)) for(scenario in scenarios) {
 z<-d[d$series_id==id,];z<-z[order(z$year),]
 if(scenario=='Omit_Anacapa_KFM' && z$MPA[1]=='Anacapa Island SMCA' && z$Source[1]=='KFM')next
 if(scenario=='Omit_Harris_lobster' && id=='Harris Point SMR|P. interruptus|KFM|Den')next
 if(scenario=='Omit_Gull_KFM' && z$MPA[1]=='Gull Island SMR' && z$Source[1]=='KFM')next
 pre<-z$t<0;post<-z$t>=1 & (if(scenario=='Full_after') TRUE else z$t<=9)
 if(scenario=='Baseline10')pre<-pre & z$t>= -10
 if(scenario=='Original_timing'){pre<-z$BA=='Before';post<-z$BA=='After' & z$t<=9}
 nb<-sum(pre);na<-sum(post);verified<-all(z$baseline_design_verified)
 for(metric in c('Change','Post_difference')) {
  eligible<-na>=3 && (metric=='Post_difference' || (nb>=3 && verified))
  reason<-if(na<3)'Fewer than3post years' else if(metric=='Change' && nb<3)'Fewer than3before years' else if(metric=='Change'&&!verified)'PISCO fixed baseline design not established' else 'Eligible'
  coverage[[length(coverage)+1]]<-data.frame(series_id=id,MPA=z$MPA[1],Taxa=z$Taxa[1],Resp=z$Resp[1],Source=z$Source[1],scenario,metric,nb,na,eligible,reason)
  if(!eligible)next
  x<-z[if(metric=='Change')pre|post else post,];x$post<-as.integer(x$year %in% z$year[post]);x$pt<-x$post*x$t
  fit<-if(metric=='Post_difference')lm(lnDiff~t,x) else if(scenario=='Trend')lm(lnDiff~t+post+pt,x) else lm(lnDiff~post+pt,x)
  stopifnot(fit$rank==length(coef(fit)))
  X<-model.matrix(fit)
  if(metric=='Post_difference') g<-colMeans(X) else {x0<-x;x0$post<-0;x0$pt<-0;X0<-model.matrix(delete.response(terms(fit)),x0);g<-colMeans((X-X0)[x$post==1,,drop=FALSE])}
  lag<-if(scenario=='HAC1')1 else if(scenario=='HAC5')5 else 3
  V<-uncertainty_hac(fit,x$year,lag);estimate<-sum(g*coef(fit));se<-sqrt(drop(g%*%V%*%g));if(!is.finite(se)||se<=1e-12){
   coverage[[length(coverage)]]$eligible<-FALSE;coverage[[length(coverage)]]$reason<-'No estimable sampling variance (constant or exact-fit series)';next
  }
  direct<-mean(x$lnDiff[x$post==1])-if(metric=='Change')mean(x$lnDiff[x$post==0]) else 0
  if(scenario!='Trend'||metric=='Post_difference')stopifnot(abs(estimate-direct)<1e-8)
  rows[[length(rows)+1]]<-data.frame(series_id=id,MPA=x$MPA[1],Taxa=x$Taxa[1],Resp=x$Resp[1],Source=x$Source[1],cluster=x$cluster[1],scenario,metric,n=nrow(x),nb=if(metric=='Change')nb else 0,na,estimate,se,vi=se^2,df=df.residual(fit),before_years=paste(x$year[x$post==0],collapse=';'),after_years=paste(x$year[x$post==1],collapse=';'))
  used[[length(used)+1]]<-data.frame(series_id=id,scenario,metric,year=x$year,lnDiff=x$lnDiff,t=x$t,post=x$post,pt=x$pt,lag)
 }
}
e<-do.call(rbind,rows);write.csv(e,file.path(out,'effects.csv'),row.names=FALSE)
write.csv(do.call(rbind,coverage),file.path(out,'coverage.csv'),row.names=FALSE);write.csv(do.call(rbind,used),file.path(out,'model_rows.csv'),row.names=FALSE)
pooled<-list();logs<-list()
pool<-function(x,view) {
 k<-nrow(x);ng<-length(unique(x$cluster));nm<-length(unique(x$MPA))
 base<-data.frame(scenario=x$scenario[1],metric=x$metric[1],Taxa=x$Taxa[1],Resp=x$Resp[1],view,k,n_mpa=nm,n_clusters=ng)
 if(k<2 || nm<2) return(cbind(base,estimate=if(k==1)x$estimate else NA,se=NA,df=NA,low=NA,high=NA,p=NA,model_low=NA,model_high=NA,sigma_mpa=NA,sigma_series=NA,status='No pooled inference: one MPA',boundary=NA))
 random<-if(anyDuplicated(x$MPA))list(~1|MPA,~1|series_id) else list(~1|MPA)
 warns<-character()
 m<-withCallingHandlers(rma.mv(estimate,V=vi,random=random,data=x,method='REML',test='t'),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')})
 stopifnot(all(is.finite(coef(m))))
 ct<-tryCatch(withCallingHandlers(coef_test(m,vcov='CR2',cluster=x$cluster,test='Satterthwaite'),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')}),error=function(err)err)
 status<-'CR2';se<-df<-lo<-hi<-p<-NA_real_
 if(inherits(ct,'error'))status<-paste('CR2 unavailable:',conditionMessage(ct)) else {
  se<-ct$SE[1];df<-ct$df_Satt[1]
  if(ng<3||!is.finite(df)||df<2)status<-'Insufficient independent information: fewer than3clusters or df<2' else {lo<-coef(m)[1]-qt(.975,df)*se;hi<-coef(m)[1]+qt(.975,df)*se;p<-ct$p_Satt[1];if(df<4)status<-'CR2: few effective degrees of freedom (<4)'}
 }
 logs[[length(logs)+1]]<<-cbind(base,warnings=paste(warns,collapse=';'),sigma=paste(m$sigma2,collapse=';'))
 cbind(base,estimate=unname(coef(m)[1]),se,df,low=lo,high=hi,p,model_low=m$ci.lb[1],model_high=m$ci.ub[1],sigma_mpa=m$sigma2[1],sigma_series=if(length(m$sigma2)>1)m$sigma2[2] else NA,status,boundary=any(m$sigma2<1e-8))
}
for(scenario in scenarios) for(metric in c('Change','Post_difference')) for(resp in c('Bio','Den')) for(taxon in unique(e$Taxa)) {
 x<-e[e$scenario==scenario & e$metric==metric & e$Resp==resp & e$Taxa==taxon,]
 stopifnot(!anyDuplicated(x$series_id))
 if(!nrow(x))next
 views<-list(All=x,Divers=x[x$Source!='Landsat',],Satellite=x[x$Source=='Landsat',])
 for(view in names(views))if(nrow(views[[view]]))pooled[[length(pooled)+1]]<-pool(views[[view]],view)
}
p<-do.call(rbind,pooled)
stopifnot(!anyDuplicated(p[c('scenario','metric','Taxa','Resp','view')]))
p$q<-ave(p$p,interaction(p$scenario,p$metric,p$view),FUN=function(x)p.adjust(x,'BH',n=max(9,sum(is.finite(x)))))
write.csv(p,file.path(out,'pooled.csv'),row.names=FALSE);write.csv(do.call(rbind,logs),file.path(out,'fit_log.csv'),row.names=FALSE)
# Common years between satellite and each available baseline-eligible diver kelp history.
matched<-list()
ke<-subset(e,scenario=='Primary'&metric=='Change'&Taxa=='M. pyrifera'&Resp=='Bio')
for(mpa in intersect(ke$MPA[ke$Source=='Landsat'],ke$MPA[ke$Source!='Landsat'])) {
 ids<-ke$series_id[ke$MPA==mpa];x<-d[d$series_id%in%ids,];years<-Reduce(intersect,split(x$year,x$series_id));x<-x[x$year%in%years & (x$t<0|x$t%in%1:9),]
 for(id in unique(x$series_id)) {z<-x[x$series_id==id,];nb<-sum(z$t<0);na<-sum(z$t>0);if(nb>=3&&na>=3)matched[[length(matched)+1]]<-data.frame(MPA=mpa,Source=z$Source[1],nb,na,before_years=paste(z$year[z$t<0],collapse=';'),after_years=paste(z$year[z$t>0],collapse=';'),change=mean(z$lnDiff[z$t>0])-mean(z$lnDiff[z$t<0]))}
}
write.csv(do.call(rbind,matched),file.path(out,'matched_kelp_years.csv'),row.names=FALSE)
labels<-c('M. pyrifera'='Kelp','M. franciscanus'='Red urchin','S. purpuratus'='Purple urchin','P. interruptus'='Lobster','S. pulcher'='Sheephead')
label<-function(x)paste(labels[x$Taxa],ifelse(x$Resp=='Bio','biomass','density'))
x<-subset(e,scenario=='Primary');x$label<-label(x)
y<-subset(p,scenario=='Primary'&view=='All');y$label<-label(y)
levels<-rev(c('Kelp biomass','Sheephead biomass','Lobster biomass','Red urchin biomass','Purple urchin biomass','Sheephead density','Lobster density','Red urchin density','Purple urchin density'))
x$label<-factor(x$label,levels);y$label<-factor(y$label,levels)
y$precision <- ifelse(is.finite(y$low),'Interval available','Descriptive only')
x$metric<-factor(x$metric,c('Change','Post_difference'),labels=c('Before–After change','After-protection difference'));y$metric<-factor(y$metric,c('Change','Post_difference'),labels=levels(x$metric))
th<-theme_minimal(base_size=12)+theme(legend.position='none',panel.grid.minor=element_blank(),plot.margin=margin(20,20,20,20),panel.grid.major.y=element_blank())
g<-ggplot(x,aes(estimate,label))+geom_vline(xintercept=0,color='grey65',linetype='dashed')+geom_point(aes(color=Source),position=position_jitter(width=0,height=.14,seed=68),alpha=.55,size=1.8)+geom_errorbar(data=y,aes(xmin=low,xmax=high),orientation='y',width=.22,linewidth=.65,na.rm=TRUE)+geom_point(data=y[y$n_mpa>1,],aes(shape=precision),size=3)+scale_shape_manual(values=c('Interval available'=18,'Descriptive only'=5))+facet_wrap(~metric,nrow=1)+scale_color_manual(values=c(Landsat='#D55E00',KFM='#0072B2',LTER='#009E73',PISCO='#CC79A7'))+labs(x='Log protected/reference ratio or change in that ratio',y=NULL)+th
ggsave(file.path(out,'common_estimands.png'),g,width=13,height=7.5,dpi=180);ggsave(file.path(out,'common_estimands.pdf'),g,width=13,height=7.5)
cat('Saved',nrow(e),'effects;',nrow(p),'pooled rows\n');print(subset(p,scenario=='Primary'&view=='All')[c('Taxa','Resp','metric','k','n_mpa','estimate','low','high','df','q','status')],row.names=FALSE)

capture.output(if(is.null(warnings())) 'No uncaptured warnings' else warnings(),file=file.path(out,'warnings.txt'))
