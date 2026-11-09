# Reuse the same pooling specification; input effects are unchanged primary contrasts.
suppressPackageStartupMessages({library(metafor);library(clubSandwich)})
out<-'outputs/model_eligibility_20261101'
frozen<-jsonlite::fromJSON(file.path(out,'freeze_manifest.json'))
stopifnot(all(vapply(names(frozen),function(p) digest::digest(file=p,algo='sha256')==frozen[[p]],logical(1))))
e<-read.csv(file.path(out,'effects.csv'));scenarios<-unique(e$scenario)
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
