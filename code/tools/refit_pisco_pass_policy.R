suppressPackageStartupMessages({library(metafor);library(clubSandwich)})
source('code/R/primary_uncertainty_utils.R')
O<-'outputs/pisco_pass_crane_20261101';d<-read.csv(file.path(O,'model_input.csv'));rows<-list();logs<-list();used<-list();coverage<-list()
for(expr in as.list(parse('code/tools/refit_common_estimand.R')))if(is.call(expr)&&identical(expr[[1]],as.name('<-'))&&identical(expr[[2]],as.name('pool')))eval(expr)
for(sc in unique(d$scenario))for(id in unique(d$series_id)){
 z<-d[d$scenario==sc & d$series_id==id & d$t!=0,];z<-z[order(z$year),];if(!nrow(z))next
 nb<-sum(z$t<0);na<-sum(z$t>0);ba<-nb>=3&&all(z$baseline_design_verified)
 for(metric in c('Change','Post_slope')){
  ok<-na>=3 && (metric=='Post_slope'||ba)
  coverage[[length(coverage)+1]]<-data.frame(series_id=id,scenario=sc,metric,nb,na,eligible=ok)
  if(!ok)next
  x<-if(metric=='Post_slope')z[z$t>0,] else z;x$post<-as.integer(x$t>0);x$pt<-x$post*x$t
  fit<-if(metric=='Post_slope')lm(lnDiff~t,x) else lm(lnDiff~post+pt,x)
  if(fit$rank!=length(coef(fit))||df.residual(fit)<1)next
  g<-if(metric=='Post_slope')c(0,1) else c(0,1,mean(x$t[x$post==1]));V<-uncertainty_hac(fit,x$year,3);est<-sum(g*coef(fit));se<-sqrt(drop(g%*%V%*%g))
  if(!is.finite(se)||se<1e-12)next
  rows[[length(rows)+1]]<-data.frame(series_id=id,MPA=x$MPA[1],Source=x$Source[1],Taxa=x$Taxa[1],Resp=x$Resp[1],cluster=x$cluster[1],scenario=sc,contrast=x$contrast[1],arm=x$arm[1],metric,n=nrow(x),nb=if(metric=='Change')nb else 0,na,estimate=est,se,vi=se^2,df=df.residual(fit),low=est-qt(.975,df.residual(fit))*se,high=est+qt(.975,df.residual(fit))*se)
  used[[length(used)+1]]<-data.frame(series_id=id,scenario=sc,metric,year=x$year,t=x$t,lnDiff=x$lnDiff)
 }
}
e<-do.call(rbind,rows);pools<-list();pool_checks<-list()
pool_bfgs<-pool;body(pool_bfgs)<-parse(text=gsub('test = "t"','test = "t", control=list(optimizer="optim",optmethod="BFGS")',paste(deparse(body(pool)),collapse='\n'),fixed=TRUE))[[1]]
for(sc in unique(e$scenario))for(mt in unique(e$metric))for(resp in c('Bio','Den'))for(view in c('All_programs','PISCO')){
 x<-e[e$scenario==sc & e$metric==mt & e$Resp==resp,];if(view=='PISCO')x<-x[x$Source=='PISCO',]
 if(nrow(x)){
  ans<-tryCatch(pool(x,view),error=function(e)e);solver<-'nlminb';first_error<-''
  if(inherits(ans,'error')){first_error<-conditionMessage(ans);solver<-'BFGS fallback';ans<-tryCatch(pool_bfgs(x,view),error=function(e)e)}
  pool_checks[[length(pool_checks)+1]]<-data.frame(scenario=sc,metric=mt,Resp=resp,view,solver,first_error,status=if(inherits(ans,'error'))conditionMessage(ans) else 'fitted')
  if(!inherits(ans,'error'))pools[[length(pools)+1]]<-ans
 }
}
p<-do.call(rbind,pools);write.csv(do.call(rbind,pool_checks),file.path(O,'pool_checks.csv'),row.names=FALSE)
write.csv(e,file.path(O,'effects.csv'),row.names=FALSE);write.csv(p,file.path(O,'pooled.csv'),row.names=FALSE);write.csv(do.call(rbind,used),file.path(O,'model_rows.csv'),row.names=FALSE);write.csv(do.call(rbind,coverage),file.path(O,'model_coverage.csv'),row.names=FALSE);write.csv(do.call(rbind,logs),file.path(O,'fit_log.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(O,'session_info.txt'));cat(nrow(e),'effects;',nrow(p),'pools\n')
