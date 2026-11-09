# Explicit-period initialization, unchanged Hill/asymptotic forms, all-start selection.
# Bounds follow the existing PISCO Hill fit. Bound-touching fits are diagnostic only.
pisco_nonlinear_candidate <- function(d,model) {
 d$time_offset<-d$time.model+.01
 base<-mean(d$lnDiff[d$BA=='Before']);change<-mean(d$lnDiff[d$BA=='After'])-base
 max_time<-max(d$time.model);fits<-list();records<-list()
 grid<-expand.grid(L=unique(pmin(c(.5,1,3,6,11,20),max_time*1.8)),
                   K=if(model=='Sigmoid')c(.5,1,2,5,10) else 1,
                   mult=c(.5,1,2))
 for(i in seq_len(nrow(grid))) {
  st<-c(M=change*grid$mult[i],B=base,L=grid$L[i])
  lo<-c(M=-50,B=-50,L=.2);hi<-c(M=50,B=50,L=max_time*2)
  if(model=='Sigmoid') {st<-c(st,K=grid$K[i]);lo<-c(lo,K=.2);hi<-c(hi,K=15)}
  formula<-if(model=='Sigmoid')lnDiff~M*(time_offset/L)^K/(1+(time_offset/L)^K)+B else lnDiff~M*time.model/(L+time.model)+B
  warning_text<-character()
  f<-tryCatch(withCallingHandlers(minpack.lm::nlsLM(formula,data=d,start=as.list(st),lower=lo,upper=hi,
       control=minpack.lm::nls.lm.control(maxiter=500)),warning=function(w){warning_text<<-c(warning_text,conditionMessage(w));invokeRestart('muffleWarning')}),error=identity)
  converged<-inherits(f,'nls')&&isTRUE(f$convInfo$isConv)
  bound<-FALSE;cov_ok<-FALSE;rss<-NA_real_
  if(converged) {
   cf<-coef(f);bound<-any((cf-lo)/(hi-lo)<1e-5 | (hi-cf)/(hi-lo)<1e-5)
   cov_ok<-tryCatch(all(is.finite(vcov(f))),error=function(e)FALSE);rss<-sum(resid(f)^2)
  }
  eligible<-converged&&!bound&&cov_ok
  records[[i]]<-data.frame(model,start=i,L_start=st['L'],K_start=if('K'%in%names(st))st['K'] else NA,
    converged,bound,cov_ok,eligible,rss,warnings=paste(unique(warning_text),collapse=';'),error=if(inherits(f,'error'))conditionMessage(f) else '')
  if(eligible)fits[[as.character(i)]]<-f
 }
 log<-do.call(rbind,records)
 if(!length(fits))return(list(fit=NULL,starts=log))
 best<-names(which.min(vapply(fits,function(f)sum(resid(f)^2),1.0)))
 list(fit=fits[[best]],starts=log)
}
