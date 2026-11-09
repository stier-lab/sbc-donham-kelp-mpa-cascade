# Primary PISCO fish rule selected by Adrian on 2026-10-31.
# Station membership/annual inputs are pinned by pisco_baseline_policy.json.
run_pisco_sheephead_analysis <- function(dat,resp_name) {
 root<-here::here();source(file.path(root,"code/R/primary_uncertainty_utils.R"),local=TRUE)
 defs<-uncertainty_definitions(root)
 source(file.path(root,"code/R/pisco_nonlinear_candidates.R"),local=TRUE)
 starts<-list()
  contrast <- function(fit,model,ends=c(0,11),year) {
    theta<-coef(fit)
    if(model=='Mean') {g<-setNames(1,names(theta));est<-unname(theta[1])}
    if(model=='Step') {g<-setNames(c(0,1),names(theta));est<-unname(theta[2])}
    if(model=='Linear') {g<-setNames(c(0,diff(ends)),names(theta));est<-sum(g*theta)}
    if(model %in% c('Sigmoid','Asymptotic')) {
      nd<-data.frame(time.model=ends,time_offset=ends+.01,time=ends)
      f<-function(par){e<-list2env(c(as.list(nd),as.list(setNames(par,names(theta)))),parent=fit$m$getEnv());diff(as.numeric(eval(fit$m$formula()[[3]],e)))}
      est<-f(theta);g<-drop(numDeriv::jacobian(f,theta));names(g)<-names(theta)
    }
    se<-sqrt(drop(t(g)%*%vcov(fit)%*%g))
    lag<-max(1L,floor(4*(length(year)/100)^(2/9)))
    vh<-uncertainty_hac(fit,year,lag)
    hs<-max(se,sqrt(drop(t(g)%*%vh%*%g)))
    stopifnot(is.finite(est),is.finite(se),se>0,is.finite(hs))
    data.frame(Mean=est,SE=se,HAC_SE=hs,Model=model)
  }
  fit_one<-function(d,model) {
    if(model %in% c("Sigmoid","Asymptotic")) {
      ans<-pisco_nonlinear_candidate(d,model);z<-ans$fit
      starts[[length(starts)+1]]<<-cbind(MPA=unique(d$CA_MPA_Name_Short),Resp=resp_name,ans$starts)
    } else z<-switch(model,Mean=lm(lnDiff~1,d),Step=lm(lnDiff~BA,d),Linear=lm(lnDiff~time.model,d))
    if(is.null(z))stop('No fitted model')
    if(model %in% c('Sigmoid','Asymptotic') && (!inherits(z,'nls')||!isTRUE(z$convInfo$isConv)))stop('Nonlinear fallback rejected')
    if(any(!is.finite(coef(z)))||any(!is.finite(vcov(z))))stop('Invalid coefficients/covariance')
    z
  }

 marker<-jsonlite::fromJSON(file.path(root,"data/harmonized/pisco_baseline_policy.json"))
 target<-marker$MPAs
 result<-run_ci_analysis(dat[!dat$CA_MPA_Name_Short %in% target,],"S. pulcher","PISCO",resp_name)
 candidates<-list()
 for(mpa in target) {
  d<-dat[dat$CA_MPA_Name_Short==mpa,];d<-d[order(d$year),]
  d$BA<-factor(d$BA,c("Before","After"));d$time.model<-pmax(d$year-2012,0)
  stopifnot(all((d$BA=="Before")== (d$year<2012)),sum(d$BA=="Before")>=3,sum(d$BA=="After")>=4)
  pre<-d[d$BA=="Before",];pt<-summary(lm(lnDiff~year,pre))$coefficients[2,4]
  if(!is.finite(pt)||pt<.05) {result<-rbind(result,run_ci_analysis(d,"S. pulcher","PISCO",resp_name));next}
  fits<-list();aic<-numeric()
  for(model in c("Step","Linear","Asymptotic","Sigmoid")) {
   f<-tryCatch(suppressWarnings(fit_one(d,model)),error=function(e)NULL)
   if(is.null(f))next
   cc<-tryCatch(contrast(f,model,defs$observed_contrast_times(d),d$year),error=function(e)NULL)
   if(is.null(cc))next
   k<-attr(logLik(f),"df");n<-nrow(d);v<-if(n>k+1) AIC(f)+2*k*(k+1)/(n-k-1) else Inf
   if(!is.finite(v))next
   fits[[model]]<-list(f=f,cc=cc);aic[model]<-v
   candidates[[length(candidates)+1]]<-data.frame(MPA=mpa,Resp=resp_name,Model=model,AICc=v,pretrend_p=pt)
  }
  stopifnot(length(aic)>0);model<-names(which.min(aic));selected<-fits[[model]]
  if(model=="Step") r<-add_step_effect_size(d,"S. pulcher",mpa,"PISCO",resp_name)
  else if(model=="Linear") r<-add_linear_effect_size(d,"S. pulcher",mpa,"PISCO",resp_name)
  else {
   cc<-selected$cc;n<-nrow(d);df<-df.residual(selected$f)
   r<-data.frame(Taxa="S. pulcher",MPA=mpa,Mean=cc$Mean,SE=cc$SE,SD=cc$SE*sqrt(n-1),
      CI=cc$SE*qt(.975,df),Model=model,Source="PISCO",Resp=resp_name,BA="Y",Primary="Y",
      Type="pBACIPS",LinearBefore="N",N=n,DW_stat=NA_real_,DW_pval=NA_real_)
  }
  stopifnot(nrow(r)==1,abs(r$Mean-selected$cc$Mean)<1e-8,abs(r$SE-selected$cc$SE)<1e-8)
  result<-rbind(result,r)
 }
 write.csv(dplyr::bind_rows(candidates),file.path(root,"tables",paste0("table_s_pisco_primary_candidates_",resp_name,".csv")),row.names=FALSE)
 write.csv(dplyr::bind_rows(starts),file.path(root,"tables",paste0("table_s_pisco_primary_starts_",resp_name,".csv")),row.names=FALSE)
 result
}
