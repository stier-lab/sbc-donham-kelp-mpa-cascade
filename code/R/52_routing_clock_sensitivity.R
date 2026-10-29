# Fixed-cohort checks of PISCO eligibility and model clocks. Does not change primary files.
# Run: Rscript code/R/52_routing_clock_sensitivity.R
local({
  suppressPackageStartupMessages({library(here);library(dplyr);library(metafor)})
  root <- here::here()
  source(file.path(root,'code/R/primary_uncertainty_utils.R'),local=TRUE)
  files <- c('data/sumstats_final.csv','data/harmonized/harmonized_response_ratios.csv',
    'data/harmonized/harmonized_landsat_rr.csv','tables/table_02_meta_analysis.csv',
    'tables/table_03_cross_taxa_meta_regression.csv','code/R/08_effect_sizes.R',
    'code/R/02_pBACIPS_function.R','code/R/primary_uncertainty_utils.R',
    'code/R/52_routing_clock_sensitivity.R')
  hashes <- tools::md5sum(file.path(root,files))
  write_out <- function(x,name) write.csv(x,file.path(root,'tables',paste0('table_s_routing_clock_',name,'.csv')),row.names=FALSE)
  primary <- read.csv(file.path(root,files[1])); primary$key <- uncertainty_key(primary)
  stopifnot(nrow(primary)==145L,!anyDuplicated(primary$key))
  annual <- uncertainty_annual(root,primary); defs <- uncertainty_definitions(root)
  warnlog <- list(); safe <- function(expr,label) {
    w <- character();result <- tryCatch(withCallingHandlers(expr,warning=function(x){w<<-c(w,conditionMessage(x));invokeRestart('muffleWarning')}),error=identity)
    warnlog[[length(warnlog)+1L]] <<- data.frame(label=label,warnings=paste(unique(w),collapse=';'),error=if(inherits(result,'error'))conditionMessage(result) else '')
    if(inherits(result,'error')) return(NULL)
    result
  }
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
    z<-switch(model,Mean=lm(lnDiff~1,d),Step=lm(lnDiff~BA,d),Linear=lm(lnDiff~time.model,d),
      Sigmoid=defs$mySIGfun_standalone(d$lnDiff,d$time.model,max(which(d$time.model==0)),seq_len(nrow(d))),
      Asymptotic=defs$myASYfun_standalone(d$lnDiff,d$time.model,max(which(d$time.model==0)),seq_len(nrow(d))))
    if(is.null(z))stop('No fitted model')
    if(model %in% c('Sigmoid','Asymptotic') && (!inherits(z,'nls')||!isTRUE(z$convInfo$isConv)))stop('Nonlinear fallback rejected')
    if(any(!is.finite(coef(z)))||any(!is.finite(vcov(z))))stop('Invalid coefficients/covariance')
    z
  }
  rebuilt<-setNames(vector('list',nrow(primary)),primary$key)
  coverage<-list();clock_audit<-list();baseline_hac<-numeric(nrow(primary))
  for(i in seq_len(nrow(primary))) {
    row<-primary[i,];r<-safe(uncertainty_reconstruct(row,annual,defs),row$key)
    if(is.null(r))stop('Primary reconstruction failed: ',row$key)
    rebuilt[[row$key]]<-r;d<-r$data;start<-row$MPA_Start
    baseline_hac[i]<-contrast(r$fit,row$Model,year=d$year)$HAC_SE
    pre<-d[d$BA=='Before',]; strict<-d[d$year<start,]
    pretest<-function(x) {if(nrow(x)<3 || sd(x$lnDiff)==0)return(NA_real_);summary(lm(lnDiff~year,x))$coefficients[2,4]}
    pre_index_p<-if(nrow(pre)>=3 && sd(pre$lnDiff)>0)summary(lm(lnDiff~seq_len(nrow(pre)),pre))$coefficients[2,4] else NA_real_
    coverage[[i]]<-data.frame(row[c('key','MPA','Taxa','Source','Resp','Model','AnalysisType','MPA_Start')],
      n_years=nrow(d),n_before=nrow(pre),n_strict_before=nrow(strict),n_after=sum(d$BA=='After'),
      before_years=paste(pre$year,collapse=';'),strict_before_years=paste(strict$year,collapse=';'),
      designation_year_present=start %in% d$year,pretrend_p=pretest(pre),pretrend_index_p=pre_index_p,strict_pretrend_p=pretest(strict),
      constant_before=nrow(pre)>1 && sd(pre$lnDiff)==0,
      observed_years=paste(d$year,collapse=';'),input_time_matches_calendar=all(abs(d$input_time-pmax(d$year-start,0))<1e-8))
    trend<-row$Model %in% c('Linear','Sigmoid','Asymptotic')
    post<-d[d$BA=='After',];first<-min(post$year);last<-max(post$year[post$year<=first+11])
    clock_audit[[i]]<-data.frame(row[c('key','MPA','Taxa','Source','Resp','Model','AnalysisType')],
      clock=if(row$Model=='Step')'period contrast; clock not used' else r$clock,first_after_year=first,last_after_year=max(post$year),n_before=nrow(pre),
      native_max=max(d$time.model),native_start_0_outside=trend && min(d$time.model)>0,
      native_end_11_outside=trend && max(d$time.model)<11,native_contrast_0_11_outside=trend && (min(d$time.model)>0 || max(d$time.model)<11),
      calendar_contrast_0_11_outside=trend && (min(pmax(d$year-start,0))>0 || max(pmax(d$year-start,0))<11),
      matched_start_year=first,matched_end_year=last,matched_calendar_span=last-first,
      native_start=post$time.model[post$year==first],native_end=post$time.model[post$year==last])
  }
  coverage<-bind_rows(coverage);clock_audit<-bind_rows(clock_audit)
  stopifnot(all(coverage$input_time_matches_calendar))
  coverage<-coverage %>% group_by(MPA,Source) %>% mutate(all_available_responses_min_before=min(n_before),
    all_available_responses_min_strict_before=min(n_strict_before)) %>% ungroup()
  write_out(coverage,'coverage');write_out(clock_audit,'clock_audit')
  scenarios<-list(retained_primary=primary);scenarios$retained_primary$HAC_SE<-baseline_hac
  candidates<-list();routefits<-list();decisions<-list()
  # Routing changes use the same supplied calendar clock as the retained PISCO fits.
  # Date boundary is a separate labeled sensitivity; no observations are removed.
  for(boundary in c('legacy','strict')) {
    for(i in which(primary$Source=='PISCO')) {
      row<-primary[i,];cv<-coverage[coverage$key==row$key,];d<-rebuilt[[row$key]]$data
      if(boundary=='strict')d$BA<-factor(ifelse(d$year<row$MPA_Start,'Before','After'),c('Before','After'))
      d$time.model<-d$input_time
      nb<-sum(d$BA=='Before');na<-sum(d$BA=='After')
      pt<-if(boundary=='legacy')cv$pretrend_p else cv$strict_pretrend_p
      if(nb<3 || na<3 || !is.finite(pt) || pt<.05)next
      fits<-list();aics<-numeric()
      for(model in c('Step','Linear','Asymptotic','Sigmoid')) {
        f<-safe(fit_one(d,model),paste('routing',boundary,row$key,model))
        if(is.null(f))next
        cc<-safe(contrast(f,model,year=d$year),paste('contrast',boundary,row$key,model));if(is.null(cc))next
        k<-attr(logLik(f),'df');n<-nrow(d);aic<-if(n>k+1) AIC(f)+2*k*(k+1)/(n-k-1) else Inf
        if(!is.finite(aic))next
        fits[[model]]<-cc;aics[model]<-aic
        candidates[[length(candidates)+1L]]<-cbind(key=row$key,boundary,cc,AICc=aic,n_before=nb,n_after=na)
      }
      if(!length(aics))stop('No valid candidates: ',row$key)
      selected<-names(which.min(aics));routefits[[paste(boundary,row$key)]]<-fits[[selected]]
    }
    for(minimum in c(3L,4L,5L)) {
      label<-paste0('PISCO_',boundary,'_min',minimum);x<-scenarios$retained_primary
      for(i in which(x$Source=='PISCO')) {
        cv<-coverage[coverage$key==x$key[i],];nb<-if(boundary=='legacy')cv$n_before else cv$n_strict_before
        pt<-if(boundary=='legacy')cv$pretrend_p else cv$strict_pretrend_p
        eligible<-nb>=minimum && cv$n_after>=3 && is.finite(pt) && pt>=.05
        if(eligible) {z<-routefits[[paste(boundary,x$key[i])]];stopifnot(!is.null(z));x[i,c('Mean','SE','HAC_SE','Model')]<-z[1,c('Mean','SE','HAC_SE','Model')];x$AnalysisType[i]<-'pBACIPS_sensitivity'}
        decisions[[length(decisions)+1L]]<-data.frame(scenario=label,key=x$key[i],n_before=nb,pretrend_p=pt,
          rerouted=eligible,selected_model=x$Model[i],reason=if(nb<minimum)'Insufficient before-years' else if(!is.finite(pt))'Before-trend not estimable (constant baseline)' else if(pt<.05)'Before-trend screen' else 'Eligible')
      }
      scenarios[[label]]<-x
    }
  }
  # Separate stricter eligibility from relabeling the designation-year observation.
  x<-scenarios$retained_primary
  for(i in which(x$Source=='PISCO')) {
    cv<-coverage[coverage$key==x$key[i],]
    if(cv$n_strict_before>=3 && cv$n_after>=3 && is.finite(cv$strict_pretrend_p) && cv$strict_pretrend_p>=.05) {
      z<-routefits[[paste('legacy',x$key[i])]]
      if(is.null(z))stop('Legacy fit unavailable for strict-eligibility check')
      x[i,c('Mean','SE','HAC_SE','Model')]<-z[1,c('Mean','SE','HAC_SE','Model')]
      x$AnalysisType[i]<-'pBACIPS_sensitivity'
    }
  }
  scenarios$PISCO_strict_eligibility_legacy_boundary<-x
  # Requiring every available species/response to meet the rule reroutes no PISCO histories.
  stopifnot(all(coverage$all_available_responses_min_before[coverage$Source=='PISCO']<3))
  scenarios$PISCO_all_responses_min3<-scenarios$retained_primary
  # Clock-only comparisons hold the retained model family fixed for all histories.
  for(label in c('calendar_fixed11','native_observed_window','calendar_observed_window')) {
    x<-scenarios$retained_primary
    for(i in which(x$Model %in% c('Linear','Sigmoid','Asymptotic'))) {
      row<-primary[i,];r<-rebuilt[[row$key]];d<-r$data;au<-clock_audit[clock_audit$key==row$key,]
      calendar<-label!='native_observed_window'
      if(calendar)d$time.model<-pmax(d$year-row$MPA_Start,0)
      f<-if(calendar && !isTRUE(all.equal(d$time.model,r$data$time.model))) safe(fit_one(d,row$Model),paste(label,row$key)) else r$fit
      if(is.null(f))stop('Clock refit failed, no silent omission: ',row$key)
      ends<-if(label=='calendar_fixed11') c(0,11) else if(calendar) c(au$matched_start_year,au$matched_end_year)-row$MPA_Start else c(au$native_start,au$native_end)
      cc<-contrast(f,row$Model,ends,d$year);x[i,c('Mean','SE','HAC_SE')]<-cc[1,c('Mean','SE','HAC_SE')]
    }
    scenarios[[label]]<-x
  }
  # Label-only correction: exactly the same Campus Point estimator, no numerical edit.
  cp<-'Campus Point SMCA|P. interruptus|LTER|Bio'
  x<-scenarios$retained_primary;x$AnalysisType[x$key==cp]<-'CI';scenarios$Campus_label_only<-x
  # Quantify the alternative mean-level CI estimand, without changing primary policy.
  x<-scenarios$retained_primary;i<-match(cp,x$key);d<-rebuilt[[cp]]$data
  cc<-contrast(lm(lnDiff~1,d), 'Mean',year=d$year);x[i,c('Mean','SE','HAC_SE','Model')]<-cc
  x$AnalysisType[i]<-'CI';scenarios$Campus_CI_mean<-x
  paths<-function(x,scenario,uncertainty) {
    x$vi<-if(uncertainty=='IID')x$SE^2 else x$HAC_SE^2
    pooled<-x%>%group_by(MPA,Taxa,Resp)%>%summarise(Mean=mean(Mean),vi=sum(vi)/n()^2,.groups='drop')
    defs<-data.frame(xTaxa=c('P. interruptus','P. interruptus','S. pulcher','S. pulcher','S. purpuratus','M. franciscanus'),
      xResp=c('Bio','Den','Bio','Den','Den','Den'),yTaxa=c(rep('S. purpuratus',4),rep('M. pyrifera',2)),yResp=c('Bio','Den','Bio','Den','Bio','Bio'))
    bind_rows(lapply(seq_len(nrow(defs)),function(i){z<-defs[i,];a<-pooled[pooled$Taxa==z$xTaxa&pooled$Resp==z$xResp,];b<-pooled[pooled$Taxa==z$yTaxa&pooled$Resp==z$yResp,]
      d<-merge(a,b,by='MPA',suffixes=c('_x','_y'));m<-rma(yi=Mean_y,vi=vi_y,mods=~Mean_x,data=d,method='REML',test='knha')
      data.frame(scenario,uncertainty,panel=letters[i],k=m$k,mpas=paste(sort(d$MPA),collapse=';'),slope=as.numeric(coef(m)[2]),ci_low=m$ci.lb[2],ci_high=m$ci.ub[2],p=m$pval[2])}))
  }
  pooled<-list();pathways<-list();series<-list()
  for(label in names(scenarios)) {
    x<-scenarios[[label]];stopifnot(identical(x$key,primary$key),identical(x$N,primary$N))
    series[[label]]<-cbind(scenario=label,x)
    for(u in c('IID','calendar_HAC')) {
      se<-if(u=='IID')x$SE else x$HAC_SE
      z<-uncertainty_joint(x,se,label);z$uncertainty<-u
      pooled[[paste(label,u)]]<-z;pathways[[paste(label,u)]]<-paths(x,label,u)
    }
  }
  pooled<-bind_rows(pooled);pathways<-bind_rows(pathways)
  base<-filter(pooled,scenario=='retained_primary')%>%select(Taxa,Resp,uncertainty,base_estimate=estimate,base_q=q)
  pooled<-left_join(pooled,base,by=c('Taxa','Resp','uncertainty'))%>%mutate(direction_changed=sign(estimate)!=sign(base_estimate),bh_changed=(q<.05)!=(base_q<.05))
  bpath<-filter(pathways,scenario=='retained_primary')%>%select(panel,uncertainty,base_slope=slope,base_p=p,base_mpas=mpas)
  pathways<-left_join(pathways,bpath,by=c('panel','uncertainty'))%>%mutate(direction_changed=sign(slope)!=sign(base_slope),nominal_changed=(p<.05)!=(base_p<.05))
  stopifnot(all(pathways$mpas==pathways$base_mpas))
  table2<-read.csv(file.path(root,'tables/table_02_meta_analysis.csv'));table2$Resp<-ifelse(table2$Response=='Biomass','Bio','Den')
  j<-merge(filter(pooled,scenario=='retained_primary',uncertainty=='IID'),table2,by=c('Taxa','Resp'))
  stopifnot(nrow(j)==9,max(abs(j$estimate-j$Estimate))<1e-7,max(abs(j$q-j$pval_fdr))<1e-7)
  table3<-read.csv(file.path(root,'tables/table_03_cross_taxa_meta_regression.csv'))
  # Panel order follows Table 3's six named trophic relationships.
  b<-filter(pathways,scenario=='retained_primary',uncertainty=='IID')
  pred<-c('P. interruptus biomass','P. interruptus density','S. pulcher biomass','S. pulcher density','S. purpuratus density','M. franciscanus density')
  resp<-c('S. purpuratus biomass','S. purpuratus density','S. purpuratus biomass','S. purpuratus density','M. pyrifera biomass','M. pyrifera biomass')
  for(i in 1:6) {
    t<-table3[table3$Predictor==pred[i] & table3$Response_Var==resp[i],];z<-b[b$panel==letters[i],]
    stopifnot(nrow(t)==1,z$k==t$k,abs(round(z$slope,4)-t$Slope_Est)<1e-8,abs(z$p-t$Slope_pval)<1e-8)
  }
  # The clock comparison below uses the same observed calendar endpoints in both fits.
  ca<-filter(pooled,scenario=='calendar_observed_window')
  na<-filter(pooled,scenario=='native_observed_window')
  dc<-merge(ca,na,by=c('Taxa','Resp','uncertainty'),suffixes=c('_calendar','_native'))
  dc$estimate_difference<-dc$estimate_calendar-dc$estimate_native
  dc$clock_direction_changed<-sign(dc$estimate_calendar)!=sign(dc$estimate_native)
  dc$clock_bh_changed<-(dc$q_calendar<.05)!=(dc$q_native<.05)
  write_out(dc,'matched_clock_difference')
  write_out(bind_rows(candidates),'candidates');write_out(bind_rows(decisions),'routing_decisions')
  write_out(bind_rows(series),'series');write_out(pooled,'pooled');write_out(pathways,'pathways');write_out(bind_rows(warnlog),'fit_log')
  stopifnot(identical(hashes,tools::md5sum(file.path(root,files))))
  write_out(data.frame(path=files,md5=unname(hashes)),'inputs')
  cat('PASS: 145 histories reconstructed, fixed annual inputs and pathway memberships; primary Table 2 reproduced.\n')
  print(filter(pooled,uncertainty=='IID')%>%group_by(scenario)%>%summarise(direction_changes=sum(direction_changed),bh_changes=sum(bh_changed)))
})
