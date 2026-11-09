# Coverage-selected, fixed-station sheephead comparisons. Primary files are read-only.
# Run from repository root: Rscript code/R/59_pisco_fixed_stations.R
local({
  suppressPackageStartupMessages({library(here);library(dplyr);library(tidyr);library(data.table);library(metafor)})
  root <- here::here()
  data_root <- Sys.getenv("MPA_DATA_REPO", file.path(dirname(root), "sbc-donham-kelp-mpa-data"))
  source(file.path(root,"code/R/primary_uncertainty_utils.R"),local=TRUE)
  defs <- uncertainty_definitions(root)
  out <- function(x,n) write.csv(x,file.path(root,"tables",paste0("table_s_pisco_fixed_",n,".csv")),row.names=FALSE)
  input <- c(file.path(root,c("data/sumstats_final.csv","data/harmonized/harmonized_raw_responses.csv",
    "data/harmonized/harmonized_response_ratios.csv","tables/table_s_selection_station_years.csv",
    "tables/table_s_selection_fish_invalid_size_years.csv","tables/table_s_routing_clock_series.csv")),
    file.path(data_root,c("data/PISCO/MLPA_kelpforest_fish_2024.csv","data/PISCO/spp_attribute_table_downloaded_9-13-22_SHSPUL.csv")))
  input <- c(input,file.path(root,c("code/R/59_pisco_fixed_stations.R","code/R/primary_uncertainty_utils.R","code/R/02_pBACIPS_function.R","code/R/08_effect_sizes.R")))
  input<-c(input,list.files(file.path(root,"data/external/pisco_policy_20261031"),full.names=TRUE),file.path(root,"code/R/pisco_nonlinear_candidates.R"),file.path(root,"tables/table_s_primary_hac_histories.csv"))
  hashes <- tools::md5sum(input)
  stopifnot(!anyNA(hashes))
  primary <- read.csv(input[1]);primary$key<-uncertainty_key(primary)
  historical <- read.csv(file.path(root,"data/external/pisco_policy_20261031/primary_before.csv"))
  historical$key <- uncertainty_key(historical)
  rr <- read.csv(file.path(root,"data/external/pisco_policy_20261031/ratios_before.csv"))
  annual <- rr %>% filter(source=="PISCO",y=="Semicossyphus pulcher") %>%
    transmute(MPA=CA_MPA_Name_Short,Source=source,Taxa="S. pulcher",Resp=resp,year,lnDiff)
  source(file.path(root,"code/R/pisco_nonlinear_candidates.R"),local=TRUE)
  mpa_names <- c("Campus Point SMCA","Naples SMCA","Point Vicente SMCA")
  target <- historical %>% filter(MPA %in% mpa_names,Source=="PISCO",Taxa=="S. pulcher")
  stopifnot(nrow(target)==6,all(target$MPA_Start==2012))
  roster <- read.csv(file.path(root,"data/external/pisco_policy_20261031/table_s_selection_station_years.csv")) %>% filter(branch=="fish",MPA %in% mpa_names)
  bad <- read.csv(input[5])
  # Reproduce the source's species-attribute eligibility before forming surveyed events.
  att <- fread(input[8]); eligible_codes <- att[!is.na(WL_input_length),pisco_classcode]
  fish <- fread(input[7]);fish<-fish[campus %in% c("UCSB","VRG") & level!="CAN" &
    site %in% roster$site & classcode %in% eligible_codes]
  cells <- unique(fish[,.(site,year,zone,transect)])
  sp <- fish[classcode=="SPUL"]
  sp[,biomass:=count*.0144*fish_tl^3.04]
  sums <- sp[,.(Den=sum(count)/60,Bio=sum(biomass)/60),by=.(site,year,zone,transect)]
  cells <- merge(cells,sums,by=c("site","year","zone","transect"),all.x=TRUE)
  # Exact source reproduction: numeric NAs are zero-filled after aggregation.
  # Invalid size/count site-years are excluded from fixed-station eligibility below.
  cells[is.na(Den),Den:=0];cells[is.na(Bio),Bio:=0]
  zones <- cells[,.(Den=mean(Den),Bio=mean(Bio)),by=.(site,year,zone)]
  sites <- as.data.frame(zones[,.(Den=mean(Den),Bio=mean(Bio)),by=.(site,year)])
  sites <- inner_join(sites,roster[,c("MPA","site","status")],by="site",relationship="many-to-many")
  mix <- sites %>% group_by(MPA,year,status) %>% summarise(Den=mean(Den),Bio=mean(Bio),.groups="drop") %>%
    pivot_longer(c(Den,Bio),names_to="Resp",values_to="value")
  raw <- read.csv(file.path(root,"data/external/pisco_policy_20261031/raw_before.csv")) %>% filter(source=="PISCO",taxon_name=="Semicossyphus pulcher",CA_MPA_Name_Short %in% mpa_names) %>%
    transmute(MPA=CA_MPA_Name_Short,year,status,Resp=resp,original=value)
  validation <- inner_join(mix,raw,by=c("MPA","year","status","Resp"))
  validation$difference<-validation$value-validation$original
  stopifnot(nrow(validation)==nrow(raw),max(abs(validation$difference))<1e-7)
  out(validation,"reconstruction")
  # One correction per MPA/response, derived from the complete primary mixture;
  # hold it fixed across station subsets and year subsets to isolate membership/model.
  corrections <- mix %>% group_by(MPA,Resp) %>% summarise(correction=min(value[value>0])/2,.groups="drop")
  ratio <- function(x) x %>% left_join(corrections,by=c("MPA","Resp")) %>%
    pivot_wider(names_from=status,values_from=value) %>% filter(is.finite(mpa),is.finite(reference)) %>%
    mutate(lnDiff=log((mpa+correction)/(reference+correction)))
  mix_rr <- ratio(mix)
  check <- inner_join(mix_rr,annual %>% filter(Source=="PISCO",Taxa=="S. pulcher",MPA %in% mpa_names),by=c("MPA","year","Resp"),suffix=c("_rebuilt","_primary"))
  stopifnot(nrow(check)==sum(target$N),max(abs(check$lnDiff_rebuilt-check$lnDiff_primary))<1e-7)
  # Coverage choice is independent of sheephead responses. Enumerate station groups;
  # maximize station count, then common-year count, then alphabetical station key.
  # Do not pool overlapping pairs as independent observations.
  station_years <- roster %>% select(MPA,site,status,years) %>% separate_rows(years,sep=";") %>%
    mutate(year=as.integer(years)) %>% select(-years) %>% anti_join(bad,by=c("site","year"))
  subsets <- function(x) unlist(lapply(seq_along(x),function(n) combn(sort(x),n,simplify=FALSE)),recursive=FALSE)
  groups<-list(); choices<-list(); matched<-list(); membership<-list()
  for(m in mpa_names) {
    sy<-filter(station_years,MPA==m)
    for(p in subsets(unique(sy$site[sy$status=="mpa"]))) for(r in subsets(unique(sy$site[sy$status=="reference"]))) {
      ss<-c(p,r);yy<-sort(Reduce(intersect,lapply(ss,function(z)sy$year[sy$site==z])))
      if(sum(yy>2012)<4)next
      groups[[length(groups)+1]]<-data.frame(MPA=m,protected=paste(p,collapse=";"),reference=paste(r,collapse=";"),
        n_stations=length(ss),n_years=length(yy),strict_before=sum(yy<2012),inclusive_before=sum(yy<=2012),
        after=sum(yy>2012),years=paste(yy,collapse=";"))
    }
  }
  groups<-bind_rows(groups);out(groups,"feasible_groups")
  for(boundary in c("inclusive","strict")) for(minimum in c(3,4)) for(m in mpa_names) {
    z<-filter(groups,MPA==m);nb<-z[[paste0(boundary,"_before")]]
    z<-z[nb>=minimum,];z<-arrange(z,desc(n_stations),desc(n_years),protected,reference)
    label<-paste(boundary,minimum,sep="_min")
    if(!nrow(z)) {choices[[length(choices)+1]]<-data.frame(scenario=label,MPA=m,eligible=FALSE);next}
    z<-z[1,];choices[[length(choices)+1]]<-cbind(scenario=label,eligible=TRUE,z)
    yy<-as.integer(strsplit(z$years,";",fixed=TRUE)[[1]]);ss<-strsplit(paste(z$protected,z$reference,sep=";"),";",fixed=TRUE)[[1]]
    fixed<-filter(sites,MPA==m,site %in% ss,year %in% yy)
    stopifnot(nrow(fixed)==length(ss)*length(yy),!anyDuplicated(fixed[c("site","year")]))
    membership[[length(membership)+1]]<-mutate(fixed,scenario=label)
    f<-fixed %>% group_by(MPA,year,status) %>% summarise(Den=mean(Den),Bio=mean(Bio),.groups="drop") %>%
      pivot_longer(c(Den,Bio),names_to="Resp",values_to="value") %>% ratio() %>% mutate(membership="fixed")
    a<-filter(mix_rr,MPA==m,year %in% yy) %>% mutate(membership="mixture_matched_years")
    matched[[length(matched)+1]]<-bind_rows(f,a) %>% mutate(scenario=label,boundary=boundary,minimum=minimum)
  }
  choices<-bind_rows(choices);dats<-bind_rows(matched)
  out(choices,"coverage");out(bind_rows(membership),"station_values");out(dats,"annual")
  warnlog<-list()
  safe<-function(expr,label) {
    w<-character();z<-tryCatch(withCallingHandlers(expr,warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart("muffleWarning")}),error=identity)
    warnlog[[length(warnlog)+1]]<<-data.frame(label=label,warnings=paste(unique(w),collapse=";"),error=if(inherits(z,"error"))conditionMessage(z) else "")
    if(inherits(z,"error"))NULL else z
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
      Sigmoid=pisco_nonlinear_candidate(d,"Sigmoid")$fit,
      Asymptotic=pisco_nonlinear_candidate(d,"Asymptotic")$fit)
    if(is.null(z))stop('No fitted model')
    if(model %in% c('Sigmoid','Asymptotic') && (!inherits(z,'nls')||!isTRUE(z$convInfo$isConv)))stop('Nonlinear fallback rejected')
    if(any(!is.finite(coef(z)))||any(!is.finite(vcov(z))))stop('Invalid coefficients/covariance')
    z
  }
  paths<-function(x,scenario,uncertainty) {
    x$vi<-if(uncertainty=='IID')x$SE^2 else x$HAC_SE^2
    pooled<-x%>%group_by(MPA,Taxa,Resp)%>%summarise(Mean=mean(Mean),vi=sum(vi)/n()^2,.groups='drop')
    defs<-data.frame(xTaxa=c('P. interruptus','P. interruptus','S. pulcher','S. pulcher','S. purpuratus','M. franciscanus'),
      xResp=c('Bio','Den','Bio','Den','Den','Den'),yTaxa=c(rep('S. purpuratus',4),rep('M. pyrifera',2)),yResp=c('Bio','Den','Bio','Den','Bio','Bio'))
    bind_rows(lapply(seq_len(nrow(defs)),function(i){z<-defs[i,];a<-pooled[pooled$Taxa==z$xTaxa&pooled$Resp==z$xResp,];b<-pooled[pooled$Taxa==z$yTaxa&pooled$Resp==z$yResp,]
      d<-merge(a,b,by='MPA',suffixes=c('_x','_y'));m<-rma(yi=Mean_y,vi=vi_y,mods=~Mean_x,data=d,method='REML',test='knha')
      data.frame(scenario,uncertainty,panel=letters[i],k=m$k,mpas=paste(sort(d$MPA),collapse=';'),slope=as.numeric(coef(m)[2]),ci_low=m$ci.lb[2],ci_high=m$ci.ub[2],p=m$pval[2])}))
  }

  fits<-list();candidates<-list();decisions<-list()
  for(idx in split(seq_len(nrow(dats)),interaction(dats$scenario,dats$MPA,dats$Resp,dats$membership,drop=TRUE))) {
    d<-dats[idx,];d<-d[order(d$year),];boundary<-d$boundary[1]
    d$BA<-factor(ifelse(if(boundary=="strict")d$year<2012 else d$year<=2012,"Before","After"),c("Before","After"))
    d$time.model<-pmax(d$year-2012,0)
    ends<-defs$observed_contrast_times(d)
    pre<-d[d$BA=="Before",];pt<-if(nrow(pre)>=3&&sd(pre$lnDiff)>0)summary(lm(lnDiff~year,pre))$coefficients[2,4] else NA_real_
    id<-d[1,c("scenario","MPA","Resp","membership")]
    lf<-fit_one(d,"Linear");ci_model<-if(summary(lf)$coefficients[2,4]<=.05)"Linear" else "Mean"
    # CI is reselected on the same observations; also retain original family to isolate membership.
    original<-target$Model[target$MPA==id$MPA & target$Resp==id$Resp]
    for(method in c("CI","original_family")) {
      model<-if(method=="CI")ci_model else original
      f<-fit_one(d,model);cc<-contrast(f,model,ends,d$year)
      fits[[length(fits)+1]]<-cbind(id,method,cc,N=nrow(d),n_before=nrow(pre),pretrend_p=pt,
        first_after=min(d$year[d$BA=="After"]),last_contrast=2012+ends[2])
    }
    eligible<-is.finite(pt)&&pt>=.05
    decisions[[length(decisions)+1]]<-cbind(id,pretrend_p=pt,eligible,reason=if(eligible)"Eligible" else "Baseline trend or unestimable baseline")
    if(!eligible)next
    aics<-numeric();cfit<-list()
    for(model in c("Step","Linear","Asymptotic","Sigmoid")) {
      f<-safe(fit_one(d,model),paste(unlist(id),model,collapse="|"));if(is.null(f))next
      cc<-safe(contrast(f,model,ends,d$year),paste(unlist(id),model,"contrast",collapse="|"));if(is.null(cc))next
      k<-attr(logLik(f),"df");n<-nrow(d);aic<-if(n>k+1) AIC(f)+2*k*(k+1)/(n-k-1) else Inf
      if(!is.finite(aic))next
      aics[model]<-aic;cfit[[model]]<-cc
      candidates[[length(candidates)+1]]<-cbind(id,cc,AICc=aic)
    }
    stopifnot(length(aics)>0);cc<-cfit[[names(which.min(aics))]]
    fits[[length(fits)+1]]<-cbind(id,method="pBACIPS",cc,N=nrow(d),n_before=nrow(pre),pretrend_p=pt,
      first_after=min(d$year[d$BA=="After"]),last_contrast=2012+ends[2])
  }
  fits<-bind_rows(fits);out(fits,"fits");out(bind_rows(candidates),"candidates");out(bind_rows(decisions),"decisions")
  # Replace only reconstructed histories. All other primary histories stay unchanged.
  base<-primary
  h<-read.csv(file.path(root,"tables/table_s_primary_hac_histories.csv")) %>% filter(scenario=="automatic")
  base$HAC_SE<-h$synthesis_se[match(base$key,uncertainty_key(h))]
  historical_base<-read.csv(file.path(root,"data/external/pisco_policy_20261031/routing_series_before.csv")) %>% filter(scenario=="retained_primary") %>% select(-scenario)
  stopifnot(identical(base$key,primary$key))
  pools<-list();pathways<-list()
  for(sc in unique(fits$scenario)) for(mem in c("fixed","mixture_matched_years")) for(method in c("CI","pBACIPS")) {
    x<-base
    # Histories that fail the alternative minimum retain their original CI comparison.
    ii<-which(x$key %in% target$key);jj<-match(x$key[ii],historical_base$key)
    x[ii,c("Mean","SE","HAC_SE","Model","N")]<-historical_base[jj,c("Mean","SE","HAC_SE","Model","N")]
    z<-filter(fits,scenario==sc,membership==mem,method==.env$method)
    # Failed baseline screens retain CI on the same matched observations, not the original mixture.
    fallback<-filter(fits,scenario==sc,membership==mem,method=="CI") %>% anti_join(z,by=c("MPA","Resp"))
    z<-bind_rows(z,fallback)
    for(i in seq_len(nrow(z))) {
      j<-which(x$Source=="PISCO" & x$Taxa=="S. pulcher" & x$MPA==z$MPA[i] & x$Resp==z$Resp[i])
      stopifnot(length(j)==1);x[j,c("Mean","SE","HAC_SE","Model","N")]<-z[i,c("Mean","SE","HAC_SE","Model","N")]
    }
    label<-paste(sc,mem,method,sep="|")
    for(u in c("IID","calendar_HAC")) {
      se<-if(u=="IID")x$SE else x$HAC_SE
      p<-uncertainty_joint(x,se,label);p$uncertainty<-u;pools[[length(pools)+1]]<-p
      pathways[[length(pathways)+1]]<-paths(x,label,u)
    }
  }
  out(bind_rows(pools),"pooled");out(bind_rows(pathways),"pathways");out(bind_rows(warnlog),"fit_log")
  stopifnot(identical(hashes,tools::md5sum(input)))
  out(data.frame(path=input,md5=unname(hashes)),"inputs")
  cat("PASS: source annual values and ratios reproduced; fixed station/year groups and models fitted; primary inputs unchanged.\n")
  print(choices)
  print(filter(fits,membership=="fixed",method %in% c("CI","pBACIPS")))
})
