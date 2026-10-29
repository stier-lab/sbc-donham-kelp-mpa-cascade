# Descriptive fixed-pair diagnostic; never changes the primary estimators.
local({
  suppressPackageStartupMessages({library(dplyr); library(tidyr); library(ggplot2)})
  source("code/R/satellite_siting_utils.R", local=TRUE)
  inputs <- c("data/harmonized/harmonized_landsat_annual_index.csv",
    "data/harmonized/harmonized_landsat_rr.csv", "data/harmonized/harmonized_site_metadata.csv",
    "data/sumstats_final.csv", "tables/table_s_historical_heatwave_pairs.csv",
    "code/R/46_satellite_siting.R", "code/R/satellite_siting_utils.R", "docs/SATELLITE_SITING_PLAN.md")
  stopifnot(all(file.exists(inputs)))
  protected <- c(list.files("data/harmonized",full.names=TRUE), "data/sumstats_final.csv",
    list.files("plots",pattern="^fig_0|^fig_kelp_resilience",full.names=TRUE))
  protected <- protected[!dir.exists(protected)]
  hashes <- tools::md5sum(protected)
  save <- function(d, suffix) write.csv(d,paste0("tables/table_s_siting_",suffix,".csv"),row.names=FALSE,na="")
  primary <- read.csv(inputs[4]); mpas <- sort(unique(primary$MPA[primary$Source=="Landsat"]))
  annual <- read.csv(inputs[1]); meta <- read.csv(inputs[3])
  annual <- annual[annual$CA_MPA_Name_Short %in% mpas,]
  stopifnot(!anyDuplicated(annual[c("CA_MPA_Name_Short","year","status")]))
  d <- annual %>% select(CA_MPA_Name_Short,year,status,index,log_offset) %>%
    pivot_wider(names_from=status,values_from=c(index,log_offset)) %>%
    rename(MPA=CA_MPA_Name_Short,inside=index_mpa,reference=index_reference) %>%
    left_join(meta %>% select(CA_MPA_Name_Short,MPA_Start,Location),by=c("MPA"="CA_MPA_Name_Short")) %>%
    arrange(MPA,year)
  stopifnot(!anyNA(d),all(d$inside>=0 & d$reference>=0),all(d$log_offset_mpa>0),
    max(abs(d$log_offset_mpa-d$log_offset_reference))<1e-10)
  d$contrast <- log(d$inside+d$log_offset_mpa)-log(d$reference+d$log_offset_reference)
  rr <- read.csv(inputs[2]) %>% filter(CA_MPA_Name_Short %in% mpas) %>% select(CA_MPA_Name_Short,year,lnDiff)
  guard <- left_join(d,rr,by=c("MPA"="CA_MPA_Name_Short","year"))
  stopifnot(!anyNA(guard$lnDiff),max(abs(guard$lnDiff-guard$contrast))<1e-8)
  save(d,"annual")
  metrics <- windows <- list()
  for(m in mpas) {
    x <- d[d$MPA==m,]; start <- unique(x$MPA_Start); stopifnot(length(start)==1)
    for(window in c("full","last10","last5")) {
      first <- if(window=="full") 1984 else start-if(window=="last10") 10 else 5
      z <- x[x$year>=first & x$year<start,]
      enough <- nrow(z)>=4
      metrics[[length(metrics)+1]] <- data.frame(MPA=m,window,start,first,
        n_years=nrow(z),years=paste(z$year,collapse=";"),
        before_log=mean(z$contrast),before_index=mean(z$inside-z$reference),
        inside_zero=sum(z$inside==0),reference_zero=sum(z$reference==0),
        slope=if(enough) unname(coef(lm(contrast~year,z))[2]) else NA_real_,
        tracking=if(enough) siting_tracking(z$inside,z$reference) else NA_real_,
        detrended_tracking=if(enough) siting_tracking(resid(lm(inside~year,z)),resid(lm(reference~year,z))) else NA_real_)
    }
    for(w in c(5L,3L,7L)) for(a in c("actual","placebo_near","placebo_far")) {
      anchor <- start-switch(a,actual=0,placebo_near=w+1,placebo_far=2*w+1)
      if(a!="actual") stopifnot(anchor+w<start)
      windows[[length(windows)+1]] <- cbind(MPA=m,start,scenario=a,siting_window(x,anchor,w))
    }
  }
  metrics <- bind_rows(metrics); windows <- bind_rows(windows)
  windows <- windows %>% group_by(MPA,width) %>% mutate(fixed_cohort=all(eligible),
    positive_fixed_cohort=all(eligible & positive_means)) %>% ungroup()
  save(metrics,"preperiod"); save(windows,"windows")
  summary <- windows %>% filter(fixed_cohort) %>% group_by(width,scenario) %>%
    summarise(n_pairs=n(),members=paste(sort(MPA),collapse=";"),
      mean_log=mean(change_log),median_log=median(change_log),positive_log=sum(change_log>0),
      mean_index=mean(change_index),positive_index=sum(change_index>0),.groups="drop")
  save(summary,"summary")
  positive <- windows %>% filter(positive_fixed_cohort) %>% group_by(width,scenario) %>%
    summarise(n_pairs=n(),members=paste(sort(MPA),collapse=";"),mean_period_log=mean(change_period_log),.groups="drop")
  save(positive,"positive_mean_summary")
  comparison <- windows %>% filter(fixed_cohort) %>% select(MPA,width,scenario,change_log) %>%
    pivot_wider(names_from=scenario,values_from=change_log) %>%
    mutate(actual_minus_near=actual-placebo_near,actual_minus_far=actual-placebo_far)
  save(comparison,"paired_comparison")
  historical <- read.csv(inputs[5]) %>% filter(specification=="primary")
  save(historical,"historical_existing")
  save(data.frame(path=inputs,md5=unname(tools::md5sum(inputs))),"inputs")
  # Unshifted index trajectories reveal zeros and different footprints; not kg.
  long <- d %>% pivot_longer(c(inside,reference),names_to="status",values_to="index")
  long <- long %>% group_by(MPA,status) %>% mutate(segment=cumsum(c(TRUE,diff(year)>1))) %>% ungroup()
  p <- ggplot(long,aes(year,index,color=status,group=interaction(status,segment)))+
    geom_line(linewidth=.35)+geom_point(size=.65)+
    geom_vline(data=distinct(d,MPA,MPA_Start),aes(xintercept=MPA_Start),inherit.aes=FALSE,linetype=2,color="grey45")+
    facet_wrap(~MPA,ncol=4)+scale_color_manual(values=c(inside="#007F73",reference="#BA4A58"),labels=c("Inside","Reference"))+
    labs(x="Year",y="Annual kelp index (not biomass units)",color=NULL)+theme_bw(base_size=9)+theme(legend.position="bottom")
  ggsave("plots/fig_s_siting_histories.png",p,width=13,height=11,dpi=160)
  ggsave("plots/fig_s_siting_histories.pdf",p,width=13,height=11)
  logdata <- d %>% group_by(MPA) %>% mutate(segment=cumsum(c(TRUE,diff(year)>1))) %>% ungroup()
  q <- ggplot(logdata,aes(year,contrast,group=segment))+geom_hline(yintercept=0,color="grey70")+
    geom_line(linewidth=.35,color="#007F73")+geom_point(size=.65,color="#007F73")+
    geom_vline(data=distinct(d,MPA,MPA_Start),aes(xintercept=MPA_Start),inherit.aes=FALSE,linetype=2,color="grey45")+
    facet_wrap(~MPA,ncol=4)+labs(x="Year",y="Log inside/reference contrast (recorded offset)")+theme_bw(base_size=9)
  ggsave("plots/fig_s_siting_log_histories.png",q,width=13,height=11,dpi=160)
  ggsave("plots/fig_s_siting_log_histories.pdf",q,width=13,height=11)
  stopifnot(identical(hashes,tools::md5sum(protected)))
  save(data.frame(status="complete",n_primary_pairs=length(mpas),primary_inputs_unchanged=TRUE,
    inference="descriptive only; no independent-year or randomization P values"),"status")
  print(summary)
})
