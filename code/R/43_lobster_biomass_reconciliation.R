# Read-only decomposition of secondary size-pool versus primary lobster biomass.
(function() {
  suppressPackageStartupMessages({library(dplyr); library(tidyr)})
  root <- here::here()
  source(file.path(root,"code/R/01_utils.R"),local=TRUE)
  source(file.path(root,"code/R/lobster_biomass_utils.R"),local=TRUE)
  protected <- file.path(root,c("data/harmonized/harmonized_raw_responses.csv",
    "data/harmonized/harmonized_response_ratios.csv","data/sumstats_final.csv",
    "tables/table_02_meta_analysis.csv","tables/table_03_cross_taxa_meta_regression.csv"))
  before <- tools::md5sum(protected)
  write <- function(d,stem) write.csv(d,file.path(root,"tables",paste0("table_s_lobster_reconciliation_",stem,".csv")),row.names=FALSE)
  paths <- file.path(root,c("data/ALL_sizefreq_2024.csv","data/PISCO/MLPA_kelpforest_swath_2024.csv",
    "data/PISCO/master_site_table_Emilyedit.csv"))
  paths <- c(paths,Sys.getenv("DONHAM_PRIMARY_LOBSTER_CACHE",unset=
    file.path(root,"data/cache/vrg_panint_bootstrap.rds")))
  if (!all(file.exists(paths))) {
    write(data.frame(status="unavailable",reason="Raw sizes, swaths, mapping and primary bootstrap cache required"),"status")
    return(invisible(NULL))
  }
  sizes <- read.csv(paths[1]) %>% filter(classcode=="PANINT",is.finite(size),size>0,is.finite(count),count>0)
  sites <- read.csv(paths[3]) %>% select(site,site_new,campus,CA_MPA_Name_Short,site_status,site_designation,BaselineRegion)
  duplicates <- sites %>% group_by(site,CA_MPA_Name_Short) %>% filter(n()>1) %>% ungroup()
  write(duplicates,"mapping_repeated_rows")
  # Reproduce the inherited primary first-row rule explicitly; do not recurate mappings.
  sites <- sites[!duplicated(sites[c("site","CA_MPA_Name_Short")]),]
  usable <- sites %>% filter(site_status %in% c("mpa","reference"),!site %in% EXCLUDED_REFERENCE_SITES)
  pool <- function(s) s %>% inner_join(usable %>% distinct(site,CA_MPA_Name_Short,site_status),by="site",relationship="many-to-many") %>%
    group_by(CA_MPA_Name_Short,year,status=site_status) %>% group_modify(~lobster_pool_moments(.x)) %>% ungroup()
  seasonal <- pool(filter(sizes,month>=5,month<=10)) %>% rename(seasonal_exact_mass=exact_mass,seasonal_rounded_mass=rounded_mass)
  all_months <- pool(sizes) %>% select(CA_MPA_Name_Short,year,status,all_month_mass=exact_mass)
  primary <- read.csv(protected[1]) %>% filter(source=="PISCO",taxon_name=="Panulirus interruptus",resp=="Bio") %>%
    select(CA_MPA_Name_Short,year,status,primary_biomass=value)
  secondary <- read.csv(file.path(root,"tables/table_s_lobster_size_coverage.csv")) %>%
    select(CA_MPA_Name_Short,year,status,density=value,secondary_biomass=all_biomass)
  audit <- inner_join(primary,secondary,by=c("CA_MPA_Name_Short","year","status")) %>%
    left_join(seasonal,by=c("CA_MPA_Name_Short","year","status")) %>%
    left_join(all_months,by=c("CA_MPA_Name_Short","year","status")) %>%
    mutate(seasonal_exact_biomass=ifelse(density==0,0,density*seasonal_exact_mass),
      seasonal_rounded_biomass=ifelse(density==0,0,density*seasonal_rounded_mass),
      all_month_biomass=ifelse(density==0,0,density*all_month_mass))
  # Reconstruct cached primary VRG annuals, retaining its zero-size and site exclusions.
  cached <- readRDS(paths[4]) %>% inner_join(usable,by="site",relationship="many-to-many")
  bad <- cached %>% filter(count>0, !is.finite(biomass)|biomass==0) %>% distinct(site,year)
  cached <- cached %>% anti_join(bad,by=c("site","year"))
  cached_annual <- cached %>% group_by(CA_MPA_Name_Short,year,status=site_status,site,zone,transect) %>%
    summarise(count=sum(count),biomass=sum(biomass),.groups="drop") %>%
    mutate(count=count/transect,biomass=biomass/transect) %>%
    group_by(CA_MPA_Name_Short,year,status,site) %>%
    summarise(density=mean(count)/60,biomass=mean(biomass)/60,.groups="drop") %>%
    group_by(CA_MPA_Name_Short,year,status) %>%
    summarise(cached_primary_biomass=mean(biomass),cached_primary_density=mean(density),
      n_primary_sites=n(),.groups="drop")
  audit <- left_join(audit,cached_annual,by=c("CA_MPA_Name_Short","year","status"))
  # Join exactly as the original primary pool does: site, renamed site and campus.
  primary_pools <- sizes %>% filter(!site %in% EXCLUDED_REFERENCE_SITES) %>%
    inner_join(sites,by=c("site","site_new","campus"),relationship="many-to-many") %>%
    group_by(CA_MPA_Name_Short,year,status=site_status) %>%
    group_modify(~lobster_pool_moments(.x)) %>% ungroup() %>%
    select(CA_MPA_Name_Short,year,status,primary_pool_mass=exact_mass,primary_pool_variance=mass_variance)
  audit <- audit %>% left_join(primary_pools,by=c("CA_MPA_Name_Short","year","status")) %>%
    mutate(same_count_pool_expectation=cached_primary_density*primary_pool_mass,
      cached_minus_primary=cached_primary_biomass-primary_biomass,
      secondary_to_primary=ifelse(primary_biomass>0,secondary_biomass/primary_biomass,NA_real_),
      density_to_primary=ifelse(cached_primary_density>0,density/cached_primary_density,NA_real_))
  write(audit,"annual")
  # Cached draw totals should approximate N times the pool's expected individual mass.
  expected <- cached %>% left_join(primary_pools,by=c("CA_MPA_Name_Short","year","site_status"="status")) %>%
    mutate(expected_biomass=count*primary_pool_mass,
      mc_se=sqrt(count*primary_pool_variance/1000),
      standardized_mc_error=ifelse(mc_se>0,(biomass-expected_biomass)/mc_se,NA_real_))
  write(expected %>% select(site,year,zone,CA_MPA_Name_Short,site_status,count,biomass,
    expected_biomass,mc_se,standardized_mc_error),"bootstrap_expectations")
  write(bad,"excluded_site_years")
  summary <- data.frame(status="completed_secondary_only",overlap_rows=nrow(audit),
    positive_primary_rows=sum(audit$primary_biomass>0),
    max_cached_primary_difference=max(abs(audit$cached_minus_primary),na.rm=TRUE),
    min_secondary_ratio=min(audit$secondary_to_primary,na.rm=TRUE),
    max_secondary_ratio=max(audit$secondary_to_primary,na.rm=TRUE),
    min_full_year_ratio=min(audit$all_month_biomass[audit$primary_biomass>0]/audit$primary_biomass[audit$primary_biomass>0],na.rm=TRUE),
    max_full_year_ratio=max(audit$all_month_biomass[audit$primary_biomass>0]/audit$primary_biomass[audit$primary_biomass>0],na.rm=TRUE),
    max_full_year_expectation_difference=max(abs(audit$all_month_biomass-audit$same_count_pool_expectation),na.rm=TRUE),
    max_zero_variance_pool_error=max(abs(expected$biomass-expected$expected_biomass)[
      is.finite(expected$mc_se)&expected$mc_se==0&is.finite(expected$expected_biomass)],na.rm=TRUE),
    max_absolute_mc_z=max(abs(expected$standardized_mc_error),na.rm=TRUE),
    primary_changed=!identical(before,tools::md5sum(protected)))
  write(summary,"status")
  write(data.frame(path=paths,md5=unname(tools::md5sum(paths))),"inputs")
  stopifnot(!summary$primary_changed,summary$max_cached_primary_difference<1e-10,
    summary$max_full_year_expectation_difference<1e-8,
    summary$max_zero_variance_pool_error<1e-8,
    !anyDuplicated(audit[c("CA_MPA_Name_Short","year","status")]))
  print(summary)
})()
