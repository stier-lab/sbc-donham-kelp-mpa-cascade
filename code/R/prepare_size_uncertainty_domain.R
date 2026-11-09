# Recover the pre-pairing annual domain that defines the adaptive zero adjustment.
# Run in a separate R process; requires valid existing source biomass caches.
analysis_root <- normalizePath('.')
dr <- Sys.getenv('DONHAM_DATA_REPO',file.path(dirname(analysis_root),'sbc-donham-kelp-mpa-data'))
private <- file.path(analysis_root,'local_archive/primary_size_uncertainty_2026-10-30')
dir.create(private,recursive=TRUE,showWarnings=FALSE)
modules <- c('00_libraries.R','00c_analysis_constants.R','01_utils.R','02_pBACIPS_function.R',
             '03_data_import.R','04_pisco_processing.R','05_kfm_processing.R')
caches <- list.files(file.path(dr,'data/cache'),pattern='\\.rds$',full.names=TRUE)
required <- file.path(dr,'data/cache',c('vrg_panint_bootstrap.rds','pisco_urchin_bootstrap.rds',
  'kfm_urchin_bootstrap.rds','kfm_macro_bootstrap_sbc11215.rds'))
stopifnot(all(file.exists(required)))
old_hash <- tools::md5sum(caches)
setwd(dr)
here::i_am('code/R/run_all.R')
# Fail closed if a source processor tries to regenerate a missing/invalid cache.
saveRDS <- function(...)stop('Cache writes prohibited in size-domain reconstruction.')
for(m in modules){
  source(file.path(dr,'code/R',m))
  if(m=='01_utils.R'){
    if(should_force_bootstrap())stop('Forced cache regeneration is not allowed here.')
    bootstrap_biomass <- function(...)stop('Cache regeneration prohibited; rebuild sources separately.')
  }
}
library(data.table)
p <- as.data.table(All_PISCO.mpa)[y %in% c('STRPURAD','PANINT'),
  .(MPA=CA_MPA_Name_Short,year,status=site_status,value=biomass/PISCO_SWATH_AREA_M2,
    Taxa=ifelse(y=='PANINT','P. interruptus','S. purpuratus'),Source='PISCO')]
k <- as.data.table(kfm.ave.ave)[y %in% c('Strongylocentrotus purpuratus','STRPURAD'),
  .(MPA=CA_MPA_Name_Short,year,status,value=biomass,Taxa='S. purpuratus',Source='KFM')]
d <- rbind(p,k);d <- d[is.finite(value) & status %in% c('mpa','reference')]
d[,series_id:=paste(MPA,Taxa,Source,'Bio',sep='|')]
stopifnot(!anyDuplicated(d[,.(series_id,year,status)]),identical(old_hash,tools::md5sum(caches)))
fwrite(d,file.path(private,'full_annual_domain.csv'))
source_files <- file.path(dr,'data',c('ALL_sizefreq_2024.csv','Site_List_All.csv','MPAfeatures_subset.csv',
  'PISCO/master_site_table_Emilyedit.csv','PISCO/MLPA_kelpforest_swath_2024.csv',
  'PISCO/MLPA_kelpforest_fish_2024.csv','PISCO/spp_attribute_table_downloaded_9-13-22_SHSPUL.csv',
  'MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv',
  'MBON/SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv',
  'MBON/KFM_Macrocystis_RawData_1984-2023.csv','MBON/SBCMBON_kelp_forest_integrated_fish_20231022.csv'))
files <- c(file.path(dr,'code/R',modules),caches,source_files,
  file.path(analysis_root,'code/R/prepare_size_uncertainty_domain.R'))
fwrite(data.table(path=files,md5=unname(tools::md5sum(files))),file.path(private,'domain_inputs.csv'))
