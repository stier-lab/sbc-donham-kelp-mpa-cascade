# Independent reconstruction of shared weights, annual offsets and retained fits.
suppressPackageStartupMessages(library(data.table))
source('code/R/primary_uncertainty_utils.R')
private <- 'local_archive/primary_size_uncertainty_2026-10-30'
manifest <- fread('tables/table_s_size_uncertainty_inputs.csv')
stopifnot(all(unname(tools::md5sum(manifest$path))==manifest$md5))
status <- fread('tables/table_s_size_uncertainty_status.csv')
d <- readRDS(file.path(private,paste0('draws_seed',status$seed,'.rds')))
ratios <- readRDS(file.path(private,paste0('pool_ratios_seed',status$seed,'.rds')))
pools <- fread('tables/table_s_size_uncertainty_pools.csv')
inventory <- fread('tables/table_s_primary_size_inventory.csv')
members <- fread(file.path(private,'complete_record_membership.csv'))
members <- merge(members,pools[,.(pool,series_id,year,status)],by=c('series_id','year','status'))
primary <- read.csv('data/sumstats_final.csv');primary$key <- uncertainty_key(primary)
annual <- uncertainty_annual(normalizePath('.'),primary);defs <- uncertainty_definitions(normalizePath('.'))
baseline <- d$contrasts[scheme=='baseline']
stopifnot(nrow(primary)==143,nrow(baseline)==143,identical(baseline$key,primary$key),
  max(abs(baseline$estimate-primary$Mean))<1e-7,max(abs(baseline$se-primary$SE))<1e-7)
published <- fread('tables/table_03_cross_taxa_meta_regression.csv')[
  Predictor=='P. interruptus biomass' & Response_Var=='S. purpuratus biomass']
path_base <- d$pathway[scheme=='baseline']
stopifnot(nrow(path_base)==1,abs(round(path_base$slope,4)-published$Slope_Est)<1e-8,
  abs(path_base$p-published$Slope_pval)<1e-8)
stochastic <- c('individual_weights','recorded_block_weights')
for(scheme_name in stochastic){
  e <- d$contrasts[scheme==scheme_name]
  stopifnot(uniqueN(e$draw)==status$draws_per_scheme,nrow(e)==143*status$draws_per_scheme)
  fixed <- e[!key %in% pools$series_id]
  i <- match(fixed$key,primary$key)
  stopifnot(all(fixed$estimate==primary$Mean[i]),all(fixed$se==primary$SE[i]))
  w <- d$weights[[scheme_name]]
  stopifnot(!anyDuplicated(w$source_row),all(w$weight>0))
  if(scheme_name=='recorded_block_weights')stopifnot(all(w[,uniqueN(weight),by=block]$V1==1))
  # Independently sum weighted individual masses in each actual consuming pool.
  joined <- merge(members,w,by='source_row',suffixes=c('_member','_source'))
  stopifnot(all(joined$count_member==joined$count_source))
  checked <- joined[,.(mass=sum(count_member*weight*mass)/sum(count_member*weight),
                       baseline=sum(count_member*mass)/sum(count_member)),by=pool]
  observed <- ratios[scheme==scheme_name & draw==1]
  checked <- merge(checked,observed,by='pool')
  stopifnot(nrow(checked)==nrow(pools),max(abs(checked$mass/checked$baseline-checked$ratio))<1e-12)
  singleton_ids <- pools[singleton==TRUE,pool]
  stopifnot(max(abs(ratios[scheme==scheme_name & pool %in% singleton_ids,ratio]-1))<1e-12)
  # Reconstruct changed log ratios and fit lm independently of production projection.
  domain <- fread(file.path(private,'full_annual_domain.csv'))[series_id %in% pools$series_id]
  z <- merge(pools[,.(pool,series_id,year,status)],observed[,.(pool,ratio)],by='pool')
  domain <- merge(domain,z,by=c('series_id','year','status'),all.x=TRUE)
  domain[is.na(ratio),ratio:=1];domain[,v:=value*ratio]
  domain[,corrected:=v+min(v[v>0])/2,by=series_id]
  for(id in unique(pools$series_id)){
    row <- primary[primary$key==id,];re <- uncertainty_reconstruct(row,annual,defs)
    a <- dcast(domain[series_id==id],year~status,value.var='corrected')
    nd <- re$data;ix <- match(nd$year,a$year);stopifnot(!anyNA(ix))
    nd$lnDiff <- log(a$mpa[ix]/a$reference[ix])
    fit <- if(row$Model=='Mean')lm(lnDiff~1,nd) else lm(lnDiff~time.model,nd)
    est <- sum(coef(fit)*re$gradient);se <- sqrt(drop(t(re$gradient)%*%vcov(fit)%*%re$gradient))
    want <- e[draw==1 & key==id]
    stopifnot(nrow(want)==1,abs(est-want$estimate)<1e-10,abs(se-want$se)<1e-10)
  }
}
# Full pre-pairing domain must retain the otherwise missing Gull Island minimum.
domain <- fread(file.path(private,'full_annual_domain.csv'))
gull <- domain[series_id=='Gull Island SMR|S. purpuratus|PISCO|Bio']
paired <- inventory[series_id=='Gull Island SMR|S. purpuratus|PISCO|Bio']
stopifnot(min(gull$value[gull$value>0]) < min(paired$value[paired$value>0]))
minimum <- gull[value==min(value[value>0])]
minimum_pool <- merge(minimum,pools,by=c('series_id','year','status'))
stopifnot(nrow(minimum_pool)==1,!minimum_pool$retained_paired_year,
  sd(ratios[scheme=='individual_weights' & pool==minimum_pool$pool,ratio])>0)
# Changing only the unpaired minimum changes the adjustment on retained years.
before_offset <- min(gull$value[gull$value>0])/2
gull[year==minimum$year & status==minimum$status,value:=value*2]
after_offset <- min(gull$value[gull$value>0])/2
stopifnot(after_offset>before_offset)
pair <- dcast(paired,year~status,value.var='value')[1]
stopifnot(abs(log((pair$mpa+before_offset)/(pair$reference+before_offset))-
  log((pair$mpa+after_offset)/(pair$reference+after_offset)))>1e-6)
# Hand calculation: partial overlap must use shared record weights, not a shared mean.
mass <- c(1,8,27);frequency <- c(2,1,1);w <- c(.5,2,1)
weighted <- function(ix)sum(mass[ix]*frequency[ix]*w[ix])/sum(frequency[ix]*w[ix])
stopifnot(abs(weighted(1:3)-11)<1e-12,abs(weighted(2:3)-43/3)<1e-12)
cat('PASS: fingerprints; fixed 143 histories; shared record/block weights; frequency counts; singleton limits; independent annual/LM reconstruction; unpaired-year offset domain.\n')
