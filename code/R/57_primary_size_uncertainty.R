# Conditional size-sample sensitivity, NOT a replacement sampling confidence interval.
# Shared original observations receive the same weights in every consuming history.
# Count surveys, membership, allometry and the two directly measured LTER series stay fixed.
local({
  suppressPackageStartupMessages({library(data.table); library(Matrix); library(metafor)})
  source('code/R/primary_uncertainty_utils.R', local=TRUE)
  root <- normalizePath('.')
  dr <- Sys.getenv('DONHAM_DATA_REPO', file.path(dirname(root), 'sbc-donham-kelp-mpa-data'))
  B <- as.integer(Sys.getenv('DONHAM_SIZE_DRAWS', '1000'))
  stopifnot(B >= 20)
  seed <- as.integer(Sys.getenv('DONHAM_SIZE_SEED','20261030'));set.seed(seed)
  output_dir <- Sys.getenv('DONHAM_SIZE_OUTPUT_DIR','tables')
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  protected <- c('data/sumstats_final.csv', 'data/harmonized/harmonized_raw_responses.csv',
    'data/harmonized/harmonized_response_ratios.csv', 'tables/table_02_meta_analysis.csv',
    'tables/table_03_cross_taxa_meta_regression.csv', 'data/cache/temporal_primary_models.rds')
  before <- tools::md5sum(protected)
  inventory_inputs <- fread('tables/table_s_primary_size_inputs.csv')
  inventory_paths <- ifelse(inventory_inputs$repository=='analysis',
    file.path(root,inventory_inputs$path),file.path(dr,inventory_inputs$path))
  stopifnot(all(unname(tools::md5sum(inventory_paths)) == inventory_inputs$md5))
  inventory <- fread('tables/table_s_primary_size_inventory.csv')
  pools <- inventory[method == 'resampled size pool' & positive == TRUE]
  pools[, retained_paired_year := TRUE]
  ids <- unique(pools$series_id)
  domain_path <- 'local_archive/primary_size_uncertainty_2026-10-30/full_annual_domain.csv'
  domain_manifest <- sub('full_annual_domain.csv','domain_inputs.csv',domain_path)
  domain_inputs <- fread(domain_manifest)
  stopifnot(all(unname(tools::md5sum(domain_inputs$path))==domain_inputs$md5))
  domain <- fread(domain_path)[series_id %in% ids]
  membership_path <- 'local_archive/primary_size_inventory_2026-10-30/record_membership.csv'
  membership <- fread(membership_path)[source_file == 'ALL_sizefreq_2024.csv']
  size_path <- file.path(dr, 'data/ALL_sizefreq_2024.csv')
  raw_sizes <- fread(size_path); raw_sizes[, source_row := .I]
  extra <- domain[value>0][!pools,on=.(series_id,year,status)]
  # Recover sizes for unpaired offset-domain observations too. Use the same
  # original-row mapping, source filters and diameter cutoffs as module 55.
  expr <- as.list(parse(file.path(dr,'code/R/01_utils.R')))
  constant <- function(n)eval(Filter(function(x)is.call(x)&&identical(x[[1]],as.name('<-'))&&
    identical(x[[2]],as.name(n)),expr)[[1]][[3]],baseenv())
  sf <- copy(raw_sizes[(campus %in% c('VRG','UCSB') & classcode %in% c('STRPUR','PANINT')) |
    (campus=='KFM' & classcode=='Strongylocentrotus purpuratus') |
    (campus=='LTER' & classcode %in% c('SPS','SPL'))])
  sf <- sf[!site %in% c(constant('EXCLUDED_REFERENCE_SITES'),constant('EXCLUDED_KFM_SITES'))]
  sm <- fread(file.path(dr,'data/PISCO/master_site_table_Emilyedit.csv'))
  sm <- sm[!duplicated(sm[,.(site,CA_MPA_Name_Short)]),
    .(site,site_new,campus,MPA=CA_MPA_Name_Short,status=site_status)]
  sf <- merge(sf,sm,by=c('site','site_new','campus'),all.x=TRUE,allow.cartesian=TRUE)
  mm <- unique(fread(file.path(dr,'data/MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv'))[
    ,.(site=site_id,MPA=CA_MPA_Name_Short,status)])
  ix <- match(sf$site,mm$site)
  i <- which(is.na(sf$MPA)&sf$classcode!='PANINT');sf[i,MPA:=mm$MPA[ix[i]]]
  i <- which(is.na(sf$status)&sf$classcode!='PANINT');sf[i,status:=mm$status[ix[i]]]
  sf[,Taxa:=ifelse(classcode=='PANINT','P. interruptus','S. purpuratus')]
  sf <- sf[!is.na(size)&!is.na(count)&count>0&status %in% c('mpa','reference')]
  stopifnot(all(sf$size>0),all(sf$count==floor(sf$count)))
  pool_rows <- rbindlist(lapply(c('PISCO','KFM'),function(src){
    z <- copy(sf[if(src=='PISCO') Taxa=='P. interruptus'|size>=25 else Taxa=='S. purpuratus'])
    z[,Source:=src];z
  }))
  extra_members <- merge(extra,pool_rows,by=c('MPA','Taxa','Source','year','status'),allow.cartesian=TRUE)
  extra_members[,mass:=ifelse(Taxa=='P. interruptus',.001352821*size^2.913963113,.00059*size^2.870)]
  extra_stats <- extra_members[,.(n_sized=sum(count),mean_mass=sum(count*mass)/sum(count)),
    by=.(series_id,year,status)]
  stopifnot(nrow(extra_stats)==nrow(extra))
  extra <- merge(extra,extra_stats,by=c('series_id','year','status'))
  extra[,`:=`(retained_paired_year=FALSE,positive=TRUE,method='resampled size pool',
    sparse_lt10=n_sized<10,singleton=n_sized==1)]
  pools <- rbindlist(list(pools,extra),fill=TRUE);pools[,pool:=.I]
  membership <- rbindlist(list(membership,extra_members[,.(series_id,year,status,
    source_file='ALL_sizefreq_2024.csv',source_row,count,site)]),use.names=TRUE)
  membership_output <- sub('full_annual_domain.csv','complete_record_membership.csv',domain_path)
  fwrite(membership,membership_output)
  sizes <- raw_sizes[source_row %in% membership$source_row]
  sizes[, mass := fifelse(classcode == 'PANINT', .001352821 * size^2.913963113,
                          .00059 * size^2.870)]
  # Missing survey identifiers are kept as a single recorded block, never invented.
  block_columns <- c('campus','method','year','site','month','day','location','zone','transect')
  sizes[, block := do.call(paste, c(.SD, sep='|')), .SDcols=block_columns]
  sizes[, block_index := match(block, unique(block))]
  members <- merge(membership, pools[,.(series_id,year,status,pool)],
                   by=c('series_id','year','status'))
  stopifnot(nrow(members)==nrow(membership), !anyDuplicated(members[,.(pool,source_row)]))
  j <- match(members$source_row, sizes$source_row)
  shared_rows <- members[,.(n_consumers=uniqueN(series_id)),by=source_row][n_consumers>1,source_row]
  pools[,shared_with_other_history:=pool %in% members[source_row %in% shared_rows,pool]]
  stopifnot(!anyNA(j), all(members$count==sizes$count[j]), all(sizes$count>0),
            all(sizes$count==floor(sizes$count)), all(is.finite(sizes$mass)))
  incidence <- sparseMatrix(i=members$pool, j=j, x=members$count,
                            dims=c(nrow(pools),nrow(sizes)))
  mean_mass <- function(weights) as.numeric(incidence %*% (sizes$mass*weights)) /
    as.numeric(incidence %*% weights)
  baseline_mass <- mean_mass(rep(1,nrow(sizes)))
  stopifnot(max(abs(baseline_mass/pools$mean_mass-1))<1e-12,
            identical(as.numeric(rowSums(incidence)), as.numeric(pools$n_sized)))
  blocks <- data.table(pool=members$pool,block=sizes$block[j])[,.(n_recorded_blocks=uniqueN(block)),by=pool]
  pools <- merge(pools,blocks,by='pool',sort=FALSE);setorder(pools,pool)
  pools[, lacks_finer_survey_ids := vapply(pool,function(p)
    any(is.na(sizes$transect[j[members$pool==p]])), logical(1))]

  primary <- read.csv(protected[1]);primary$key <- uncertainty_key(primary)
  annual <- uncertainty_annual(root,primary);defs <- uncertainty_definitions(root)
  selected <- primary[match(ids,primary$key),]
  stopifnot(length(ids)==13, all(selected$Model %in% c('Mean','Linear')))
  reconstructed <- setNames(lapply(seq_along(ids),function(i)
    uncertainty_reconstruct(selected[i,],annual,defs)), ids)
  all_raw <- fread(protected[2])[resp=='Bio']
  all_raw[, Taxa := fcase(taxon_name=='Panulirus interruptus','P. interruptus',
    taxon_name=='Strongylocentrotus purpuratus','S. purpuratus',default=NA_character_)]
  all_raw[, key := paste(CA_MPA_Name_Short,Taxa,source,resp,sep='|')]
  all_raw <- all_raw[key %in% ids]
  setnames(domain,'series_id','key')
  check_domain <- merge(all_raw,domain,by=c('key','year','status'),suffixes=c('_released','_domain'))
  stopifnot(nrow(check_domain)==nrow(all_raw),
    max(abs(check_domain$value_released-check_domain$value_domain))<1e-8)
  all_raw <- domain
  raw_pool <- match(paste(all_raw$key,all_raw$year,all_raw$status),
                    paste(pools$series_id,pools$year,pools$status))
  # Recalculate the adjustment over its full original domain, including unpaired
  # years. Their size variation affects retained contrasts through that adjustment.
  structure <- setNames(lapply(ids,function(id){
    r <- reconstructed[[id]]; ii <- which(all_raw$key==id)
    a <- all_raw[ii]; idx_m <- match(paste(r$data$year,'mpa'),paste(a$year,a$status))
    idx_r <- match(paste(r$data$year,'reference'),paste(a$year,a$status))
    stopifnot(!anyNA(idx_m),!anyNA(idx_r))
    X <- model.matrix(r$fit); inv <- solve(crossprod(X))
    list(rows=ii,mpa=idx_m,ref=idx_r,X=X,inv=inv,projection=inv%*%t(X),
         g=r$gradient,n=nrow(X),p=ncol(X),primary_row=match(id,primary$key))
  }),ids)
  effects <- function(ratios){
    stopifnot(length(ratios)==nrow(pools), all(is.finite(ratios)),all(ratios>0))
    v <- all_raw$value;use <- which(!is.na(raw_pool));v[use] <- v[use]*ratios[raw_pool[use]]
    x <- primary
    for(id in ids){
      s <- structure[[id]];z <- v[s$rows];z <- z+min(z[z>0])/2
      y <- log(z[s$mpa]/z[s$ref]);beta <- s$projection%*%y
      sigma2 <- sum((y-s$X%*%beta)^2)/(s$n-s$p)
      x$Mean[s$primary_row] <- sum(s$g*beta)
      x$SE[s$primary_row] <- sqrt(as.numeric(t(s$g)%*%s$inv%*%s$g)*sigma2)
    }
    x
  }
  baseline <- effects(rep(1,nrow(pools)))
  if(max(abs(baseline$Mean-primary$Mean))>=1e-7) print(as.data.table(data.frame(key=primary$key,
    original=primary$Mean,reconstructed=baseline$Mean))[abs(original-reconstructed)>1e-7])
  stopifnot(max(abs(baseline$Mean-primary$Mean))<1e-7,
            max(abs(baseline$SE-primary$SE))<1e-7)
  baseline_joint <- uncertainty_joint(primary,primary$SE,'baseline')
  published <- fread(protected[4]);published[, Resp:=ifelse(Response=='Biomass','Bio','Den')]
  check <- merge(as.data.table(baseline_joint),published,by=c('Taxa','Resp'))
  stopifnot(nrow(check)==9,max(abs(check$estimate-check$Estimate))<1e-7,
            max(abs(check$q-check$pval_fdr))<1e-7)
  joint <- function(x,label){
    d <- x[x$Resp=='Bio',];d$Taxa <- factor(d$Taxa);d$Source <- factor(d$Source)
    m <- rma.mv(Mean,V=SE^2,mods=~Taxa-1,random=list(~1|MPA,~1|Source),
                method='REML',test='t',data=d)
    z <- data.frame(scenario=label,Taxa=sub('^Taxa','',names(coef(m))),Resp='Bio',
      estimate=as.numeric(coef(m)),se=m$se,ci_low=m$ci.lb,ci_high=m$ci.ub,p=m$pval)
    den <- baseline_joint[baseline_joint$Resp=='Den',names(z)];den$scenario <- label
    z <- rbind(z,den);z$q <- p.adjust(z$p,'BH');z
  }
  pathway <- function(x,label){
    z <- as.data.table(x)[Resp=='Bio' & Taxa %in% c('P. interruptus','S. purpuratus'),
      .(Mean=mean(Mean),vi=sum(SE^2)/.N^2),by=.(MPA,Taxa)]
    d <- merge(z[Taxa=='P. interruptus'],z[Taxa=='S. purpuratus'],by='MPA',suffixes=c('_x','_y'))
    stopifnot(nrow(d)==4)
    m <- rma(Mean_y,vi_y,mods=~Mean_x,data=d,method='REML',test='knha')
    data.table(scenario=label,k=4,slope=unname(coef(m)[2]),ci_low=m$ci.lb[2],
               ci_high=m$ci.ub[2],p=m$pval[2])
  }
  estimates <- list();pooled <- list();paths <- list();draw_ratios <- list();idx <- 0L
  validation_weights <- list()
  record <- function(ratios,scheme,draw,detail=''){
    x <- effects(ratios);label <- paste(scheme,draw,sep=':');idx <<- idx+1L
    estimates[[idx]] <<- as.data.table(data.frame(scheme,draw,detail,key=x$key,estimate=x$Mean,se=x$SE))
    z <- as.data.table(joint(x,label));z[,`:=`(scheme=scheme,draw=draw,detail=detail)]
    pooled[[idx]] <<- z
    z <- pathway(x,label);z[,`:=`(scheme=scheme,draw=draw,detail=detail)];paths[[idx]] <<- z
    draw_ratios[[idx]] <<- data.table(scheme,draw,pool=pools$pool,ratio=ratios)
  }
  record(rep(1,nrow(pools)),'baseline',0)
  for(scheme in c('individual_weights','recorded_block_weights')){
    for(b in seq_len(B)){
      # Gamma(count)/count is the shared multiplier for a frequency-bin record:
      # equivalent to summing independent Exp(1) weights for its counted animals.
      w <- if(scheme=='individual_weights') rgamma(nrow(sizes),shape=sizes$count)/sizes$count else
        rexp(max(sizes$block_index))[sizes$block_index]
      if(b==1) validation_weights[[scheme]] <- data.table(source_row=sizes$source_row,
        count=sizes$count,mass=sizes$mass,block=sizes$block,weight=w)
      record(mean_mass(w)/baseline_mass,scheme,b)
      if(b%%100==0) message(scheme,': ',b,'/',B)
    }
  }
  # Deliberately imposed mass stress, not estimated variance or imputed lengths.
  for(scope in c('singleton','sparse_lt10')){
    hit <- which(pools$Taxa=='P. interruptus' & pools[[scope]] & pools$retained_paired_year)
    for(direction in c(-1,1)){
      ratios <- rep(1,nrow(pools));ratios[hit] <- 2^(direction*ifelse(pools$status[hit]=='mpa',1,-1))
      record(ratios,paste0(scope,'_joint_stress'),direction,'protected/reference masses changed in opposite directions')
    }
  }
  hit <- which(pools$Taxa=='P. interruptus' & pools$sparse_lt10 & pools$retained_paired_year)
  for(p in hit)for(multiplier in c(.5,2)){
    ratios <- rep(1,nrow(pools));ratios[p] <- multiplier
    record(ratios,'sparse_one_pool_stress',p*10+as.integer(multiplier*2),
      paste(pools$series_id[p],pools$year[p],pools$status[p],multiplier,sep=';'))
  }
  out <- function(d,stem)fwrite(d,file.path(output_dir,paste0('table_s_size_uncertainty_',stem,'.csv')))
  e <- rbindlist(estimates);p <- rbindlist(pooled);a <- rbindlist(paths)
  # Keep compact results public; draw-level and source-row diagnostics remain local.
  stochastic <- c('individual_weights','recorded_block_weights')
  target_ids <- primary$key[primary$Resp=='Bio' & primary$Taxa %in% c('P. interruptus','S. purpuratus')]
  out(e[!scheme %in% stochastic & key %in% target_ids],'contrasts')
  out(p[!scheme %in% stochastic],'pooled');out(a[!scheme %in% stochastic],'pathway')
  out(pools[,.(pool,series_id,year,status,n_sized,n_recorded_blocks,lacks_finer_survey_ids,
    singleton,sparse_lt10,shared_with_other_history,retained_paired_year)],'pools')
  quant <- function(v)as.list(setNames(as.numeric(quantile(v,c(.025,.5,.975))),c('q025','median','q975')))
  out(e[scheme %in% c('individual_weights','recorded_block_weights') & key %in% ids,
    c(quant(estimate),list(size_sd=sd(estimate),n_draws=.N)),by=.(scheme,key)],'contrast_summary')
  out(p[scheme %in% c('individual_weights','recorded_block_weights'),
    c(quant(estimate),list(fraction_q_lt05=mean(q<.05),fraction_positive=mean(estimate>0),
      n_draws=.N)),by=.(scheme,Taxa,Resp)],'pooled_summary')
  out(a[scheme %in% c('individual_weights','recorded_block_weights'),
    c(quant(slope),list(fraction_negative=mean(slope<0),fraction_nominal_ci_below_zero=mean(ci_high<0),
      n_draws=.N)),by=scheme],'pathway_summary')
  # Shared-sample covariance of conditional point-estimate perturbations.
  covs <- lapply(c('individual_weights','recorded_block_weights'),function(s){
    wide <- dcast(e[scheme==s & key %in% ids],draw~key,value.var='estimate')
    C <- cov(as.matrix(wide[,-1]));z <- as.data.table(as.table(C));setnames(z,c('key1','key2','covariance'))
    z[,scheme:=s];z
  });out(rbindlist(covs),'covariance')
  private <- 'local_archive/primary_size_uncertainty_2026-10-30'
  dir.create(private,recursive=TRUE,showWarnings=FALSE)
  saveRDS(list(contrasts=e,pooled=p,pathway=a,weights=validation_weights),
          file.path(private,paste0('draws_seed',seed,'.rds')))
  saveRDS(rbindlist(draw_ratios),file.path(private,paste0('pool_ratios_seed',seed,'.rds')))
  stopifnot(identical(before,tools::md5sum(protected)))
  inputs <- c(protected,size_path,membership_path,membership_output,domain_path,domain_manifest,'tables/table_s_primary_size_inventory.csv',
    'tables/table_s_primary_size_inputs.csv','code/R/57_primary_size_uncertainty.R','code/R/primary_uncertainty_utils.R')
  out(data.table(path=inputs,md5=unname(tools::md5sum(inputs))),'inputs')
  out(data.table(draws_per_scheme=B,seed=seed,affected_histories=length(ids),
    primary_histories=nrow(primary),positive_pools=nrow(pools),primary_unchanged=TRUE), 'status')
  message('Conditional size sensitivity completed; primary files unchanged.')
})
