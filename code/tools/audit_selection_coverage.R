# Reconstruct survey-site continuity without changing abundance or model inputs.
# Rscript code/tools/audit_selection_coverage.R /path/to/data-processing-repo
local({
  suppressPackageStartupMessages({library(data.table);library(here)})
  args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1L)
  root<-here::here();dr<-normalizePath(args[1],mustWork=TRUE)
  paths<-file.path(dr,c('data/PISCO/MLPA_kelpforest_swath_2024.csv',
    'data/PISCO/MLPA_kelpforest_fish_2024.csv','data/PISCO/master_site_table_Emilyedit.csv','code/R/01_utils.R'))
  paths<-c(paths,file.path(root,c('data/sumstats_final.csv','data/harmonized/reference_primary_keys.csv',
    'code/tools/audit_selection_coverage.R')))
  hashes<-tools::md5sum(paths)
  expr<-Filter(function(x)is.call(x)&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name('EXCLUDED_REFERENCE_SITES')),as.list(parse(paths[4])))
  stopifnot(length(expr)==1);excluded<-eval(expr[[1]][[3]],baseenv())
  sw<-fread(paths[1],select=c('site','year','month','campus','zone'))
  sw<-unique(sw[(campus=='VRG'|(campus=='UCSB' & !is.na(zone) & !zone %in% c('DEEP','SPECIALO'))) & month>=5 & month<=10,.(site,year)])
  sw[,branch:='benthic']
  fish_records<-fread(paths[2],select=c('site','year','campus','level','classcode','count','fish_tl'))
  fish_records<-fish_records[campus %in% c('UCSB','VRG') & level!='CAN']
  bad_fish<-unique(fish_records[classcode=='SPUL' & (!is.finite(count)|count<0|(count>0 & (!is.finite(fish_tl)|fish_tl<=0))),.(site,year)])
  fwrite(bad_fish,file.path(root,'tables/table_s_selection_fish_invalid_size_years.csv'))
  fi<-unique(fish_records[,.(site,year)]);fi[,branch:='fish']
  obs<-rbind(sw,fi)
  # A surveyed event is evidence of coverage, not proof of a valid biomass value.
  eligible_sites<-obs[,.(n=uniqueN(year)),by=site][n>=5,site]
  primary<-fread(file.path(root,'data/sumstats_final.csv'))
  map<-unique(fread(paths[3])[BaselineRegion=='SOUTH',.(site,MPA=CA_MPA_Name_Short,status=site_status)])
  conflicts<-map[map[, .N,by=.(site,MPA)][N>1],on=.(site,MPA)]
  fwrite(conflicts,file.path(root,'tables/table_s_selection_mapping_conflicts.csv'))
  map<-map[MPA %in% primary[Source=='PISCO',MPA]]
  stopifnot(!anyDuplicated(map[,.(site,MPA)]))
  obs<-merge(obs[site %in% eligible_sites & !site %in% excluded],map,by='site',allow.cartesian=TRUE)
  keys<-fread(file.path(root,'data/harmonized/reference_primary_keys.csv'))
  starts<-unique(primary[Source=='PISCO',.(MPA,MPA_Start)]);stopifnot(!anyDuplicated(starts$MPA))
  obs<-merge(obs,starts,by='MPA');obs<-obs[status %in% c('mpa','reference')]
  out<-list();membership<-list()
  years<-function(s)as.integer(strsplit(s,',',fixed=TRUE)[[1]])
  for(m in starts$MPA)for(b in c('fish','benthic')) {
    k<-keys[MPA==m & Source=='PISCO' & (if(b=='fish')Taxa=='S. pulcher' else Taxa!='S. pulcher')]
    if(!nrow(k))next
    common<-Reduce(intersect,lapply(k$observed_years,years));o<-obs[MPA==m & branch==b]
    for(s in unique(o$site)) {
      z<-o[site==s];membership[[length(membership)+1]]<-data.frame(MPA=m,branch=b,site=s,status=z$status[1],
        years=paste(sort(unique(z$year)),collapse=';'))
    }
    for(a in unique(o[status=='mpa',site]))for(c in unique(o[status=='reference',site])) {
      paired<-intersect(o[site==a,year],o[site==c,year]);start<-o$MPA_Start[1]
      for(scope in c('survey_coverage','all_retained_branch_responses',if(b=='fish')'fish_count_length_checked')) {
        yy<-if(scope=='survey_coverage')paired else intersect(paired,common)
        if(scope=='fish_count_length_checked')yy<-setdiff(yy,bad_fish[site %in% c(a,c),year])
        out[[length(out)+1]]<-data.frame(MPA=m,branch=b,scope,protected_site=a,reference_site=c,
          designation_year=start,n_pre=sum(yy<start),n_through_designation=sum(yy<=start),n_after=sum(yy>start),
          strict_before_years=paste(sort(yy[yy<start]),collapse=';'),
          through_designation_years=paste(sort(yy[yy<=start]),collapse=';'),
          all_paired_years=paste(sort(yy),collapse=';'))
      }
    }
  }
  out<-rbindlist(out);membership<-rbindlist(membership)
  summary<-out[,.(max_strict_before=max(n_pre),max_through_designation=max(n_through_designation),
    pairs_three_strict=sum(n_pre>=3 & n_after>=3),pairs_four_strict=sum(n_pre>=4 & n_after>=4),
    pairs_three_through=sum(n_through_designation>=3 & n_after>=3),pairs_four_through=sum(n_through_designation>=4 & n_after>=4)),by=.(MPA,branch,scope)]
  fwrite(out,file.path(root,'tables/table_s_selection_station_pairs.csv'))
  fwrite(summary,file.path(root,'tables/table_s_selection_branch_coverage.csv'))
  fwrite(membership,file.path(root,'tables/table_s_selection_station_years.csv'))
  stopifnot(identical(hashes,tools::md5sum(paths)))
  fwrite(data.frame(repository=ifelse(startsWith(paths,paste0(dr,'/')),'data','analysis'),
    input=ifelse(startsWith(paths,paste0(dr,'/')),substring(paths,nchar(dr)+2),substring(paths,nchar(root)+2)),
    md5=unname(hashes)),file.path(root,'tables/table_s_selection_source_inputs.csv'))
  print(summary[scope=='fish_count_length_checked' & (pairs_three_through>0|pairs_four_through>0)])
  cat('Initial grids are survey-coverage upper bounds; the fish-only check removes invalid-count/missing-positive-length site-years without refitting biomass.\n')
})
