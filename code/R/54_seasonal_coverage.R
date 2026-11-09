# Recorded survey months for the retained histories; no abundance estimates changed.
# DONHAM_DATA_REPO overrides the sibling source-data workspace.
local({
  suppressPackageStartupMessages({library(data.table);library(here)})
  root<-here::here();dr<-Sys.getenv('DONHAM_DATA_REPO',file.path(dirname(root),'sbc-donham-kelp-mpa-data'))
  if(!dir.exists(dr)) {
    manifest<-file.path(root,'tables/table_s_seasonal_inputs.csv')
    outputs<-file.path(root,'tables/table_s_seasonal_output_hashes.csv')
    if(!file.exists(manifest)||!file.exists(outputs))stop('Seasonal audit requires source data or verified derived exports')
    m<-fread(manifest)[repository=='analysis'];o<-fread(outputs)
    stopifnot(all(unname(tools::md5sum(file.path(root,m$path)))==m$md5),all(unname(tools::md5sum(file.path(root,o$path)))==o$md5))
    message('Seasonal audit: verified saved derived exports; raw inputs unavailable, so source coverage was not recomputed.')
  } else {
  used<-character();read<-function(p,repo=dr){f<-file.path(repo,p);used<<-c(used,f);fread(f)}
  primary<-read('data/sumstats_final.csv',root)
  keys<-read('data/harmonized/reference_primary_keys.csv',root)
  tax<-c('Macrocystis pyrifera'='M. pyrifera','Panulirus interruptus'='P. interruptus',
         'Semicossyphus pulcher'='S. pulcher','Strongylocentrotus purpuratus'='S. purpuratus',
         'Mesocentrotus franciscanus'='M. franciscanus')
  expr<-as.list(parse(file.path(dr,'code/R/01_utils.R')))
  constant<-function(n)eval(Filter(function(x)is.call(x)&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name(n)),expr)[[1]][[3]],baseenv())
  used<-c(used,file.path(dr,'code/R/01_utils.R'))
  pm<-unique(read('data/PISCO/master_site_table_Emilyedit.csv')[BaselineRegion=='SOUTH',.(site,MPA=CA_MPA_Name_Short,status=site_status)])
  pm<-pm[MPA %in% primary[Source=='PISCO',MPA]];stopifnot(!anyDuplicated(pm[,.(site,MPA)]))
  sw<-read('data/PISCO/MLPA_kelpforest_swath_2024.csv')
  sw<-unique(sw[(campus=='VRG'|(campus=='UCSB' & !is.na(zone)&!zone %in% c('DEEP','SPECIALO')))&month>=5&month<=10,.(site,year,month)])
  fi<-read('data/PISCO/MLPA_kelpforest_fish_2024.csv')
  fi<-unique(fi[campus %in% c('UCSB','VRG') & level!='CAN',.(site,year,month)])
  eligible<-unique(rbind(sw,fi)[,.(site,year)])[,.(n=.N),by=site][n>=5,site]
  survey<-list();add<-function(x){survey[[length(survey)+1]]<<-unique(x[,.(Source,Taxa,component,site,MPA,status,year,month)])}
  for(tt in c('M. pyrifera','P. interruptus','S. purpuratus','M. franciscanus','S. pulcher')){
    x<-copy(if(tt=='S. pulcher')fi else sw);x<-merge(x[site %in% eligible & !site %in% constant('EXCLUDED_REFERENCE_SITES')],pm,by='site',allow.cartesian=TRUE)
    if(tt=='S. pulcher') {
      policy<-read('tables/table_s_pisco_fixed_coverage.csv',root)[scenario=='strict_min3' & eligible==TRUE]
      allowed<-rbindlist(lapply(seq_len(nrow(policy)),function(i) {
        z<-policy[i];ss<-strsplit(paste(z$protected,z$reference,sep=';'),';',fixed=TRUE)[[1]]
        data.table(MPA=z$MPA,site=ss)
      }))
      x<-rbind(x[!MPA %in% policy$MPA],merge(x[MPA %in% policy$MPA],allowed,by=c('MPA','site')))
    }
    x[,`:=`(Source='PISCO',Taxa=tt,component='count_survey')];add(x)
  }
  mm<-unique(read('data/MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv')[,.(site=site_id,MPA=CA_MPA_Name_Short,status)])
  mb<-read('data/MBON/SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv')
  mb[,`:=`(year=as.integer(substr(date,1,4)),month=as.integer(substr(date,6,7)),site=site_id,Taxa=unname(tax[taxon_name]))]
  for(src in c('kfm','lter')){
    x<-mb[data_source==src & !is.na(Taxa)]
    if(src=='kfm')x<-x[!site %in% constant('EXCLUDED_KFM_SITES') & !(taxon_name=='Macrocystis pyrifera' & proj_taxon_id=='t-k-002')]
    if(src=='lter')x<-x[!site %in% constant('EXCLUDED_LTER_SITES') & Taxa %in% c('S. purpuratus','M. franciscanus') & sample_method=='quad']
    x<-merge(x,mm,by='site',allow.cartesian=TRUE);x[,`:=`(Source=toupper(src),component='count_survey')];add(x)
  }
  kf<-read('data/MBON/SBCMBON_kelp_forest_integrated_fish_20231022.csv')
  kf<-kf[data_source=='kfm' & proj_taxon_id %in% c('t-k-127','t-k-233') & sample_method=='visualfish']
  kf[,`:=`(site=site_id,year=as.integer(substr(date,1,4)),month=as.integer(substr(date,6,7)))]
  kf<-merge(kf,mm,by='site',allow.cartesian=TRUE);kf[,`:=`(Source='KFM',Taxa='S. pulcher',component='count_survey')];add(kf)
  for(z in list(c('Lobster_Abundance_All_Years_20230831.csv','P. interruptus'),c('Annual_Kelp_All_Years_20240305.csv','M. pyrifera'),c('Annual_fish_comb_20240307.csv','S. pulcher'))){
    x<-read(paste0('data/LTER/',z[1]));if(z[2]=='M. pyrifera')x<-x[FRONDS!= -99999 & !SITE %in% constant('EXCLUDED_LTER_SITE_CODES')]
    if(z[2]=='S. pulcher')x<-x[SCIENTIFIC_NAME=='Semicossyphus pulcher' & COUNT!= -99999 & !SITE %in% constant('EXCLUDED_LTER_SITE_CODES')]
    x[,`:=`(site=SITE,year=YEAR,month=MONTH)];x<-merge(x,mm,by='site',allow.cartesian=TRUE)
    x[,`:=`(Source='LTER',Taxa=z[2],component='count_survey')];add(x)
  }
  st<-read('data/MBON/KFM_Macrocystis_RawData_1984-2023.csv');st[,`:=`(site=sprintf('a-k-%02d',SiteNumber),date=as.Date(SurveyDate,'%m/%d/%Y'))]
  st<-st[!site %in% constant('EXCLUDED_KFM_SITES') & ScientificName=='Macrocystis pyrifera']
  st[,`:=`(year=as.integer(format(date,'%Y')),month=as.integer(format(date,'%m')))]
  st<-merge(st,mm,by='site',allow.cartesian=TRUE);st[,`:=`(Source='KFM',Taxa='M. pyrifera',component='kelp_size_sample')];add(st)
  # Separate size pools are matched by MPA, status, year and species, not survey month.
  sf<-read('data/ALL_sizefreq_2024.csv')
  sf<-sf[(campus %in% c('UCSB','VRG') & classcode %in% c('MESFRA','STRPUR','PANINT')) | (campus=='KFM' & classcode %in% c('Strongylocentrotus franciscanus','Strongylocentrotus purpuratus')) | (campus=='LTER' & classcode %in% c('SFL','SFS','SPS','SPL'))]
  sf<-sf[!site %in% c(constant('EXCLUDED_KFM_SITES'),constant('EXCLUDED_REFERENCE_SITES'))]
  sm<-unique(read('data/PISCO/master_site_table_Emilyedit.csv')[,.(site,site_new,campus,MPA=CA_MPA_Name_Short,status=site_status)])
  sf<-merge(sf,sm,by=c('site','site_new','campus'),all.x=TRUE,allow.cartesian=TRUE)
  ix<-match(sf$site,mm$site);missing_mpa<-which(is.na(sf$MPA));sf[missing_mpa,MPA:=mm$MPA[ix[missing_mpa]]]
  missing_status<-which(is.na(sf$status));sf[missing_status,status:=mm$status[ix[missing_status]]]
  sf[,Taxa:=fcase(classcode %in% c('MESFRA','Strongylocentrotus franciscanus','SFL','SFS'),'M. franciscanus',classcode %in% c('STRPUR','Strongylocentrotus purpuratus','SPS','SPL'),'S. purpuratus',classcode=='PANINT','P. interruptus',default=NA_character_)]
  sf<-sf[!is.na(Taxa) & status %in% c('mpa','reference')]
  # Size dates are separate metadata; preserve source-program identity and unknown months.
  size_events<-rbindlist(lapply(c('PISCO','KFM','LTER'),function(src){x<-sf[Taxa!='P. interruptus' | (src=='PISCO' & campus %in% c('UCSB','VRG'))];if(src=='PISCO')x<-x[Taxa=='P. interruptus'|size>=25];unique(x[,.(Source=src,Taxa,MPA,status,year,month,size_source=campus,site)])}))
  raw<-unique(rbindlist(survey));raw<-raw[status %in% c('mpa','reference')]
  stopifnot(all(raw$month %in% 1:12),!anyNA(raw$year))
  # Join on every response's actual paired calendar years. Month coverage does not
  # imply that all raw biomass values passed size/measurement checks.
  membership<-rbindlist(lapply(seq_len(nrow(keys[Source!='Landsat'])),function(i){k<-keys[Source!='Landsat'][i];data.table(MPA=k$MPA,Taxa=k$Taxa,Source=k$Source,Resp=k$Resp,year=as.integer(strsplit(k$observed_years,',',fixed=TRUE)[[1]]))}))
  events<-merge(raw,membership,by=c('MPA','Taxa','Source','year'),allow.cartesian=TRUE)
  expected<-merge(membership,primary[,.(MPA,Taxa,Source,Resp,MPA_Start)],by=c('MPA','Taxa','Source','Resp'))
  size_events<-merge(size_events,membership[Resp=='Bio'],by=c('MPA','Taxa','Source','year'),allow.cartesian=TRUE)
  pairs<-rbindlist(lapply(seq_len(nrow(expected)),function(i){k<-expected[i];d<-events[MPA==k$MPA & Taxa==k$Taxa & Source==k$Source & Resp==k$Resp & year==k$year & component=='count_survey'];a<-sort(unique(d[status=='mpa',month]));b<-sort(unique(d[status=='reference',month]));common<-intersect(a,b)
    data.table(k,protected_months=paste(a,collapse=';'),reference_months=paste(b,collapse=';'),n_common_months=length(common),same_month_set=identical(a,b),date_coverage_complete=length(a)>0 & length(b)>0,
      minimum_month_gap=if(length(a)&&length(b))min(abs(outer(a,b,'-')))else NA_real_,
      protected_mean_month=if(length(a))mean(d[status=='mpa',month])else NA_real_,reference_mean_month=if(length(b))mean(d[status=='reference',month])else NA_real_)
  }))
  pairs[,period:=ifelse(year<=MPA_Start,'Before','After')]
  summary<-pairs[,.(paired_years=.N,dates_both_sides=sum(date_coverage_complete),same_month_set=sum(same_month_set & date_coverage_complete),shared_month=sum(n_common_months>0),no_shared_month=sum(n_common_months==0 & date_coverage_complete),gap_at_least_two_months=sum(minimum_month_gap>=2,na.rm=TRUE)),by=.(Source,Taxa,Resp)]
  months<-events[,.(site_months=.N),by=.(Source,Taxa,Resp,component,status,month)]
  # Trend in timing: site-month weighted centroid, per year, not species-record counts.
  shifts<-pairs[,.(paired_years=.N,mean_month_gap=mean(protected_mean_month-reference_mean_month,na.rm=TRUE),no_shared_month=sum(n_common_months==0)),by=.(MPA,Source,Taxa,Resp,period)]
  write<-function(x,n)fwrite(x,file.path(root,'tables',paste0('table_s_seasonal_',n,'.csv')))
  write(size_events,'size_months');
  write(events,'site_months');write(pairs,'paired_years');write(summary,'summary');write(months,'months');write(shifts,'periods')
  used<-unique(c(used,file.path(root,'code/R/54_seasonal_coverage.R'),file.path(dr,paste0('code/R/',c('03_data_import.R','04_pisco_processing.R','05_kfm_processing.R','06_lter_processing.R')))))
  write(data.table(repository=ifelse(startsWith(used,paste0(dr,'/')),'data','analysis'),path=ifelse(startsWith(used,paste0(dr,'/')),substring(used,nchar(dr)+2),substring(used,nchar(root)+2)),md5=unname(tools::md5sum(used))),'inputs')
  output_paths<-list.files(file.path(root,'tables'),pattern='^table_s_seasonal_.*[.]csv$',full.names=TRUE)
  output_paths<-output_paths[basename(output_paths)!='table_s_seasonal_output_hashes.csv']
  write(data.table(path=file.path('tables',basename(output_paths)),md5=unname(tools::md5sum(output_paths))),'output_hashes')
  print(summary);cat('Coverage only: separate size dates are audited where recorded; annual biomass validity and seasonally matched responses are not reconstructed.\n')
  }
})
