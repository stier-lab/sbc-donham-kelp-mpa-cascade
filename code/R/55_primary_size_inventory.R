# Inventory only: no refits, new pooling rules, seasonal restrictions or primary edits.
local({
 suppressPackageStartupMessages(library(data.table))
 root<-normalizePath('.');dr<-Sys.getenv('DONHAM_DATA_REPO',file.path(dirname(root),'sbc-donham-kelp-mpa-data'))
 used<-character();read<-function(p,repo=root){f<-file.path(repo,p);used<<-c(used,f);fread(f)}
 primary_files<-file.path(root,c('data/sumstats_final.csv','data/harmonized/harmonized_raw_responses.csv','tables/table_02_meta_analysis.csv','tables/table_03_cross_taxa_meta_regression.csv'));before<-tools::md5sum(primary_files)
 keys<-read('data/harmonized/reference_primary_keys.csv')[Resp=='Bio' & Taxa %in% c('P. interruptus','S. purpuratus')]
 years<-keys[,.(year=as.integer(strsplit(observed_years,',',fixed=TRUE)[[1]])),by=.(series_id,MPA,Taxa,Source)]
 raw<-read('data/harmonized/harmonized_raw_responses.csv')[resp=='Bio']
 raw[,Taxa:=fcase(taxon_name=='Panulirus interruptus','P. interruptus',taxon_name=='Strongylocentrotus purpuratus','S. purpuratus',default=NA_character_)]
 setnames(raw,c('CA_MPA_Name_Short','source'),c('MPA','Source'))
 annual<-merge(years,raw[,.(MPA,Taxa,Source,year,status,value)],by=c('MPA','Taxa','Source','year'))
 stopifnot(nrow(annual)==2*nrow(years),!anyDuplicated(annual[,.(series_id,year,status)]),all(is.finite(annual$value)))
 expr<-as.list(parse(file.path(dr,'code/R/01_utils.R')))
 constant<-function(n)eval(Filter(function(x)is.call(x)&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name(n)),expr)[[1]][[3]],baseenv())
 sf<-read('data/ALL_sizefreq_2024.csv',dr);sf[,source_row:=.I]
 sf<-sf[(campus %in% c('VRG','UCSB') & classcode %in% c('STRPUR','PANINT'))|(campus=='KFM' & classcode=='Strongylocentrotus purpuratus')|(campus=='LTER' & classcode %in% c('SPS','SPL'))]
 sf<-sf[!site %in% c(constant('EXCLUDED_REFERENCE_SITES'),constant('EXCLUDED_KFM_SITES'))]
 sm<-read('data/PISCO/master_site_table_Emilyedit.csv',dr)
 sm<-sm[!duplicated(sm[,.(site,CA_MPA_Name_Short)]),.(site,site_new,campus,MPA=CA_MPA_Name_Short,status=site_status)]
 sf<-merge(sf,sm,by=c('site','site_new','campus'),all.x=TRUE,allow.cartesian=TRUE)
 mm<-unique(read('data/MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv',dr)[,.(site=site_id,MPA=CA_MPA_Name_Short,status)])
 # Backfill only urchins, exactly as primary processing; no lobster fallback.
 ix<-match(sf$site,mm$site);i<-which(is.na(sf$MPA)&sf$classcode!='PANINT');sf[i,MPA:=mm$MPA[ix[i]]];i<-which(is.na(sf$status)&sf$classcode!='PANINT');sf[i,status:=mm$status[ix[i]]]
 sf[,Taxa:=ifelse(classcode=='PANINT','P. interruptus','S. purpuratus')]
 sf<-sf[!is.na(size)&!is.na(count)&count>0 & status %in% c('mpa','reference')]
 stopifnot(all(sf$count==floor(sf$count)),all(sf$size>0))
 poolrows<-rbindlist(lapply(c('PISCO','KFM'),function(src){x<-copy(sf[if(src=='PISCO') Taxa=='P. interruptus'|size>=25 else Taxa=='S. purpuratus']);x[,Source:=src];x}))
 members<-merge(annual[value>0 & !(Source=='LTER' & Taxa=='P. interruptus')],poolrows,by=c('MPA','Taxa','Source','year','status'),allow.cartesian=TRUE)
 stopifnot(!anyDuplicated(members[,.(series_id,year,status,source_row)]))
 members[,mass:=ifelse(Taxa=='P. interruptus',0.001352821*size^2.913963113,0.00059*size^2.870)]
 pools<-members[,.(n_sized=sum(count),size_records=.N,n_size_sites=uniqueN(site),n_size_classes=uniqueN(size),size_sites=paste(sort(unique(site)),collapse=';'),size_programs=paste(sort(unique(campus)),collapse=';'),mean_mass=sum(count*mass)/sum(count),largest_individual_mass_fraction=max(mass)/sum(count*mass)),by=.(series_id,year,status)]
 inventory<-merge(annual,pools,by=c('series_id','year','status'),all.x=TRUE)
 inventory[,method:=ifelse(Source=='LTER','direct measured lengths','resampled size pool')]
 stopifnot(!anyNA(inventory[method=='resampled size pool'&value>0,n_sized]))
 # Direct LTER lobster: report measured count and reuse, not a reconstructed size pool.
 dl<-read('data/LTER/Lobster_Abundance_All_Years_20230831.csv',dr);dl[,source_row:=.I]
 setnames(dl,c('SITE','YEAR'),c('site','year'));dl<-merge(dl,mm,by='site',allow.cartesian=TRUE)
 dmembers<-merge(annual[Source=='LTER'],dl,by=c('MPA','year','status'),allow.cartesian=TRUE)
 direct<-dmembers[,.(n_sized=sum(COUNT[COUNT>0 & is.finite(SIZE_MM)&SIZE_MM>0]),size_records=sum(COUNT>0 & is.finite(SIZE_MM)&SIZE_MM>0),n_size_sites=uniqueN(site[COUNT>0 & is.finite(SIZE_MM)&SIZE_MM>0]),n_size_classes=uniqueN(SIZE_MM[COUNT>0 & is.finite(SIZE_MM)&SIZE_MM>0]),size_sites=paste(sort(unique(site[COUNT>0 & is.finite(SIZE_MM)&SIZE_MM>0])),collapse=';'),positive_count_missing_size=sum(COUNT[COUNT>0 & (!is.finite(SIZE_MM)|SIZE_MM<=0)])),by=.(series_id,year,status)]
 for(n in c('n_sized','size_records','n_size_sites','n_size_classes','size_sites'))inventory[direct,on=.(series_id,year,status),(n):=get(paste0('i.',n))]
 inventory[direct,on=.(series_id,year,status),positive_count_missing_size:=i.positive_count_missing_size]
 # Validate direct hierarchy against the published primary annual values.
 dmembers[,biomass:=COUNT*0.001352821*SIZE_MM^2.913963113];dmembers[!is.finite(biomass),biomass:=0]
 if(file.exists(file.path(root,'data/harmonized/lobster_missing_lengths.json'))) {
   donors<-dmembers[COUNT>0 & SIZE_MM>0,.(mean_mass=sum(biomass)/sum(COUNT)),by=.(series_id,site,DATE)]
   dmembers<-merge(dmembers,donors,by=c('series_id','site','DATE'),all.x=TRUE)
   stopifnot(all(is.finite(dmembers[COUNT>0 & (!is.finite(SIZE_MM)|SIZE_MM<=0),mean_mass])))
   dmembers[COUNT>0 & (!is.finite(SIZE_MM)|SIZE_MM<=0),biomass:=COUNT*mean_mass]
 }
 rec<-dmembers[,.(b=sum(biomass)),by=.(series_id,year,status,site,TRANSECT,REPLICATE)][,.(b=mean(b)),by=.(series_id,year,status,site,TRANSECT)][,.(b=mean(b)),by=.(series_id,year,status,site)][,.(reconstructed=mean(b)/300),by=.(series_id,year,status)]
 rec<-merge(rec,annual,by=c('series_id','year','status'));stopifnot(max(abs(rec$reconstructed-rec$value))<1e-10)
 # Shared original records across retained histories (count grouped animals once).
 use<-rbindlist(list(members[,.(series_id,year,status,source_file='ALL_sizefreq_2024.csv',source_row,count,site)],dmembers[COUNT>0 & SIZE_MM>0,.(series_id,year,status,source_file='Lobster_Abundance_All_Years_20230831.csv',source_row,count=COUNT,site)]))
 use<-unique(use)
 overlap<-merge(use,use,by=c('source_file','source_row'),allow.cartesian=TRUE,suffixes=c('_a','_b'))[series_id_a<series_id_b]
 stopifnot(all(overlap$count_a==overlap$count_b),all(overlap$year_a==overlap$year_b))
 shared<-overlap[,.(shared_records=.N,shared_measured_animals=sum(count_a),shared_sites=paste(sort(unique(site_a)),collapse=';')),by=.(series_id_a,series_id_b,year=year_a,status_a,status_b)]
 flags<-unique(rbind(shared[,.(series_id=series_id_a,year,status=status_a)],shared[,.(series_id=series_id_b,year,status=status_b)]));flags[,shared_with_other_history:=TRUE]
 inventory<-merge(inventory,flags,by=c('series_id','year','status'),all.x=TRUE);inventory[is.na(shared_with_other_history),shared_with_other_history:=FALSE]
 inventory[,`:=`(positive=value>0,sparse_lt10=value>0 & !is.na(n_sized)&n_sized<10,singleton=value>0 & !is.na(n_sized)&n_sized==1)]
 summary<-inventory[,.(paired_years=uniqueN(year),positive_side_years=sum(positive),min_measured_positive=min(n_sized[positive]),median_measured_positive=as.numeric(median(n_sized[positive])),singleton_side_years=sum(singleton),under10_side_years=sum(sparse_lt10),under20_side_years=sum(positive&!is.na(n_sized)&n_sized<20),shared_side_years=sum(shared_with_other_history),sparse_years=paste(sort(unique(year[sparse_lt10])),collapse=',')),by=.(series_id,MPA,Taxa,Source,method)]
 summary[,priority:=fcase(singleton_side_years>0,'Singleton size samples',under10_side_years>0,'Some small samples',shared_side_years>0,'Shared measurements',default='No sparse/shared flag at these thresholds')]
 write<-function(d,n)fwrite(d,file.path(root,paste0('tables/table_s_primary_size_',n,'.csv')))
 private<-file.path(root,'local_archive/primary_size_inventory_2026-10-30');dir.create(private,recursive=TRUE,showWarnings=FALSE)
 write_private<-function(d,n)fwrite(d,file.path(private,paste0(n,'.csv')))
 write_private(dmembers[COUNT>0 & (!is.finite(SIZE_MM)|SIZE_MM<=0),.(series_id,year,status,site,TRANSECT,REPLICATE,source_row,COUNT,SIZE_MM)],'missing_lengths');
 write(inventory,'inventory');write(summary,'comparisons');write(shared,'sharing');write_private(use,'record_membership');write(rec,'direct_validation')
 used<-unique(c(used,file.path(dr,paste0('code/R/',c('01_utils','03_data_import','04_pisco_processing','05_kfm_processing','06_lter_processing'),'.R')),file.path(root,'code/R/55_primary_size_inventory.R')))
 write(data.table(repository=ifelse(startsWith(used,paste0(dr,'/')),'data','analysis'),path=ifelse(startsWith(used,paste0(dr,'/')),substring(used,nchar(dr)+2),substring(used,nchar(root)+2)),md5=unname(tools::md5sum(used))),'inputs')
 stopifnot(identical(before,tools::md5sum(primary_files)),nrow(summary)==nrow(keys))
 print(summary);cat('Inventory only; primary data/results unchanged. Thresholds are descriptive flags, not exclusion criteria.\n')
})
