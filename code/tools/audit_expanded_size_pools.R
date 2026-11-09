# Expanded inventory; retain raw size records and primary calculations.
suppressPackageStartupMessages(library(data.table))
O<-'outputs/biomass_support_20261102';D<-'/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026';private<-file.path(O,'private');dir.create(private,showWarnings=FALSE);writeLines('private/',file.path(O,'.gitignore'))
v<-fread('outputs/all_program_station_refits_20261102/verified_station_values.csv')[Resp=='Bio' & Taxa %in% c('P. interruptus','S. purpuratus','M. franciscanus')]
annual<-v[,.(count_sites=paste(sort(unique(site[is.finite(value)])),collapse=';'),value=mean(value,na.rm=TRUE)),by=.(series_id,MPA,Taxa,Source,year,status)]
expr<-as.list(parse(file.path(D,'code/R/01_utils.R')));constant<-function(n)eval(Filter(function(x)is.call(x)&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name(n)),expr)[[1]][[3]],baseenv())
sf<-fread(file.path(D,'data/ALL_sizefreq_2024.csv'));sf[,source_row:=.I]
sf<-sf[(campus %in% c('VRG','UCSB') & classcode %in% c('STRPUR','MESFRA','PANINT'))|(campus=='KFM' & classcode %in% c('Strongylocentrotus purpuratus','Strongylocentrotus franciscanus'))|(campus=='LTER' & classcode %in% c('SPS','SPL','SFS','SFL'))]
sf<-sf[!site %in% c(constant('EXCLUDED_REFERENCE_SITES'),constant('EXCLUDED_KFM_SITES'))]
sf[,Taxa:=fcase(classcode=='PANINT','P. interruptus',classcode %in% c('MESFRA','Strongylocentrotus franciscanus','SFS','SFL'),'M. franciscanus',default='S. purpuratus')]
sm<-fread(file.path(D,'data/PISCO/master_site_table_Emilyedit.csv'));sm<-sm[!duplicated(sm[,.(site,CA_MPA_Name_Short)]),.(site,site_new,campus,MPA=CA_MPA_Name_Short,status=site_status)]
sf<-merge(sf,sm,by=c('site','site_new','campus'),all.x=TRUE,allow.cartesian=TRUE)
mm<-unique(fread(file.path(D,'data/MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv'))[,.(site=site_id,MPA=CA_MPA_Name_Short,status)])
ix<-match(sf$site,mm$site);idx<-which(is.na(sf$MPA)&sf$Taxa!='P. interruptus');sf[idx,MPA:=mm$MPA[ix[idx]]];idx<-which(is.na(sf$status)&sf$Taxa!='P. interruptus');sf[idx,status:=mm$status[ix[idx]]]
sf<-sf[is.finite(size)&size>0 & is.finite(count)&count>0 & status %in% c('mpa','reference')];stopifnot(all(sf$count==floor(sf$count)))
poolrows<-rbindlist(lapply(c('PISCO','KFM','LTER'),function(src){x<-copy(sf[if(src=='PISCO') Taxa=='P. interruptus'|size>=25 else Taxa!='P. interruptus']);x[,Source:=src];x}))
members<-merge(annual[value>0 & !(Source=='LTER'&Taxa=='P. interruptus')],poolrows,by=c('MPA','Taxa','Source','year','status'),allow.cartesian=TRUE)
stopifnot(!anyDuplicated(members[,.(series_id,year,status,source_row)]))
members[,mass:=fcase(Taxa=='P. interruptus',0.001352821*size^2.913963113,Taxa=='S. purpuratus',0.00059*size^2.870,default=NA_real_)]
# The exact red-urchin allometry is not needed to count support; do not invent it.
members[,full_date:=!is.na(month)&month>=1&month<=12&!is.na(day)&day>=1&day<=31]
members[,same_program:=fifelse(Source=='PISCO',campus %in% c('UCSB','VRG'),campus==Source)]
members[,same_count_site:=mapply(function(s,ss)s %in% strsplit(ss,';',fixed=TRUE)[[1]],site,count_sites)]
members[,survey_group:=do.call(paste,c(.SD,sep='|')),.SDcols=c('campus','method','site','year','month','day','location','zone','transect')]
pools<-members[,.(n_sized=sum(count),size_records=.N,n_size_sites=uniqueN(site),n_size_classes=uniqueN(size),n_recorded_groups=uniqueN(survey_group),size_sites=paste(sort(unique(site)),collapse=';'),size_programs=paste(sort(unique(campus)),collapse=';'),fraction_same_program=sum(count[same_program])/sum(count),fraction_same_count_sites=sum(count[same_count_site])/sum(count),fraction_complete_dates=sum(count[full_date])/sum(count),size_months=paste(sort(unique(month[month>=1&month<=12])),collapse=';'),mean_mass=if(all(is.finite(mass)))sum(count*mass)/sum(count) else NA_real_),by=.(series_id,year,status)]
inventory<-merge(annual,pools,by=c('series_id','year','status'),all.x=TRUE);inventory[,method:=ifelse(Source=='LTER'&Taxa=='P. interruptus','direct measured lengths','annual MPA/role/species pooled sizes')]
# Direct LTER lobster measurement support.
dl<-fread(file.path(D,'data/LTER/Lobster_Abundance_All_Years_20230831.csv'));dl[,source_row:=.I];setnames(dl,c('SITE','YEAR'),c('site','year'));dl<-merge(dl,mm,by='site',allow.cartesian=TRUE)
dmembers<-merge(annual[Source=='LTER'&Taxa=='P. interruptus'],dl,by=c('MPA','year','status'),allow.cartesian=TRUE)
direct<-dmembers[,.(n_sized=sum(COUNT[COUNT>0 & SIZE_MM>0 & is.finite(SIZE_MM)]),size_records=sum(COUNT>0 & SIZE_MM>0 & is.finite(SIZE_MM)),n_size_sites=uniqueN(site[COUNT>0&is.finite(SIZE_MM)&SIZE_MM>0]),n_size_classes=uniqueN(SIZE_MM[COUNT>0&is.finite(SIZE_MM)&SIZE_MM>0]),positive_count_missing_size=sum(COUNT[COUNT>0&(!is.finite(SIZE_MM)|SIZE_MM<=0)]),size_sites=paste(sort(unique(site[COUNT>0&is.finite(SIZE_MM)&SIZE_MM>0])),collapse=';')),by=.(series_id,year,status)]
for(n in setdiff(names(direct),c('series_id','year','status')))inventory[direct,on=.(series_id,year,status),(n):=get(paste0('i.',n))]
use<-unique(rbindlist(list(members[,.(series_id,year,status,source_file='ALL_sizefreq_2024.csv',source_row,count,site)],dmembers[COUNT>0&is.finite(SIZE_MM)&SIZE_MM>0,.(series_id,year,status,source_file='Lobster_Abundance_All_Years_20230831.csv',source_row,count=COUNT,site)])))
shared_ids<-use[,.(n_histories=uniqueN(series_id)),by=.(source_file,source_row)][n_histories>1];flags<-unique(merge(use,shared_ids,by=c('source_file','source_row'))[,.(series_id,year,status)]);flags[,shared_with_other_history:=TRUE];inventory<-merge(inventory,flags,by=c('series_id','year','status'),all.x=TRUE);inventory[is.na(shared_with_other_history),shared_with_other_history:=FALSE]
inventory[,`:=`(positive=value>0,sparse_lt10=value>0 & !is.na(n_sized)&n_sized<10,singleton=value>0 & !is.na(n_sized)&n_sized==1,missing_size_support=value>0 & (is.na(n_sized)|n_sized==0))]
summary<-inventory[,.(paired_years=uniqueN(year),positive_side_years=sum(positive),singleton_side_years=sum(singleton),under10_side_years=sum(sparse_lt10),missing_support_side_years=sum(missing_size_support),shared_side_years=sum(shared_with_other_history&positive),cross_program_side_years=sum(positive&fraction_same_program<1,na.rm=TRUE),outside_count_footprint_side_years=sum(positive&fraction_same_count_sites<1,na.rm=TRUE),incomplete_size_date_side_years=sum(positive&fraction_complete_dates<1,na.rm=TRUE),min_n_sized=if(any(positive&!is.na(n_sized)))min(n_sized[positive&!is.na(n_sized)]) else NA_real_),by=.(series_id,MPA,Taxa,Source,method)]
fwrite(inventory,file.path(O,'animal_pool_inventory.csv'));fwrite(summary,file.path(O,'animal_pool_summary.csv'));fwrite(members,file.path(private,'pool_members.csv'));fwrite(use,file.path(private,'record_membership.csv'));fwrite(dmembers,file.path(private,'direct_lobster_members.csv'))
# Compare overlap to earlier fifteen-history inventory; discrepancies are reported, not overwritten.
old<-fread('tables/table_s_primary_size_inventory.csv');check<-merge(inventory,old,by=c('series_id','year','status'),suffixes=c('_new','_old'));check[,n_sized_equal:=(is.na(n_sized_new)&is.na(n_sized_old))|(!is.na(n_sized_new)&!is.na(n_sized_old)&n_sized_new==n_sized_old)];fwrite(check[,.(series_id,year,status,n_sized_new,n_sized_old,n_sized_equal)],file.path(O,'prior_inventory_check.csv'))
print(summary);cat('Rows:',nrow(inventory),'histories:',nrow(summary),'prior mismatches:',sum(!check$n_sized_equal),'\n')
