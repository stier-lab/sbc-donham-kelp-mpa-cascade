# Add source-mapped Anacapa SMCA size memberships absent from the original metadata membership join.
# Input SizeFreq is the cleaned pre-join source table from data module03.
assignments<-Sites2%>%filter(CA_MPA_Name_Short=='Anacapa Island SMCA',data_source=='kfm')%>%select(site_id,CA_MPA_Name_Short,status)%>%distinct()
stopifnot(!anyDuplicated(assignments$site_id))
SizeFreq.Urch.OG%>%filter(campus=='KFM',site%in%assignments$site_id)%>%count(site,CA_MPA_Name_Short,site_status,classcode)%>%write.csv('../size_mapping_before.csv',row.names=FALSE)
z<-SizeFreq%>%filter(campus=='KFM',site%in%assignments$site_id,classcode%in%c('MESFRA','STRPUR'))
z$classcode<-ifelse(z$classcode=='MESFRA','MESFRAAD','STRPURAD')
z$.source_row_id<-seq_len(nrow(z))
z<-merge(z,assignments,by.x='site',by.y='site_id',all.x=TRUE,sort=FALSE)
names(z)[names(z)=='status']<-'site_status'
stopifnot(!anyDuplicated(z$.source_row_id),all(!is.na(z$CA_MPA_Name_Short)),all(z$count>=0),all(z$size>0))
write.csv(z,'../restored_size_assignments.csv',row.names=FALSE)
z$.source_row_id<-NULL
SizeFreq.Urch.OG<-SizeFreq.Urch.OG%>%filter(!(campus=='KFM'&site%in%assignments$site_id&!is.na(CA_MPA_Name_Short)&CA_MPA_Name_Short=='Anacapa Island SMCA'))
SizeFreq.Urch.OG<-bind_rows(SizeFreq.Urch.OG,z)
saveRDS(SizeFreq.Urch.OG,'../kfm_urchin_size_pool.rds')
cat('Restored Anacapa source size assignments:',nrow(z),'rows at',length(unique(z$site)),'stations\n')
