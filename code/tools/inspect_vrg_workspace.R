b <- '/Users/adrianstier/Library/CloudStorage/GoogleDrive-<account>/My Drive/Stier Lab/People/Emily Donham/Projects/Kelp MPA/Github/MPA_Kelp'
e<-new.env(parent=emptyenv());n<-load(file.path(b,'.RData'),envir=e)
z<-do.call(rbind,lapply(n,function(k){x<-get(k,e);data.frame(object=k,class=class(x)[1],rows=if(is.data.frame(x))nrow(x) else NA,columns=if(is.data.frame(x))paste(names(x),collapse=';') else '')}))
write.csv(z,'outputs/vrg_source_search_20261102/workspace_inventory.csv',row.names=FALSE)
print(z[grepl('Size|Freq|VRG',z$object),c('object','class','rows')])
for(k in c('SizeFreq','SizeFreq.south','SizeFreq.PANINT.VRG','SizeFreq.Urch.OG')) {
 x<-get(k,e); if('campus'%in%names(x))x<-x[!is.na(x$campus)&x$campus=='VRG',,drop=FALSE]
 write.csv(x,paste0('outputs/vrg_source_search_20261102/workspace_',k,'.csv'),row.names=FALSE)
}
t<-new.env(parent=emptyenv());err<-tryCatch({nn<-load(file.path(b,'.RDataTmp'),envir=t);paste('Loaded temporary workspace:',paste(nn,collapse=', '))},error=function(e)paste('Temporary workspace not readable:',conditionMessage(e)))
writeLines(err,'outputs/vrg_source_search_20261102/temporary_workspace_status.txt')
