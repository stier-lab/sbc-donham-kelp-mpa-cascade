suppressPackageStartupMessages({library(metafor);library(clubSandwich)})
for(expr in as.list(parse('code/tools/refit_common_estimand.R')))if(is.call(expr)&&identical(expr[[1]],as.name('<-'))&&identical(expr[[2]],as.name('pool')))eval(expr)
out<-'outputs/seasonal_zone_audit_20261101';e<-read.csv(file.path(out,'effects.csv'));logs<-list();rows<-list()
for(s in c('Original','Fixed_all_program_units')) {
 x<-e[e$scenario==s&e$metric=='Post_slope'&e$Taxa=='P. interruptus'&e$Resp=='Den',]
 # Diagnostic only: remove the one previously identified three-year history.
 x<-x[!(x$MPA=='Point Dume SMCA'&x$Source=='PISCO'),];x$scenario<-paste0(s,'_without_Dume_lobster');rows[[length(rows)+1]]<-pool(x,'All')
}
write.csv(do.call(rbind,rows),file.path(out,'lobster_influence_diagnostic.csv'),row.names=FALSE)
