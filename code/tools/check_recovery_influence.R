suppressPackageStartupMessages({library(metafor);library(clubSandwich)})
for(expr in as.list(parse('code/tools/refit_common_estimand.R')))if(is.call(expr)&&identical(expr[[1]],as.name('<-'))&&identical(expr[[2]],as.name('pool')))eval(expr)
out<-'outputs/full_followup_recovery_20261101';e<-read.csv(file.path(out,'effects.csv'));logs<-list();rows<-list()
for(mt in c('Post_slope','Change'))for(tx in unique(e$Taxa))for(resp in c('Bio','Den')){
 x<-e[e$scenario=='Full'&e$metric==mt&e$Taxa==tx&e$Resp==resp,];if(!nrow(x))next
 for(m in unique(x$MPA)){z<-x[x$MPA!=m,];if(!nrow(z))next;z$scenario<-'Full_leave_one_MPA';r<-pool(z,'All');r$omitted_MPA<-m;rows[[length(rows)+1]]<-r}
}
write.csv(do.call(rbind,rows),file.path(out,'leave_one_MPA.csv'),row.names=FALSE)
rows<-list()
for(tx in unique(e$Taxa))for(resp in c('Bio','Den'))for(src in c('KFM','LTER','PISCO','Landsat')){
 x<-e[e$scenario=='Full'&e$metric=='Post_slope'&e$Taxa==tx&e$Resp==resp&e$Source==src,];if(!nrow(x))next;r<-pool(x,src);rows[[length(rows)+1]]<-r
}
write.csv(do.call(rbind,rows),file.path(out,'program_slopes.csv'),row.names=FALSE)
cat('Saved leave-one-MPA and separate-program descriptive checks; no BH claim from these diagnostics.\n')
