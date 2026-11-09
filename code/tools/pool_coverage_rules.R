suppressPackageStartupMessages({library(metafor);library(clubSandwich)})
for(expr in as.list(parse('code/tools/refit_common_estimand.R')))if(is.call(expr)&&identical(expr[[1]],as.name('<-'))&&identical(expr[[2]],as.name('pool')))eval(expr)
out<-'outputs/point_dume_coverage_20261101';e<-read.csv(file.path(out,'effects.csv'));logs<-list();rows<-list()
for(s in unique(e$scenario))for(mt in unique(e$metric))for(tx in unique(e$Taxa))for(r in c('Bio','Den')){
 x<-e[e$scenario==s&e$metric==mt&e$Taxa==tx&e$Resp==r,];if(nrow(x))rows[[length(rows)+1]]<-pool(x,'All')
}
p<-do.call(rbind,rows);p$q<-ave(p$p,interaction(p$scenario,p$metric),FUN=function(x)p.adjust(x,'BH',n=max(9,sum(is.finite(x)))))
write.csv(p,file.path(out,'pooled.csv'),row.names=FALSE);write.csv(do.call(rbind,logs),file.path(out,'fit_log.csv'),row.names=FALSE)
f<-read.csv(file.path(out,'footprint_ratios.csv'));r<-list()
for(s in unique(f$scenario))for(tx in unique(f$Taxa)){
 x<-f[f$scenario==s&f$Taxa==tx,];m<-lm(ln_ratio~year,x);r[[length(r)+1]]<-data.frame(scenario=s,Taxa=tx,n=nrow(x),slope=unname(coef(m)[2]),ratio_2021=exp(x$ln_ratio[x$year==2021]),ratio_2023=exp(x$ln_ratio[x$year==2023]))
}
write.csv(do.call(rbind,r),file.path(out,'footprint_slopes.csv'),row.names=FALSE)
old<-read.csv('outputs/full_followup_recovery_20261101/pooled.csv');old<-old[old$scenario=='Full',];j<-merge(p[p$scenario=='All_three_year',],old,by=c('metric','Taxa','Resp'));stopifnot(nrow(j)==nrow(old),max(abs(j$estimate.x-j$estimate.y))<1e-7,max(abs(j$q.x-j$q.y),na.rm=TRUE)<1e-7)
cat('Original full-follow-up results reproduced. Saved',nrow(p),'pooled rule results.\n')
