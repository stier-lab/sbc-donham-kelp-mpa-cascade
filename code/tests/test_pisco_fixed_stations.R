# Read-only invariants for the fixed-station sensitivity products.
suppressPackageStartupMessages({library(dplyr);library(tidyr)})
r <- function(n) read.csv(paste0('tables/table_s_pisco_fixed_',n,'.csv'))
cov<-r('coverage');annual<-r('annual');sites<-r('station_values');fits<-r('fits');groups<-r('feasible_groups')
stopifnot(nrow(cov)==12,sum(cov$eligible)==10,max(abs(r('reconstruction')$difference))<1e-7)
for(i in which(cov$eligible)) {
 c<-cov[i,];yy<-as.integer(strsplit(c$years,';')[[1]]);ss<-strsplit(paste(c$protected,c$reference,sep=';'),';')[[1]]
 d<-filter(sites,scenario==c$scenario,MPA==c$MPA)
 stopifnot(nrow(d)==length(yy)*length(ss),!anyDuplicated(d[c('site','year')]))
 stopifnot(all(table(d$year)==length(ss)),all(table(d$site)==length(yy)))
 a<-filter(annual,scenario==c$scenario,MPA==c$MPA)
 stopifnot(nrow(a)==length(yy)*4,all(is.finite(a$lnDiff)))
 for(resp in c('Bio','Den')) stopifnot(identical(sort(a$year[a$Resp==resp&a$membership=='fixed']),sort(a$year[a$Resp==resp&a$membership=='mixture_matched_years'])))
 boundary<-strsplit(c$scenario,'_min')[[1]][1];minimum<-as.integer(strsplit(c$scenario,'_min')[[1]][2])
 g<-groups[groups$MPA==c$MPA&groups[[paste0(boundary,'_before')]]>=minimum,]
 g<-arrange(g,desc(n_stations),desc(n_years),protected,reference)
 stopifnot(c$protected==g$protected[1],c$reference==g$reference[1],c$after>=4)
}
# CI and pBACIPS use exactly the same observed years/endpoints; mean-level estimates
# are explicitly distinguished from temporal contrasts by model labels.
for(i in seq_len(nrow(fits))) {
 f<-fits[i,];a<-filter(annual,scenario==f$scenario,MPA==f$MPA,Resp==f$Resp,membership==f$membership)
 stopifnot(f$N==nrow(a),f$last_contrast %in% a$year,f$last_contrast<=f$first_after+11,
   f$HAC_SE>=f$SE-1e-10,is.finite(f$Mean),f$SE>0)
 if(f$method=='pBACIPS')stopifnot(f$pretrend_p>=.05)
 a<-a[order(a$year),];a$t<-pmax(a$year-2012,0)
 a$BA<-factor(ifelse(if(grepl('^strict',f$scenario))a$year<2012 else a$year<=2012,'Before','After'),c('Before','After'))
 if(f$Model %in% c('Mean','Step','Linear')) {
   m<-switch(f$Model,Mean=lm(lnDiff~1,a),Step=lm(lnDiff~BA,a),Linear=lm(lnDiff~t,a))
   g<-switch(f$Model,Mean=1,Step=c(0,1),Linear=c(0,f$last_contrast-f$first_after))
   stopifnot(abs(sum(g*coef(m))-f$Mean)<1e-8,abs(sqrt(drop(t(g)%*%vcov(m)%*%g))-f$SE)<1e-8)
 }
 if(f$method=='CI') {
   lm0<-lm(lnDiff~t,a);expected<-if(summary(lm0)$coefficients[2,4]<=.05)'Linear' else 'Mean'
   stopifnot(f$Model==expected)
 }
}
# All four alternatives preserve the sheephead biomass result; neither sheephead
# pathway is resolved. This guards the short SI claims against future result drift.
p<-r('pooled');z<-p[grepl('\\|fixed\\|',p$scenario)&p$Taxa=='S. pulcher'&p$Resp=='Bio',]
stopifnot(nrow(z)==16,all(z$estimate>0),all(z$q<.05))
z<-r('pathways');z<-z[grepl('\\|fixed\\|',z$scenario)&z$panel %in% c('c','d'),]
stopifnot(nrow(z)==32,all(z$ci_low<0&z$ci_high>0))
b<-read.csv('tables/table_02_meta_analysis.csv') %>% transmute(Taxa,Resp=ifelse(Response=='Biomass','Bio','Den'),base=Estimate)
j<-inner_join(p,b,by=c('Taxa','Resp'));stopifnot(all(sign(j$estimate)==sign(j$base)))
b<-read.csv('tables/table_s_routing_clock_pathways.csv') %>% filter(scenario=='retained_primary',panel %in% c('a','b')) %>% select(panel,uncertainty,base_slope=slope,base_p=p)
j<-inner_join(r('pathways'),b,by=c('panel','uncertainty'))
stopifnot(max(abs(j$slope-j$base_slope))<1e-10,max(abs(j$p-j$base_p))<1e-10)
h<-r('inputs');stopifnot(identical(unname(tools::md5sum(h$path)),h$md5))
cat('PASS: reconstruction, complete station/year balance, coverage-only choice, matched models, observed endpoints, inference claims and unchanged source inputs.\n')
