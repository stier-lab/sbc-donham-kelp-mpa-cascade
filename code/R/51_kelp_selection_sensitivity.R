# Targeted contrast-selection check for two near-threshold PISCO kelp trends.
# Uses only current inputs; the two mean substitutions are diagnostic, not a
# replacement primary analysis. Shared variance can shift other biomass taxa.
source('code/R/00_libraries.R')
new <- read.csv('data/sumstats_final.csv')
snap <- readRDS('data/cache/figures_snapshot.rds')
stopifnot(identical(snap$harmonized_md5, unname(tools::md5sum(file.path('data/harmonized',
 c('harmonized_raw_responses.csv','harmonized_response_ratios.csv','harmonized_landsat_rr.csv','harmonized_landsat_annual_index.csv'))))))
dat <- snap$All.RR.sub.trans
switched <- which(new$MPA %in% c('Naples SMCA','South Point SMR') & new$Taxa=='M. pyrifera' & new$Source=='PISCO' & new$Resp=='Bio');stopifnot(length(switched)==2);fixed <- new
checks <- list()
for(i in switched){
 stopifnot(new$Model[i]=='Linear',new$Source[i]=='PISCO')
 d <- subset(dat,CA_MPA_Name_Short==new$MPA[i] & y=='Macrocystis pyrifera' & source=='PISCO' & resp=='Bio')
 # Loader uses abbreviated taxon names in this object.
 if(!nrow(d)) d <- subset(dat,CA_MPA_Name_Short==new$MPA[i] & y=='M. pyrifera' & source=='PISCO' & resp=='Bio')
 stopifnot(nrow(d)==new$N[i])
 fixed$Mean[i] <- mean(d$lnDiff);fixed$SE[i] <- sd(d$lnDiff)/sqrt(nrow(d))
 checks[[length(checks)+1]] <- data.frame(MPA=new$MPA[i],n=nrow(d),slope_p=coef(summary(lm(lnDiff~time,d)))[2,4],mean_contrast=fixed$Mean[i],mean_se=fixed$SE[i],selected_contrast=new$Mean[i])
}
checks <- bind_rows(checks)
stopifnot(all(checks$slope_p > .045 & checks$slope_p <= .05))
write.csv(checks,'tables/table_s_fixed_kelp_selection.csv',row.names=FALSE)
results <- list()
for(scenario in c('primary','two_kelp_means')){
 d <- if(scenario=='primary')new else fixed
 for(resp in c('Bio','Den')){
  x <- subset(d,Resp==resp);x$Taxa <- factor(x$Taxa)
  m <- rma.mv(Mean,SE^2,mods=~Taxa-1,random=list(~1|MPA,~1|Source),data=x,method='REML',test='t')
  results[[length(results)+1]] <- data.frame(scenario,Taxa=sub('^Taxa','',rownames(m$b)),Resp=resp,estimate=c(m$b),se=m$se,ci_low=m$ci.lb,ci_high=m$ci.ub,p=m$pval)
 }
}
r <- bind_rows(results) %>% group_by(scenario) %>% mutate(q=p.adjust(p,'BH')) %>% ungroup()
primary <- read.csv('tables/table_02_meta_analysis.csv')
observed <- subset(r,scenario=='primary')
key <- paste(primary$Taxa,ifelse(primary$Response=='Biomass','Bio','Den'))
ix <- match(paste(observed$Taxa,observed$Resp),key)
stopifnot(!anyNA(ix),max(abs(observed$estimate-primary$Estimate[ix])) < 1e-7,
 max(abs(observed$q-primary$pval_fdr[ix])) < 1e-7)
write.csv(r,'tables/table_s_fixed_kelp_selection_pooled.csv',row.names=FALSE)
cat('Near-threshold kelp contrast-selection check complete.\n')
