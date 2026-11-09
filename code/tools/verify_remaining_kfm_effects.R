p<-'outputs/kfm_remaining_sources_20261102';a<-read.csv(file.path(p,'scorpion_primary_annual_sensitivity.csv'));e<-read.csv(file.path(p,'scorpion_primary_step_sensitivity.csv'))
for(i in 1:2){a$yy<-a[[c('lnDiff','legacy')[i]]];a$after<-as.integer(a$BA=='After');f<-lm(yy~after,a);stopifnot(abs(coef(f)[2]-e$estimate[i])<1e-10,abs(summary(f)$coef[2,2]-e$se[i])<1e-10,abs(confint(f)[2,1]-e$low[i])<1e-10)}
cat('Both Scorpion step estimates, SEs and intervals independently verified.\n')
