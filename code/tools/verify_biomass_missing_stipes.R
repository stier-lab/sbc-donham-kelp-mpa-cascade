source('code/R/primary_uncertainty_utils.R');a<-read.csv('outputs/pisco_fish_rules_20261101/approved_annual_input.csv');c<-read.csv('outputs/biomass_support_20261102/missing_stipe_sensitivity.csv')
for(i in seq_len(nrow(c))){z<-a[a$series_id==c$series_id[i]&a$t>0,];if(c$scenario[i]!='Current')z<-z[z$year!=2023,];z<-z[order(z$year),];f<-lm(lnDiff~t,z);se<-sqrt(uncertainty_hac(f,z$year,3)[2,2]);stopifnot(abs(coef(f)[2]-c$estimate[i])<1e-10,abs(se-c$se[i])<1e-10)}
cat('Two missing-stipe fits independently verified\n')
