source('code/R/primary_uncertainty_utils.R')
a<-read.csv('outputs/biomass_resolution_20261102/review_annual_input.csv');r<-read.csv('outputs/biomass_protocols_20261102/sbi_size_replacements.csv');e<-read.csv('outputs/biomass_protocols_20261102/sbi_reconciliation_effects.csv')
for(k in seq_len(nrow(e))){z<-a[a$series_id==e$series_id[k]&a$t>0,];if(e$scope[k]=='2005_2011_same_years')z<-z[z$year>=2005&z$year<=2011,];if(e$scenario[k]!='current'){j<-match(paste(z$series_id,z$year),paste(r$series_id,r$year));keep<-!is.na(j);z$lnDiff[keep]<-r$new_lnDiff[j[keep]]};f<-lm(lnDiff~t,z);stopifnot(abs(coef(f)[2]-e$estimate[k])<1e-10,abs(sqrt(uncertainty_hac(f,z$year,3)[2,2])-e$se[k])<1e-10)}
cat('8 SBI reconciliation fits independently verified in R.\n')
