source('code/R/primary_uncertainty_utils.R')
a<-read.csv('outputs/biomass_resolution_20261102/review_annual_input.csv');r<-read.csv('outputs/kfm_scott_reconciliation_20261102/all_corrections.csv');e<-read.csv('outputs/kfm_scott_reconciliation_20261102/effects.csv')
for(k in seq_len(nrow(e))){z<-a[a$series_id==e$series_id[k]&a$t>0,];rr<-r[r$series_id==e$series_id[k],];if(e$scenario[k]=='historical_only')rr<-rr[rr$year<=2011,];if(e$scenario[k]!='original_review'){j<-match(z$year,rr$year);keep<-!is.na(j);z$lnDiff[keep]<-rr$new_lnDiff[j[keep]]};f<-lm(lnDiff~t,z);stopifnot(abs(coef(f)[2]-e$estimate[k])<1e-10,abs(sqrt(uncertainty_hac(f,z$year,3)[2,2])-e$se[k])<1e-10)}
cat(nrow(e),'Scott reconciliation fits independently verified in R.\n')
