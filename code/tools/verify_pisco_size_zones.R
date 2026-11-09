source('code/R/primary_uncertainty_utils.R')
o<-'outputs/pisco_size_zone_resolution_20261102/'
a<-read.csv('outputs/biomass_resolution_20261102/review_annual_input.csv');e<-read.csv(paste0(o,'effects.csv'));z<-read.csv(paste0(o,'annual_sensitivities.csv'))
for(i in seq_len(nrow(e))){
 if(e$scenario[i]=='current_all_years'){q<-a[a$series_id==e$series_id[i]&a$t>0,];q$y<-q$lnDiff}else{
  sc<-strsplit(e$scenario[i],'__',fixed=TRUE)[[1]];q<-z[z$series_id==e$series_id[i]&z$scenario==sc[1],];q$y<-q[[if(sc[2]=='current_matched_years')'original' else 'revised']]
 }
 f<-lm(y~t,q);stopifnot(nrow(q)==e$n_years[i],abs(coef(f)[2]-e$estimate[i])<1e-9,abs(sqrt(uncertainty_hac(f,q$year,3)[2,2])-e$se[i])<1e-9)
}
cat('Independently verified',nrow(e),'matched-year slopes and temporal standard errors.\n')
