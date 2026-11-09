source('code/R/primary_uncertainty_utils.R')
o<-'outputs/pisco_all_months_20261102/'
e<-read.csv(paste0(o,'effects.csv'));a<-read.csv(paste0(o,'annual_values.csv'))
for(i in seq_len(nrow(e))){z<-a[a$series_id==e$series_id[i]&a$scenario==e$scenario[i]&a$t>0,];if(e$comparison[i]!='current_all_years')z<-z[is.finite(z$revised),];z$y<-z[[if(e$comparison[i]=='revised_matched_years')'revised' else 'original']];f<-lm(y~t,z);stopifnot(nrow(z)==e$n_years[i],abs(coef(f)[2]-e$estimate[i])<1e-9,abs(sqrt(uncertainty_hac(f,z$year,3)[2,2])-e$se[i])<1e-9)}
v<-read.csv(paste0(o,'visit_values.csv'));s<-read.csv(paste0(o,'station_years.csv'))
key<-function(x)do.call(paste,c(x,sep='~'))
groups<-split(v,key(v[c('series_id','scenario','site','year')]))
for(i in seq_len(nrow(s))){z<-groups[[key(s[i,c('series_id','scenario','site','year')])]];den<-mean(tapply(z$density,z$zone,mean));bio<-mean(tapply(z$biomass,z$zone,mean));stopifnot(isTRUE(all.equal(den,s$density[i],tolerance=1e-9)),isTRUE(all.equal(bio,s$biomass[i],tolerance=1e-9)))}
cat('Verified',nrow(e),'model coefficients/SEs and',nrow(s),'station-year aggregations from visit values.\n')
