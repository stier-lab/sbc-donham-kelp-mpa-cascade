# Diagnostic reconstruction of Tom Bell's delivered arrays; no primary writes.
local({
 library(ncdf4);library(readxl)
 base <- "data/external/tom_bell_2026-10-20"
 ncpath <- "local_archive/metadata_audit_2026-10-16/CAkelpCanopyEnv_2021.nc"
 d <- nc_open(ncpath);on.exit(nc_close(d))
 lat <- ncvar_get(d,"lat");lon <- ncvar_get(d,"lon");yr <- ncvar_get(d,"year")
 key <- function(a,b)paste(round(a,6),round(b,6),sep="|")
 keys <- key(lat,lon);stopifnot(!anyDuplicated(keys))
 w <- read_excel(file.path(base,"03_MPA_Ref_Locations.xlsx"),sheet="Pixel Locations",col_names=FALSE,.name_repair="minimal")
 sets <- lapply(seq(2,ncol(w),2),function(j){
  a<-suppressWarnings(as.numeric(w[[j]][-(1:2)]));b<-suppressWarnings(as.numeric(w[[j+1]][-(1:2)]));ok<-is.finite(a)&is.finite(b)
  idx<-match(key(a[ok],b[ok]),keys);stopifnot(!anyNA(idx),!anyDuplicated(idx))
  list(site=as.integer(w[[j]][1]),group=as.character(w[[j]][2]),idx=idx)
 })
 q <- matrix(NA_real_,length(sets),length(yr));missing_fraction<-q
 for(t in seq_along(yr)){
  v<-ncvar_get(d,"biomass",start=c(1,t),count=c(-1,1))
  for(j in seq_along(sets)){
   x<-v[sets[[j]]$idx]
   # Empty MATLAB nansum input gives zero; missing-fraction 0/0 is NaN.
   q[j,t]<-sum(x,na.rm=TRUE)
   if(length(x)){
    missing_fraction[j,t]<-mean(is.na(x))
    if(missing_fraction[j,t]>.25)q[j,t]<-NA_real_
   }
  }
 }
 years<-1984:2021
 annual<-sapply(years,function(y)apply(q[,yr==y,drop=FALSE],1,function(x){if(sum(is.na(x))>1)NA_real_ else max(x,na.rm=TRUE)}))
 csvs<-lapply(c("MPAs","MPAC1","MPAC2"),function(s)read.csv(file.path(base,paste0("04_MPA_Runs_new_",s,".csv")),check.names=FALSE))
 names(csvs)<-c("MPA","Ref_1","Ref_2")
 raw<-read.csv("data/LANDSAT/MPA_Runs_new.csv",check.names=FALSE)
 sites<-read.csv(file.path(base,"03_Site_List.csv"),check.names=FALSE)
 rec<-list();matching<-list()
 for(j in seq_along(sets)){
  s<-sets[[j]];tab<-csvs[[s$group]];a<-annual[j,];b<-as.numeric(tab[tab$Site==s$site,-1])
  stopifnot(length(b)==38)
  eq<-is.na(a)==is.na(b);finite<-is.finite(a)&is.finite(b)
  rec[[j]]<-data.frame(site=s$site,group=s$group,n_pixels=length(s$idx),missing_pattern_equal=all(eq),paired_finite=sum(finite),max_abs_difference=if(any(finite))max(abs(a[finite]-b[finite]))else NA_real_,n_different=sum(!eq|ifelse(finite,abs(a-b)>1e-6,FALSE)))
  # Align using delivered centroid, not ambiguous MPA name spelling.
  m<-which(raw$status=="mpa"&abs(raw$lat-sites$Lat[s$site])<1e-8&abs(raw$lon-sites$Lon[s$site])<1e-8)
  if(length(m)==1){
   cand<-which(raw$MPA==raw$MPA[m] & raw$status==if(s$group=="MPA")"mpa" else "reference")
   for(i in cand){old<-as.numeric(raw[i,as.character(years)]);matching[[length(matching)+1]]<-data.frame(site=s$site,MPA=raw$MPA[m],group=s$group,rep=raw$rep[i],exact=all(is.na(old)==is.na(b))&&all(abs(old[is.finite(old)&is.finite(b)]-b[is.finite(old)&is.finite(b)])<1e-6))}
  }
 }
 res<-do.call(rbind,rec);matches<-do.call(rbind,matching)
 write.csv(res,file.path(base,"reconstruction_check.csv"),row.names=FALSE)
 write.csv(matches,file.path(base,"current_input_matches.csv"),row.names=FALSE)
 cross<-list()
 for(j in seq_along(sets)){
  s<-sets[[j]]
  for(g in names(csvs)){
   b<-as.numeric(csvs[[g]][s$site,-1]);a<-annual[j,];finite<-is.finite(a)&is.finite(b)
   cross[[length(cross)+1]]<-data.frame(site=s$site,mask_group=s$group,delivered_group=g,n_pixels=length(s$idx),exact=all(is.na(a)==is.na(b))&&all(abs(a[finite]-b[finite])<1e-6))
  }
 }
 write.csv(do.call(rbind,cross),file.path(base,"mask_array_crossmatches.csv"),row.names=FALSE)
 # Pre-2012 correlations are diagnostics, not a change to averaged primary references.
 cors<-lapply(1:24,function(i){inside<-as.numeric(csvs$MPA[i,-1]);r1<-as.numeric(csvs$Ref_1[i,-1]);r2<-as.numeric(csvs$Ref_2[i,-1]);pre<-years<2012
  safe<-function(x,y){ok<-is.finite(x)&is.finite(y);if(sum(ok)>2&&sd(x[ok])>0&&sd(y[ok])>0)cor(x[ok],y[ok]) else NA_real_}
  a<-safe(inside[pre],r1[pre]);b<-safe(inside[pre],r2[pre]);data.frame(site=i,annual_pre2012_r1=a,annual_pre2012_r2=b,best=if(all(is.na(c(a,b))))NA_integer_ else which.max(replace(c(a,b),is.na(c(a,b)),-Inf)))})
 write.csv(do.call(rbind,cors),file.path(base,"pre2012_reference_correlations.csv"),row.names=FALSE)
 cat("Reconstructed groups:",nrow(res),"; exact groups:",sum(res$n_different==0),"; discrepant groups:\n");print(res[res$n_different>0,])
 cat("Current raw rows:",nrow(raw),"; exact matched rows:",length(unique(paste(matches$MPA[matches$exact],matches$group[matches$exact]=="MPA",matches$rep[matches$exact]))),"\n")
})
