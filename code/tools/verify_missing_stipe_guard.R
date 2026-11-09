suppressPackageStartupMessages(library(dplyr))
D<-'/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026'
# Exercise exact source expressions against the real reconstructed station table.
expr<-as.list(parse(file.path(D,'code/R/04_pisco_processing.R')))
assignment<-function(n)Filter(function(x)is.call(x)&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name(n)),expr)
Swath.edit<-read.csv(file.path(D,'data/PISCO/MLPA_kelpforest_swath_2024.csv'));SURVEY_MONTH_START<-5;SURVEY_MONTH_END<-10
eval(assignment('pisco_missing_stipe_sites')[[1]])
stopifnot(any(pisco_missing_stipe_sites$site=='UNDERWATER_ARCH' & pisco_missing_stipe_sites$year==2023))
All_PISCO.short.sub<-read.csv('outputs/all_program_station_refits_20261102/pisco_sites.csv');before<-All_PISCO.short.sub
start<-which(vapply(expr,function(x)is.call(x)&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name('.pisco_unknown_stipes')),logical(1)))
end<-which(vapply(expr,function(x)is.call(x)&&grepl('All_PISCO.mpa\\$biomass',paste(deparse(x[[2]]),collapse=''))&&grepl('pisco_role_key',paste(deparse(x),collapse='')),logical(1)))
stopifnot(length(start)==1,length(end)==1)
for(k in seq.int(start,end))eval(expr[[k]])
stopifnot(all(is.na(All_PISCO.short.sub$biomass[.pisco_unknown_stipes])))
stopifnot(identical(before$density,All_PISCO.short.sub$density),identical(before$count,All_PISCO.short.sub$count))
z<-All_PISCO.mpa[All_PISCO.mpa$CA_MPA_Name_Short=='Point Vicente SMCA' & All_PISCO.mpa$year==2023 & All_PISCO.mpa$y=='MACPYRAD' & All_PISCO.mpa$site_status=='reference',]
stopifnot(nrow(z)==1,is.na(z$biomass))
# Unaffected biomass stays identical; observed zero stipes are not flagged.
stopifnot(identical(before$biomass[!.pisco_unknown_stipes],All_PISCO.short.sub$biomass[!.pisco_unknown_stipes]))
cat('Production missing-stipe safeguard verified on actual source and station records; unknown reference biomass remains NA; other biomass/count fields unchanged.\n')
