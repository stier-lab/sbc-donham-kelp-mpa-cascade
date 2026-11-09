# Read only: extract the previously reconstructed processor objects to private audit files.
out <- 'outputs/biomass_joint_checks_20261102/private'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
e <- new.env()
load('outputs/all_program_station_refits_20261102/stage/reconstruction.RData',envir=e)
for(n in c('kfm.sum','kfm.site','kfm.stipe.site','Swath.site.PISCO.sum.int','lter.sub')){
 stopifnot(is.data.frame(e[[n]]))
 write.csv(e[[n]],file.path(out,paste0(n,'.csv')),row.names=FALSE)
}
