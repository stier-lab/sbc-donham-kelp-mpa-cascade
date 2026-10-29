# Follow the verified scoped correction without rewriting earlier release history.
lobster_release_hash <- function(root,name,parent_md5) {
  path <- file.path(root,"data/harmonized/lobster_bootstrap_correction.json")
  if (!file.exists(path)) return(parent_md5)
  z <- jsonlite::fromJSON(path,simplifyVector=FALSE)
  stopifnot(identical(z$release,"lobster_singleton_sampling_2026_10_18"),isTRUE(z$applied),
    identical(z$scope$programme,"PISCO"),identical(z$scope$response,"Bio"),
    identical(z$scope$taxon,"Panulirus interruptus"),
    identical(z$files[[name]]$before_md5,parent_md5),z$a==0.001352821,z$b==2.913963113)
  result <- z$files[[name]]$after_md5
  next_path <- file.path(root,"data/harmonized/urchin_bootstrap_correction.json")
  if (file.exists(next_path)) {
    u <- jsonlite::fromJSON(next_path,simplifyVector=FALSE)
    stopifnot(identical(u$release,"urchin_singleton_sampling_2026_10_18"),isTRUE(u$applied),
      identical(u$files[[name]]$before_md5,result),identical(u$scope$response,"Bio"),
      setequal(unlist(u$scope$programmes),c("PISCO","LTER")),
      setequal(unlist(u$scope$taxa),c("Mesocentrotus franciscanus","Strongylocentrotus purpuratus")))
    result <- u$files[[name]]$after_md5
  }
  fixed_path <- file.path(root,"data/harmonized/kelp_fixed_calibration.json")
  if (file.exists(fixed_path)) {
    k <- jsonlite::fromJSON(fixed_path,simplifyVector=FALSE)
    calibration <- file.path(root,"data/external/sbc1275/Macrocystis_Biomass_Relationships_20261016.csv")
    stopifnot(identical(k$release,"kelp_fixed_sbc1275_2026_10_29"),isTRUE(k$applied),
      identical(k$package,"knb-lter-sbc.127.5"),k$model==1,
      identical(k$files[[name]]$before_md5,result),
      identical(k$calibration_md5,"8ddaf77458d21adfb6d0398cc1c51fda"),
      unname(tools::md5sum(calibration))==k$calibration_md5)
    result <- k$files[[name]]$after_md5
  }
  result
}
