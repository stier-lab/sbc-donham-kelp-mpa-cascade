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
  missing_path <- file.path(root,"data/harmonized/lobster_missing_lengths.json")
  if(file.exists(missing_path)) {
    m <- jsonlite::fromJSON(missing_path,simplifyVector=FALSE)
    stopifnot(identical(m$release,"lobster_missing_lengths_2026_10_30"),isTRUE(m$applied),
      identical(m$scope$programme,"LTER"),identical(m$scope$taxon,"Panulirus interruptus"),
      identical(m$scope$response,"Bio"),identical(m$files[[name]]$before_md5,result))
    result <- m$files[[name]]$after_md5
  }
  pisco_path <- file.path(root,"data/harmonized/pisco_baseline_policy.json")
  if(file.exists(pisco_path)) {
    p <- jsonlite::fromJSON(pisco_path,simplifyVector=FALSE)
    stopifnot(identical(p$release,"pisco_fixed_stations_strict3_2026_10_31"),isTRUE(p$applied),
      p$minimum_before_years==3,identical(p$files[[name]]$before_md5,result))
    result <- p$files[[name]]$after_md5
  }
  result
}

# Changed membership is explicit in the latest release; older row counts stay historical.
current_release_rows <- function(root,name,previous_rows) {
 p<-file.path(root,"data/harmonized/pisco_baseline_policy.json")
 if(!file.exists(p)) return(previous_rows)
 z<-jsonlite::fromJSON(p,simplifyVector=FALSE)
 if(is.null(z$files[[name]])) return(previous_rows)
 z$files[[name]]$rows
}
