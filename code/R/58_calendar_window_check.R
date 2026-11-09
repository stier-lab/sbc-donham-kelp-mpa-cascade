# Verify that all primary trend contrasts use surveyed calendar endpoints.
local({
  root <- here::here()
  source(file.path(root,'code/R/primary_uncertainty_utils.R'),local=TRUE)
  paths <- c('data/sumstats_final.csv','data/harmonized/harmonized_response_ratios.csv',
    'data/harmonized/harmonized_landsat_rr.csv','code/R/08_effect_sizes.R',
    'code/R/primary_uncertainty_utils.R','code/R/58_calendar_window_check.R')
  hashes <- tools::md5sum(file.path(root,paths))
  primary <- read.csv(file.path(root,paths[1])); primary$key <- uncertainty_key(primary)
  annual <- uncertainty_annual(root,primary); defs <- uncertainty_definitions(root)
  windows <- do.call(rbind,lapply(seq_len(nrow(primary)),function(i){
    row <- primary[i,]; r <- uncertainty_reconstruct(row,annual,defs); d <- r$data
    post <- d[d$BA=='After',]; first <- min(post$year)
    last <- max(post$year[post$year<=first+11])
    trend <- row$Model %in% c('Linear','Sigmoid','Asymptotic')
    stopifnot(first %in% d$year,last %in% d$year, last>first,last-first<=11,
      all(d$time.model==pmax(d$year-row$MPA_Start,0)),
      all(r$ends==c(first,last)-row$MPA_Start))
    data.frame(row[c('key','MPA','Taxa','Source','Resp','Model','AnalysisType','N')],
      quantity=if(trend)'Fitted change between surveyed post-designation years' else
        if(row$Model=='Step')'After-minus-before mean' else 'Mean over retained survey years',
      contrast_start_year=if(trend)first else NA_real_,
      contrast_end_year=if(trend)last else NA_real_,
      elapsed_years=if(trend)last-first else NA_real_,
      extrapolated=FALSE, clock=if(trend)'calendar years since designation' else 'not applicable',
      observed_years=paste(d$year,collapse=';'),
      reconstructed_estimate=r$estimate,reconstructed_se=r$iid_se)
  }))
  stopifnot(nrow(windows)==143,!anyDuplicated(windows$key),!any(windows$extrapolated),
    identical(hashes,tools::md5sum(file.path(root,paths))))
  write.csv(windows,file.path(root,'tables/table_s_primary_calendar_windows.csv'),row.names=FALSE)
  write.csv(data.frame(path=paths,md5=unname(hashes)),
    file.path(root,'tables/table_s_primary_calendar_inputs.csv'),row.names=FALSE)
  cat('PASS: all143 primary fits reconstructed; all trend endpoints are surveyed calendar years.\n')
  print(table(windows$elapsed_years,useNA='ifany'))
})
