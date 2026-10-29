suppressPackageStartupMessages(library(data.table))
source("code/R/nonfished_comparator_utils.R")
source("code/R/review_concern_utils.R")

# A known baseline, pretrend and before/after step are not conflated.
d <- data.frame(year=1991:2010)
d$contrast <- 2 + .1*(d$year-1991) + ifelse(d$year>=2004,1,0)
z <- nonfished_metrics(d,2004)
stopifnot(abs(z$pretrend-.1)<1e-12,
          abs(z$before-mean(d$contrast[d$year<2004]))<1e-12,
          abs(z$change-(mean(d$contrast[d$year>=2004])-mean(d$contrast[d$year<2004])))<1e-12,
          is.null(nonfished_metrics(d[d$year<2007,],2004)))

# Missing one reference is not an observed zero or a complete fixed pair.
r <- data.frame(site_id=c("I","C1","C2"),status=c("mpa","reference","reference"))
s <- expand.grid(site_id=r$site_id,year=2000:2002,stringsAsFactors=FALSE)
s$proj_taxon_id <- "X";s$density <- 2
s$density[s$site_id=="I"] <- 0
s <- s[!(s$site_id=="C2" & s$year==2001),]
p <- nonfished_pair(s,r,"X",2000:2002)
stopifnot(identical(p$year,c(2000L,2002L)),all(p$inside==0),all(p$reference==2))
stopifnot(nrow(nonfished_pair(s,r,"X",2015:2020))==0L)

# Three/four MPA sign flips cannot yield a conventional 0.05 test.
stopifnot(review_signflip(c(1,2,3))==.25,review_signflip(c(1,2,3,4))==.125)
summary <- review_group_summary(c(1,2,3,4),c("A","A","B","C"),TRUE)
stopifnot(summary$n_reference_groups==3,summary$reference_group_df==2)

prefix <- "tables/table_s_nonfished_"
read <- function(s) fread(paste0(prefix,s,".csv"))
status <- read("status");stopifnot(status$status=="complete")
guard <- read("primary_density_guard")
stopifnot(nrow(guard)>0,all(abs(guard$difference)<1e-9))
protected <- read("protected");stopifnot(all(protected$unchanged))
stopifnot(identical(unname(tools::md5sum(protected$path)),protected$after))
inputs <- read("inputs");stopifnot(identical(unname(tools::md5sum(inputs$path)),inputs$md5))
annual <- read("annual");metrics <- read("pair_metrics");pooled <- read("pooled")
eligible <- read("eligibility")
stopifnot(all(eligible$eligible),all(eligible$n_before>=6),all(eligible$n_after>=6))
kelp_pre <- read("kelp_prewindows")
stopifnot(nrow(kelp_pre)==9,all(kelp_pre$raw_available_years<=kelp_pre$n_years),
          all(kelp_pre$both_zero_years<=kelp_pre$raw_available_years))
stopifnot(!anyDuplicated(annual[,.(scope,panel,MPA,proj_taxon_id,offset_count,year)]),
          all(metrics$n_before>=6),all(metrics$n_after>=6),
          all(pooled[role=="core",family_n]==4),all(pooled[role=="secondary",family_n]==3),
          nrow(pooled[role=="focal_primary_scale"])>0)
# Every candidate and focal count uses an identical calendar within each panel.
cal <- annual[,.(years=paste(year,collapse=";")),by=.(scope,panel,MPA,proj_taxon_id,offset_count)]
stopifnot(all(cal[,.(n=uniqueN(years)),by=.(scope,panel,MPA)]$n==1))
for (i in seq_len(nrow(pooled))) {
  row <- pooled[i]
  pair <- metrics[scope==row$scope & panel==row$panel & proj_taxon_id==row$proj_taxon_id &
                    first_after==row$first_after &
                    ((is.na(offset_count) & is.na(row$offset_count)) | offset_count==row$offset_count)]
  stopifnot(nrow(pair)==row$n_mpa,abs(mean(pair[[row$endpoint]])-row$estimate)<1e-12,
            row$reference_group_df==row$n_mpa-1)
}
cat("Non-fished tests passed: estimands, missing-reference exclusion, matching, pooling, multiplicity, source densities and primary invariance\n")
