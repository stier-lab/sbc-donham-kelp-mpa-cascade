# Optional raw KFM diagnostic. Run from repo root; does not modify primary fits.
local({
  suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
  source("code/R/nonfished_comparator_utils.R", local = TRUE)
  source("code/R/review_concern_utils.R", local = TRUE)
  out <- "tables"
  save_table <- function(x, suffix) fwrite(x, file.path(out, paste0("table_s_nonfished_", suffix, ".csv")))
  raw_path <- Sys.getenv("NONFISHED_KFM_COUNTS", "data/MBON/SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv")
  map_path <- Sys.getenv("NONFISHED_KFM_SITES", "data/MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv")
  if (!all(file.exists(c(raw_path, map_path)))) {
    save_table(data.table(status="skipped", reason="Optional raw KFM counts/site metadata unavailable; no current result was generated"), "status")
    if (Sys.getenv("NONFISHED_REQUIRED") == "1") stop("Required non-fished KFM inputs unavailable")
    message("Non-fished comparison skipped: optional KFM raw inputs unavailable")
  } else {
    protected <- c(list.files("data/harmonized", full.names=TRUE), "data/sumstats_final.csv",
                   "tables/table_02_meta_analysis.csv", "tables/table_03_cross_taxa_meta_regression.csv",
                   "tables/table_resistance_recovery.csv", "tables/table_s_historical_heatwave_summary.csv",
                   list.files("plots", "^fig_0[1-4].*\\.(png|pdf)$", full.names=TRUE),
                   "plots/fig_kelp_resilience.png", "plots/fig_kelp_resilience.pdf")
    protected <- protected[file.exists(protected)]
    protected_before <- tools::md5sum(protected)
    input_paths <- c(raw_path, map_path, "data/external/nonfished_comparators.csv",
                     "data/external/nonfished_legacy_pairs.csv", "data/harmonized/reference_sites.csv",
                     "data/harmonized/harmonized_response_ratios.csv",
                     "data/harmonized/harmonized_raw_responses.csv", "data/sumstats_final.csv",
                     "code/R/45_nonfished_comparators.R", "code/R/nonfished_comparator_utils.R",
                     "docs/NONFISHED_COMPARISON_PLAN.md")
    save_table(data.table(path=input_paths, md5=unname(tools::md5sum(input_paths))), "inputs")
    candidates <- fread("data/external/nonfished_comparators.csv")
    roster <- fread("data/external/nonfished_legacy_pairs.csv")
    meta <- fread(map_path)
    registry <- fread("data/harmonized/reference_sites.csv")
    stopifnot(!anyDuplicated(candidates$proj_taxon_id), !anyDuplicated(roster[,.(MPA, site_id)]))
    check_map <- merge(roster, unique(meta[data_source=="kfm", .(MPA=CA_MPA_Name_Short,site_id,status)]),
                       by=c("MPA","site_id","status"))
    stopifnot(nrow(check_map)==nrow(roster))
    for (m in unique(roster$MPA)) {
      actual <- sort(roster[MPA==m & status=="reference",site_id])
      expected <- sort(unique(sub("site|", "", registry[MPA==m & Source=="KFM" & included &
        Resp=="Den" & Taxa %in% c("M. franciscanus","S. purpuratus"),reference_id], fixed=TRUE)))
      stopifnot(identical(actual, expected))
    }
    focal <- data.table(proj_taxon_id=c("t-k-026","t-k-027"),
                        taxon_name=c("Mesocentrotus franciscanus","Strongylocentrotus purpuratus"),
                        sample_method="1mquad", nominal_area=2, role="focal")
    taxa <- rbindlist(list(candidates[,.(proj_taxon_id,taxon_name,sample_method,nominal_area,role)],focal))
    raw <- fread(raw_path, na.strings=c("NA",".",""))
    raw <- raw[data_source=="kfm" & site_id %in% roster$site_id & proj_taxon_id %in% taxa$proj_taxon_id]
    raw[,year:=as.integer(format(as.Date(date),"%Y"))]
    stopifnot(!anyNA(raw$year))
    raw <- raw[year>=1991 & year<=2022]
    n_before_taxon_join <- nrow(raw)
    raw <- merge(raw, taxa, by=c("proj_taxon_id","taxon_name","sample_method"))
    stopifnot(nrow(raw)==n_before_taxon_join)
    key <- c("site_id","date","subsite_id","transect_id","replicate_id","proj_taxon_id","sample_method")
    if (anyDuplicated(raw[,..key])) stop("Duplicate KFM measurement keys; no automatic deduplication")
    stopifnot(all(is.finite(raw$count) & raw$count>=0), all(is.finite(raw$area) & raw$area>0),
              !anyNA(raw[,..key]), setequal(unique(raw$proj_taxon_id),taxa$proj_taxon_id))
    audit <- raw[,.(rows=.N, zero_rows=sum(count==0), missing_counts=sum(!is.finite(count)),
                    area_values=paste(sort(unique(area)),collapse=";"),
                    sites=uniqueN(site_id), first_year=min(year),last_year=max(year)),
                 by=.(proj_taxon_id,taxon_name,sample_method,role)]
    save_table(audit,"raw_audit")
    raw[,density:=count/area]
    visit <- raw[,.(density=mean(density),n_replicates=.N),
                 by=.(proj_taxon_id,site_id,year,date,subsite_id,transect_id)]
    visit <- visit[,.(density=mean(density),n_replicates=sum(n_replicates)),
                   by=.(proj_taxon_id,site_id,year,date)]
    site <- visit[,.(density=mean(density),n_visits=.N,n_replicates=sum(n_replicates)),
                  by=.(proj_taxon_id,site_id,year)]
    stopifnot(!anyDuplicated(site[,.(proj_taxon_id,site_id,year)]))
    save_table(site,"site_annual")
    save_table(visit,"visits")
    annual_primary <- fread("data/harmonized/harmonized_response_ratios.csv")
    primary <- fread("data/sumstats_final.csv")
    primary_raw <- fread("data/harmonized/harmonized_raw_responses.csv")
    scope_defs <- list(
      kelp_matched=data.table(MPA=c("Gull Island SMR","Scorpion SMR","Santa Barbara Island SMR"),
                              y="Macrocystis pyrifera",taxon="M. pyrifera",resp="Bio"),
      urchin_matched=data.table(MPA=c("Harris Point SMR","Gull Island SMR","Scorpion SMR","Santa Barbara Island SMR"),
                                y=c(rep("Mesocentrotus franciscanus",3),"Strongylocentrotus purpuratus"),
                                taxon=c(rep("M. franciscanus",3),"S. purpuratus"),resp="Den"))
    all_annual <- all_metrics <- all_coverage <- all_rosters <- all_guards <- all_eligibility <- list()
    kelp_windows <- kelp_annual <- kelp_preyear_omissions <- list()
    n <- 0L; j <- 0L
    for (scope in names(scope_defs)) for (panel in c("core","extended")) {
      defs <- scope_defs[[scope]]
      panel_ids <- candidates[role=="core" | panel=="extended",proj_taxon_id]
      required_ids <- c(panel_ids,focal$proj_taxon_id)
      for (m in defs$MPA) {
        def <- defs[MPA==m]
        stopifnot(nrow(primary[MPA==m & Source=="KFM" & Taxa==def$taxon & Resp==def$resp])==1L)
        anchor <- annual_primary[source=="KFM" & CA_MPA_Name_Short==m & y==def$y & resp==def$resp & year>=1991 & year<=2022]
        stopifnot(!anyDuplicated(anchor$year), all(anchor$BA==ifelse(anchor$year<=2003,"Before","After")))
        sites_m <- as.data.frame(roster[MPA==m])
        per_taxon <- lapply(required_ids,function(id) nonfished_pair(as.data.frame(site),sites_m,id,anchor$year))
        names(per_taxon) <- required_ids
        years <- sort(Reduce(intersect,lapply(per_taxon,`[[`,"year")))
        j <- j+1L
        all_coverage[[j]] <- data.table(scope,panel,MPA=m,anchor_years=paste(anchor$year,collapse=";"),
          matched_years=paste(years,collapse=";"),n_anchor=nrow(anchor),n_matched=length(years),
          n_before=sum(years<=2003),n_after=sum(years>=2004),
          n_sites=nrow(sites_m),n_candidates=length(panel_ids),
          dropped_years=paste(setdiff(anchor$year,years),collapse=";"))
        for(first_after in c(2004L,2005L)) {
          n_pre <- sum(years<first_after); n_post <- sum(years>=first_after)
          all_eligibility[[length(all_eligibility)+1L]] <- data.table(scope,panel,MPA=m,first_after,
            n_before=n_pre,n_after=n_post,eligible=n_pre>=6L & n_post>=6L,
            reason=if(n_pre>=6L & n_post>=6L) "eligible" else "fewer than six matched years in a period")
        }
        if(!length(years)) next
        all_rosters[[j]] <- as.data.table(merge(data.frame(scope,panel,year=years),
                                               as.data.frame(roster[MPA==m]),by=NULL))
        for (id in required_ids) {
          d <- per_taxon[[id]]; d <- d[d$year %in% years,]
          if (!nrow(d)) next
          # Verify focal source measurement against the primary raw-density table.
          if (id %in% focal$proj_taxon_id) {
            sp <- focal[proj_taxon_id==id,taxon_name]
            target <- primary_raw[source=="KFM" & CA_MPA_Name_Short==m & taxon_name==sp & resp=="Den"]
            for (stat in c("mpa","reference")) {
              vals <- target[status==stat][match(d$year,year),value]
              reconstructed <- d[[if(stat=="mpa") "inside" else "reference"]]
              guard <- data.table(scope,panel,MPA=m,proj_taxon_id=id,year=d$year,status=stat,
                                  primary_density=vals,reconstructed_density=reconstructed,
                                  difference=reconstructed-vals)
              stopifnot(all(is.finite(vals)),max(abs(guard$difference))<1e-9)
              all_guards[[length(all_guards)+1L]] <- guard
            }
          }
          spec <- taxa[proj_taxon_id==id]
          for (offset_count in c(.5,.1,1)) {
            a <- offset_count/spec$nominal_area
            d$log_inside <- log(d$inside+a); d$log_reference <- log(d$reference+a)
            d$contrast <- d$log_inside-d$log_reference
            d$MPA <- m; d$proj_taxon_id <- id; d$taxon_name <- spec$taxon_name
            d$scope <- scope; d$panel <- panel; d$offset_count <- offset_count; d$offset_density <- a
            d$role <- spec$role
            n <- n+1L; all_annual[[n]] <- as.data.table(d)
            for(first_after in c(2004L,2005L)) {
              z <- nonfished_metrics(d,first_after)
              if(!is.null(z)) all_metrics[[length(all_metrics)+1L]] <- cbind(
                data.table(scope,panel,MPA=m,proj_taxon_id=id,taxon_name=spec$taxon_name,
                           role=spec$role,offset_count,first_after),as.data.table(z),
                data.table(fraction_both_zero=mean(d$inside==0 & d$reference==0),
                           fraction_both_zero_before=mean(d$inside[d$year<first_after]==0 &
                                                           d$reference[d$year<first_after]==0),
                           fraction_positive_inside=mean(d$inside>0),fraction_positive_reference=mean(d$reference>0)))
            }
          }
        }
        if (scope=="kelp_matched") {
          d <- data.frame(year=years,contrast=anchor$lnDiff[match(years,anchor$year)])
          if(panel=="core") {
            raw_kelp <- primary_raw[source=="KFM" & CA_MPA_Name_Short==m &
                                     taxon_name=="Macrocystis pyrifera" & resp=="Bio"]
            d$inside <- raw_kelp[status=="mpa"][match(years,year),value]
            d$reference <- raw_kelp[status=="reference"][match(years,year),value]
            d$raw_available <- is.finite(d$inside) & is.finite(d$reference)
            corrected <- anchor[match(years,year)]
            stopifnot(max(abs(log(corrected$mpa)-log(corrected$reference)-d$contrast))<1e-9)
            kelp_annual[[length(kelp_annual)+1L]] <- cbind(MPA=m,as.data.table(d))
            windows <- list(full_1991_2003=d[d$year<=2003,], recent_1996_2003=d[d$year>=1996 & d$year<=2003,],
                            exclude_2003=d[d$year<=2002,])
            for(w in names(windows)) {
              z <- windows[[w]]
              kelp_windows[[length(kelp_windows)+1L]] <- data.table(MPA=m,window=w,n_years=nrow(z),
                slope=if(nrow(z)>=6) unname(coef(lm(contrast~year,z))["year"]) else NA_real_,
                raw_available_years=sum(z$raw_available),
                both_zero_years=sum(z$inside==0 & z$reference==0,na.rm=TRUE),
                analysis_status="post-hoc prewindow and zero-support diagnostic; no annual-iid inference")
            }
            pre <- d[d$year<=2003,]
            for(omit in pre$year) {
              z <- pre[pre$year!=omit,]
              kelp_preyear_omissions[[length(kelp_preyear_omissions)+1L]] <- data.table(MPA=m,
                omitted_year=omit,n_years=nrow(z),slope=unname(coef(lm(contrast~year,z))["year"]),
                analysis_status="post-hoc leave-one-preyear-out; no annual-iid inference")
            }
          }
          for(first_after in c(2004L,2005L)) {
            z <- nonfished_metrics(d,first_after)
            if(!is.null(z)) all_metrics[[length(all_metrics)+1L]] <- cbind(
              data.table(scope,panel,MPA=m,proj_taxon_id="primary_kelp_lnRR",taxon_name="Macrocystis pyrifera",
                         role="focal_primary_scale",offset_count=NA_real_,first_after),as.data.table(z))
          }
        }
      }
    }
    annual <- rbindlist(all_annual,fill=TRUE); metrics <- rbindlist(all_metrics,fill=TRUE)
    coverage <- rbindlist(all_coverage); guards <- rbindlist(all_guards)
    stopifnot(nrow(annual)>0,nrow(metrics)>0)
    save_table(annual,"annual"); save_table(metrics,"pair_metrics"); save_table(coverage,"coverage")
    save_table(rbindlist(all_eligibility),"eligibility")
    save_table(rbindlist(kelp_windows),"kelp_prewindows")
    save_table(rbindlist(kelp_annual),"kelp_annual")
    save_table(rbindlist(kelp_preyear_omissions),"kelp_preyear_omissions")
    save_table(rbindlist(all_rosters),"roster"); save_table(guards,"primary_density_guard")
    endpoints <- c("before","pretrend","change")
    pooled <- omissions <- list()
    keys <- c("scope","panel","proj_taxon_id","taxon_name","role","offset_count","first_after")
    metric_groups <- metrics[,.(rows=list(.I)),by=keys]
    for(idx in metric_groups$rows) {
      z <- metrics[idx]
      if (!nrow(z)) next
      for(endpoint in endpoints) {
        ref_ids <- vapply(z$MPA,function(m) paste(sort(roster[MPA==m & status=="reference",site_id]),collapse=";"),character(1))
        groups <- review_reference_components(z$MPA,ref_ids)
        summary <- review_group_summary(z[[endpoint]],groups,TRUE)
        summary$p_signflip <- review_signflip(z[[endpoint]])
        pooled[[length(pooled)+1L]] <- cbind(z[1,..keys],endpoint,as.data.table(summary))
        for(m in z$MPA) omissions[[length(omissions)+1L]] <- cbind(z[1,..keys],endpoint,
          data.table(omitted_mpa=m,n_mpa=nrow(z)-1L,estimate=mean(z[[endpoint]][z$MPA!=m])))
      }
    }
    pooled <- rbindlist(pooled,fill=TRUE)
    pooled[,p_fdr:=NA_real_]
    pooled[,family_n:=0L]
    pooled[role %in% c("core","secondary"),c("p_fdr","family_n"):=
      list(p.adjust(p,method="BH"),sum(is.finite(p))),
      by=.(scope,panel,role,offset_count,first_after,endpoint)]
    save_table(pooled,"pooled");save_table(rbindlist(omissions),"omissions")
    report <- pooled[panel=="core" & offset_count==.5 & first_after==2004 & role=="core"]
    markdown_rows <- vapply(seq_len(nrow(report)),function(i) {
      z <- report[i]
      sprintf("| %s | %s | %s | %d | %.3f [%.3f, %.3f] | %.3f |",z$scope,z$taxon_name,z$endpoint,
              z$n_mpa,z$estimate,z$ci_low,z$ci_high,z$p_fdr)
    },character(1))
    writeLines(c("# Non-target KFM Reef Comparability Results", "",
      "Generated by module 45 from checksum-bound inputs. KFM-only descriptive diagnostic; no primary data, fits or figures changed.","",
      "The candidate set and rules were fixed in NONFISHED_COMPARISON_PLAN.md before response contrasts were calculated. Four core and three secondary sessile comparators were retained irrespective of effect estimates. These are not validated unaffected negative controls.","",
      "## Core results", "",
      paste0("Primary zero correction: 0.5 individuals per nominal 2 m2 quadrat or 60 m2 band; first after year 2004. Entries are log-offset density contrasts (pretrend per calendar year), not unshifted abundance ratios. Equal-MPA t intervals use delete-reference-group jackknife uncertainty. Core BH families contain four taxa per endpoint/scenario/scope. ",sum(report$ci_low<=0 & report$ci_high>=0)," of ",nrow(report)," core intervals span zero in this specification; neither null pretrends nor null changes establish equivalence."),"",
      "| Scope | Taxon | Endpoint | MPAs | Estimate [95% interval] | Adjusted P |",
      "|---|---|---|---:|---|---:|",markdown_rows,"",
      "## Interpretation", "",
      "Pre-designation levels and tracking vary among reef pairs. Sparse or all-zero candidate histories limit information; inspect table_s_nonfished_pair_metrics.csv and both raw-density and ratio plots. Opposing candidate differences cannot be averaged into a generic reef-quality score. The test does not determine the fraction of the kelp result due to siting.","",
      "The kelp-matched scope has three MPAs and the urchin-matched scope four. Fixed legacy reference rosters reproduce primary raw urchin densities on matched years; kelp uses its current primary lnRR without reconstructing ambiguous kelp taxon IDs. All candidate records share the full site roster and each scope's matched calendar. New post-designation stations, satellite-only MPAs and mainland programmes are not covered.","",
      "The current primary kelp lnRRs have positive point slopes over the full 1991-2003 matched pre-period at all three kelp MPAs. This is not a claim of consistent linear growth: a post-hoc shorter 1996-2003 check is flat at Scorpion, whose available raw biomass entries are jointly zero in 1997-2003. Gull and Santa Barbara Island retain positive recent-window slopes. The raw-response export lacks some earlier kelp entries even though their lnRRs exist; unavailable raw values are not replaced by zero. The lnRR identity is verified against the stored corrected inside/reference indices, not an independent reconstruction of all kelp biomass. table_s_nonfished_kelp_prewindows.csv and table_s_nonfished_kelp_preyear_omissions.csv report all window/omission checks; table_s_nonfished_kelp_annual.csv retains observations and raw-availability flags. A later before/after increase cannot by itself establish a designation effect or parallel untreated trajectories.","",
      "Annual samples do not increase independent MPA replication. The group-aware t intervals are approximate with only three/four groups and shared environmental variation. Exact sign-flip P values cannot reach 0.05 with these sample sizes and assume symmetric/exchangeable MPA contrasts. Pretrend slopes and tracking are descriptive; annual-iid P values are not used.","",
      "## Reproduce", "",
      "Run `Rscript code/R/45_nonfished_comparators.R` and `Rscript code/tests/test_nonfished_comparators.R` from the repository root. Module 45 is optional in run_all.R; set NONFISHED_KFM_COUNTS and NONFISHED_KFM_SITES for alternative raw-input locations, and NONFISHED_REQUIRED=1 to reject absent inputs. The default raw files are local optional inputs, not newly redistributed raw records.","",
      "Offsets 0.1 and 1, a 2005 boundary, extended taxa, omission estimates, row-level effort audits and primary-file invariance checks are exported under tables/table_s_nonfished_*.csv. See the plan for the NPS protocol source and taxon-name limitations."),
      "docs/NONFISHED_COMPARISON_RESULTS.md")
    stopifnot(!anyDuplicated(annual[,.(scope,panel,MPA,proj_taxon_id,offset_count,year)]))
    # No titles/subtitles: the report supplies the biological question and legends.
    display <- annual[scope=="urchin_matched" & panel=="core" & offset_count==.5 & role=="core"]
    display[,taxon_label:=sub(" ","\n",taxon_name)]
    p <- ggplot(display,aes(year,contrast,colour=MPA)) +
      geom_hline(yintercept=0,colour="grey70",linewidth=.3) +
      geom_vline(xintercept=2003.5,colour="grey50",linetype="dashed",linewidth=.3) +
      geom_point(size=.7,alpha=.65) + facet_wrap(~taxon_label,ncol=2) +
      scale_colour_manual(values=c("Gull Island SMR"="#0072B2","Harris Point SMR"="#D55E00",
                                   "Santa Barbara Island SMR"="#009E73","Scorpion SMR"="#CC79A7")) +
      labs(x="Year",y="Log density contrast (inside / reference)",colour=NULL) +
      theme_classic(base_size=10)+theme(legend.position="bottom",legend.text=element_text(size=8))
    ggsave("plots/fig_s_nonfished_histories.png",p,width=7,height=4.8,dpi=300,bg="white")
    ggsave("plots/fig_s_nonfished_histories.pdf",p,width=7,height=4.8)
    # Show actual pre/post densities as well as ratios, preserving sampling gaps.
    long <- melt(display,id.vars=c("MPA","taxon_label","year"),measure.vars=c("inside","reference"),
                 variable.name="status",value.name="density")
    p2 <- ggplot(long,aes(year,density,colour=status)) + geom_point(size=.6,alpha=.75) +
      geom_vline(xintercept=2003.5,colour="grey65",linetype="dashed",linewidth=.25) +
      facet_grid(taxon_label~MPA,scales="free_y") +
      scale_colour_manual(values=c(inside="#0072B2",reference="#D55E00")) +
      labs(x="Year",y=expression("Density (individuals "*m^{-2}*")"),colour=NULL) +
      theme_classic(base_size=9)+theme(legend.position="bottom",strip.text=element_text(size=7))
    ggsave("plots/fig_s_nonfished_densities.png",p2,width=9,height=6,dpi=300,bg="white")
    ggsave("plots/fig_s_nonfished_densities.pdf",p2,width=9,height=6)
    protected_after <- tools::md5sum(protected)
    stopifnot(identical(protected_before,protected_after))
    save_table(data.table(path=protected,before=unname(protected_before),after=unname(protected_after),
                          unchanged=unname(protected_before)==unname(protected_after)),"protected")
    save_table(data.table(status="complete",scope="KFM legacy reefs only; descriptive comparability diagnostic",
                         n_raw=nrow(raw),n_configured_pairs=uniqueN(roster$MPA),
                         n_contributing_pairs=uniqueN(metrics$MPA),
                         n_kelp_mpas=uniqueN(metrics[scope=="kelp_matched",MPA]),
                         n_urchin_mpas=uniqueN(metrics[scope=="urchin_matched",MPA]),n_core=sum(candidates$role=="core"),
                         n_secondary=sum(candidates$role=="secondary"),n_guard_rows=nrow(guards),
                         max_density_reconstruction_error=max(abs(guards$difference)),
                         eligible_for_primary=FALSE),"status")
    save_table(data.table(package=c("R","data.table","ggplot2"),version=c(as.character(getRversion()),
                  as.character(packageVersion("data.table")),as.character(packageVersion("ggplot2")))),"versions")
    message("Non-fished KFM comparator analysis complete; primary files unchanged")
  }
})
