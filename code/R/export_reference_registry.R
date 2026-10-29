# One-time reference identity/footprint audit. Never sourced by routine analyses.
# Rscript code/R/export_reference_registry.R /path/to/data-processing-repo

local({
  suppressPackageStartupMessages({ library(here); library(dplyr) })
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 1L) stop("Supply the source data-processing repository path.")
  data_repo <- normalizePath(args[1], mustWork = TRUE)
  source(here::here("code", "R", "00c_analysis_constants.R"), local = TRUE)
  source(here::here("code", "R", "review_concern_utils.R"), local = TRUE)
  source_relative <- c(
    pisco_mapping = "data/PISCO/master_site_table_Emilyedit.csv",
    mbon_mapping = "data/MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv",
    pisco_swath_footprint = "data/PISCO/MLPA_kelpforest_swath_2024.csv",
    pisco_fish_footprint = "data/PISCO/MLPA_kelpforest_fish_2024.csv",
    exclusion_definitions = "code/R/01_utils.R",
    pisco_processing_definitions = "code/R/04_pisco_processing.R",
    kfm_processing_definitions = "code/R/05_kfm_processing.R",
    lter_processing_definitions = "code/R/06_lter_processing.R")
  source_files <- setNames(file.path(data_repo, source_relative), names(source_relative))
  if (any(!file.exists(source_files))) stop("Missing source audit files: ",
    paste(source_relative[!file.exists(source_files)], collapse = ", "))
  source_files <- setNames(normalizePath(source_files, mustWork = TRUE), names(source_files))
  analysis_relative <- review_reference_analysis_files()
  analysis_files <- setNames(file.path(here::here(), analysis_relative), names(analysis_relative))
  source_md5 <- tools::md5sum(source_files); analysis_md5 <- tools::md5sum(analysis_files)
  primary <- read.csv(analysis_files[["primary_results"]])
  input <- review_design_input(bind_rows(read.csv(analysis_files[["harmonized_response_ratios"]]),
    read.csv(analysis_files[["harmonized_landsat_rr"]])), primary, RESILIENCE_TAXA_SHORT)
  keys <- review_reference_key_manifest(input)
  assignment <- function(name) {
    expressions <- parse(source_files[["exclusion_definitions"]])
    selected <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
      identical(x[[2]], as.name(name)), as.list(expressions))
    if (length(selected) != 1L) stop("Missing or ambiguous source exclusion list: ", name)
    value <- eval(selected[[1]][[3]], envir = baseenv())
    if (!is.character(value) || anyNA(value)) stop("Invalid source exclusion list: ", name)
    value
  }
  pisco_exclusions <- assignment("EXCLUDED_REFERENCE_SITES")
  kfm_exclusions <- assignment("EXCLUDED_KFM_SITES")
  # Both MBON IDs and standalone LTER codes occur in the source mapping.
  # Apply the same exclusions as the count, kelp and fish processing branches.
  lter_exclusions <- union(assignment("EXCLUDED_LTER_SITES"),
                          assignment("EXCLUDED_LTER_SITE_CODES"))
  pisco_map <- read.csv(source_files[["pisco_mapping"]]) %>%
    filter(site_status == "reference", BaselineRegion == "SOUTH", nzchar(site)) %>%
    transmute(Source = "PISCO", MPA = CA_MPA_Name_Short, site_id = as.character(site),
      reference_id = paste0("site|", site), reference_name = as.character(site_new),
      mapping_file = source_files[["pisco_mapping"]], identity_verified = TRUE) %>% distinct()
  mbon_map <- read.csv(source_files[["mbon_mapping"]]) %>% filter(status == "reference", nzchar(site_id)) %>%
    transmute(Source = toupper(data_source), MPA = CA_MPA_Name_Short, site_id = as.character(site_id),
      reference_id = paste0("site|", site_id), reference_name = as.character(site_id),
      mapping_file = source_files[["mbon_mapping"]], identity_verified = TRUE) %>% distinct()
  landsat_map <- input %>% filter(source == "Landsat") %>% distinct(CA_MPA_Name_Short) %>%
    transmute(Source = "Landsat", MPA = CA_MPA_Name_Short, site_id = NA_character_,
      reference_id = paste0("polygon-label|", MPA), reference_name = MPA,
      mapping_file = analysis_relative[["harmonized_landsat_rr"]], identity_verified = FALSE)
  registry <- bind_rows(pisco_map, mbon_map, landsat_map)
  if (!requireNamespace("data.table", quietly = TRUE)) stop("The one-time exporter requires data.table.")
  swath <- as.data.frame(data.table::fread(source_files[["pisco_swath_footprint"]],
    select = c("site", "year", "month", "campus", "zone"))) %>%
    filter(campus == "VRG" | (campus == "UCSB" & !zone %in% c("DEEP", "SPECIALO"))) %>%
    filter(month >= 5L, month <= 10L) %>% distinct(site, year) %>% mutate(survey = "swath")
  fish <- as.data.frame(data.table::fread(source_files[["pisco_fish_footprint"]],
    select = c("site", "year", "campus", "level"))) %>%
    filter(campus %in% c("UCSB", "VRG"), level != "CAN") %>% distinct(site, year) %>% mutate(survey = "fish")
  footprint <- bind_rows(swath, fish)
  series_years <- split(input$year, input$series_id)
  reference_sites <- keys %>% select(series_id, MPA, Taxa, Source, Resp) %>%
    left_join(registry, by = c("MPA", "Source"), relationship = "many-to-many") %>%
    mutate(reference_id = ifelse(is.na(reference_id), paste0("unmapped|", series_id), reference_id),
      identity_verified = coalesce(identity_verified, FALSE),
      excluded = (Source == "PISCO" & site_id %in% pisco_exclusions) |
        (Source == "KFM" & Taxa != "S. pulcher" & site_id %in% kfm_exclusions) |
        (Source == "LTER" & site_id %in% lter_exclusions),
      exclusion_reason = ifelse(excluded, "upstream reference-site exclusion", ""))
  reference_sites$reference_years <- vapply(seq_len(nrow(reference_sites)), function(i) {
    z <- reference_sites[i, ]; years <- series_years[[z$series_id]]
    if (grepl("^unmapped\\|", z$reference_id)) years <- integer()
    if (z$Source == "PISCO") years <- intersect(years, footprint$year[
      footprint$site == z$site_id & footprint$survey == if (z$Taxa == "S. pulcher") "fish" else "swath"])
    paste(sort(unique(years)), collapse = ",")
  }, character(1))
  file_md5 <- setNames(c(unname(source_md5), unname(analysis_md5)), c(source_files, analysis_relative))
  reference_sites <- reference_sites %>% mutate(
    registry_version = 1L,
    n_reference_years = ifelse(nzchar(reference_years), lengths(strsplit(reference_years, ",", fixed = TRUE)), 0L),
    included = !excluded & nzchar(reference_years),
    coverage_kind = ifelse(Source == "PISCO", "observed survey site/year footprint", "declared membership; harmonized series coverage"),
    contribution_audit = case_when(grepl("^unmapped\\|", reference_id) ~ "No declared reference identity available",
      Source == "PISCO" ~ "mapped reference with survey footprint; biomass missingness not reconstructed",
      Source == "Landsat" ~ "MPA pair label only; shared polygon identity/overlap unverified",
      TRUE ~ "declared MBON reference; individual site contributions not reconstructed"),
    mapping_md5 = unname(file_md5[mapping_file]),
    footprint_file = ifelse(Source == "PISCO", ifelse(Taxa == "S. pulcher",
      source_files[["pisco_fish_footprint"]], source_files[["pisco_swath_footprint"]]), NA_character_),
    footprint_md5 = unname(file_md5[footprint_file]),
    exclusion_definitions_md5 = unname(source_md5[match(source_files[["exclusion_definitions"]], names(source_md5))])) %>%
    arrange(series_id, reference_id)
  if (anyDuplicated(reference_sites[c("series_id", "reference_id")]))
    stop("Ambiguous repeated reference identity in source mappings; resolve explicitly before exporting.")
  if (!identical(unname(source_md5), unname(tools::md5sum(source_files))) ||
      !identical(unname(analysis_md5), unname(tools::md5sum(analysis_files))))
    stop("An audit/input file changed during export; rerun from a stable snapshot.")
  derived_relative <- c(reference_sites = "data/harmonized/reference_sites.csv",
                       reference_primary_keys = "data/harmonized/reference_primary_keys.csv")
  write.csv(reference_sites, here::here(derived_relative[["reference_sites"]]), row.names = FALSE)
  write.csv(keys, here::here(derived_relative[["reference_primary_keys"]]), row.names = FALSE)
  mapping_provenance <- bind_rows(
    data.frame(role = names(analysis_relative), kind = "analysis_input", path = unname(analysis_relative), md5 = unname(analysis_md5)),
    data.frame(role = names(source_files), kind = "source_audit", path = unname(source_files), md5 = unname(source_md5)),
    data.frame(role = names(derived_relative), kind = "derived_registry", path = unname(derived_relative),
      md5 = unname(tools::md5sum(here::here(derived_relative)))),
    data.frame(role = "registry_exporter", kind = "generator", path = "code/R/export_reference_registry.R",
      md5 = unname(tools::md5sum(here::here("code/R/export_reference_registry.R")))))
  write.csv(mapping_provenance, here::here("data/harmonized/mapping_provenance.csv"), row.names = FALSE)
  csv_inputs <- sort(list.files(here::here("data", "harmonized"), pattern = "\\.csv$", full.names = TRUE))
  writeLines(paste(unname(tools::md5sum(csv_inputs)), basename(csv_inputs)),
             here::here("data", "harmonized", "checksums.md5"))
  checked <- review_read_reference_registry(input)
  stopifnot(nrow(checked$reference_sites) == nrow(reference_sites))
  writeLines(c("# Reference identity registry", "",
    "These small, analysis-ready derivatives are the only reference-registry inputs to routine script 27. The ordinary pipeline never reads the original PISCO fish/swath files, site mappings or data-processing code. Rebuilds are an explicit one-time audit, not part of run_all.R.", "",
    "- reference_sites.csv: one declared reference identity per primary series, including excluded sites. reference_years records surveyed PISCO site/year footprints intersected with the finite harmonized primary history. For MBON and Landsat, years describe harmonized coverage under declared membership, not reconstructed individual reference contributions. included, excluded, exclusion_reason, identity_verified, contribution_audit, coverage_kind and registry_version document these distinctions. Mapping, footprint and exclusion-definition MD5 columns record the source snapshot.",
    "- reference_primary_keys.csv: every current primary series, including MPA, Taxa, Source, Resp, Region, AnalysisType, Model and exact observed/before/after year lists. Missing declarations remain explicit unmapped identities, not dropped series.",
    "- mapping_provenance.csv: role, kind, path, md5. analysis_input rows bind the registry to both harmonized lnRR files and sumstats_final.csv. derived_registry rows checksum the vendored site/footprint and key files. source_audit rows retain original mapping, survey and exclusion/processing-code checksums for provenance only; routine analyses do not open these paths.", "",
    "Script 27 stops before fitting if a required input or derived-registry MD5 differs, a current primary key/model/year classification is absent or changed, or reference-year flags/coverage are inconsistent. It cannot detect new changes in an unavailable upstream workspace: an upstream data refresh requires an explicit re-export, followed by updating the repository input-checksum manifest. No raw survey rows are vendored.", "",
    "PISCO uses exact declared site IDs, swath campus/zone and May-October filters, and fish campus/non-CAN filters (no fish seasonal filter). Source exclusion constants are parsed individually; the upstream pipeline is never sourced. KFM exclusions apply to non-sheephead series; LTER uses its declared exclusion list. Exact shared IDs connect groups; coordinates or outcome similarity never create identities. Landsat MPA labels are not global polygon IDs, so sharing/overlap remains unverified. PISCO biomass-specific missingness and MBON individual site contributions are not reconstructed.", "",
    paste0("Registry rows: ", nrow(reference_sites), "; primary series: ", nrow(keys), "."), "",
    "```sh", "Rscript code/R/export_reference_registry.R /path/to/Donham-Stier-CA-MPA-Data-2026", "```"),
    here::here("data", "harmonized", "REFERENCE_REGISTRY.md"))
  cat("Exported", nrow(reference_sites), "reference rows for", nrow(keys), "primary series; tracked-input validation passed.\n")
})
