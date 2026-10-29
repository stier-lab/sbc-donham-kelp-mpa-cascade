suppressPackageStartupMessages({library(dplyr); library(tidyr)})
source("code/R/primary_uncertainty_utils.R")
source("code/R/fish_siting_utils.R")
source("code/R/fish_recovery_utils.R")
read <- function(x) read.csv(paste0("tables/table_s_fish_sampling_", x, ".csv"))
stopifnot(read("status")$status == "complete", read("reconstruction")$maximum_difference < 1e-10)
for (name in c("protected", "inputs")) {
  d <- read(name); expected <- if (name == "protected") d$before_md5 else d$md5
  stopifnot(identical(unname(tools::md5sum(d$path)), expected))
}
coverage <- read("coverage"); annual <- read("annual"); fits <- read("fits")
used <- read("used_zones"); zones <- read("zone_roster")
depth <- read("depth_audit")
stopifnot(all(depth$missing_depth == 0), all(depth$max_deviation <= depth$width + 1e-10, na.rm = TRUE))
old_cover <- read.csv("tables/table_s_fish_recovery_coverage.csv") %>%
  filter(scope == "primary_years", frame == "bottom_mid", rule == "core80", eligible)
old_roster <- read.csv("tables/table_s_fish_recovery_roster.csv")
for (i in seq_len(nrow(coverage))) {
  c <- coverage[i, ]; yy <- as.integer(strsplit(c$years, ";", fixed = TRUE)[[1]])
  stopifnot(identical(c$eligible, fish_recovery_eligible(yy, old_cover$start[match(paste(c$endpoint,c$MPA),paste(old_cover$endpoint,old_cover$MPA))], c$endpoint)))
  if (!c$eligible) next
  scenario <- paste("primary_years|bottom_mid|core80", c$endpoint, sep = "|")
  rr <- filter(old_roster, .data$scenario == .env$scenario, MPA == c$MPA)
  u <- filter(used, endpoint == c$endpoint, MPA == c$MPA, restriction == c$restriction)
  zr <- filter(zones, endpoint == c$endpoint, MPA == c$MPA)
  stopifnot(setequal(u$site, rr$site), setequal(u$year, yy), !anyDuplicated(u[c("site","year","zone")]))
  for (s in rr$site) for (y in yy)
    stopifnot(setequal(u$zone[u$site == s & u$year == y], zr$zone[zr$site == s]))
}
keys <- c("endpoint","MPA","restriction","version","size","year")
w <- annual %>% pivot_wider(names_from = series, values_from = contrast)
stopifnot(all(abs(w$Difference_warm - (w$Targeted - w$Nontargeted_warm)) < 1e-10))
matched <- annual %>% select(all_of(keys), series, contrast) %>%
  pivot_wider(names_from = version, values_from = contrast)
stopifnot(!anyNA(matched$restricted), !anyNA(matched$original_matched))
for (d in split(annual, interaction(annual[c("endpoint","MPA","restriction","version","size","series")], drop=TRUE))) {
  f <- fish_siting_fit(d, unique(d$endpoint), unique(d$start), 2)$result
  saved <- filter(fits, endpoint==d$endpoint[1], MPA==d$MPA[1], restriction==d$restriction[1],
    version==d$version[1],size==d$size[1],series==d$series[1])
  stopifnot(nrow(saved)==1, abs(f$estimate-saved$estimate)<1e-10, identical(f$years,saved$years))
}
r <- read("retention"); a <- read("retention_annual")
stopifnot(all(abs(r$difference - (r$excluded_mean - r$retained_mean)) < 1e-10))
stopifnot(all(is.finite(a$excluded_value)), all(is.finite(a$retained_value)))
audit <- read("retention_audit")
stopifnot(all(audit$n_years[audit$window=="early3"]<=3), all(audit$n_retained[audit$usable]>0))
map <- read.csv("tables/table_s_fish_siting_roster.csv")
panel <- read.csv("tables/table_s_fish_siting_panel.csv")
site <- read.csv("tables/table_s_fish_siting_site_annual.csv") %>% filter(frame=="bottom_mid",size=="all_sizes")
for (i in which(audit$usable)) {
  d <- audit[i, ]; rr <- strsplit(d$retained,";",fixed=TRUE)[[1]]
  stopifnot(!d$site %in% rr, all(map$status[map$MPA==d$MPA & map$site %in% rr]==d$status))
  yy <- as.integer(strsplit(d$years,";",fixed=TRUE)[[1]])
  for (metric in c("Targeted","Nontargeted_warm")) {
    codes <- if (metric=="Targeted") panel$species[panel$Targeted=="Targeted"] else
      panel$species[panel$Targeted=="Nontargeted" & panel$thermal=="Warm"]
    x <- site %>% filter(site %in% c(d$site,rr),year %in% yy,species %in% codes) %>%
      group_by(site,year) %>% summarise(value=mean(log(density*60+.5)),.groups="drop")
    expected <- mean(x$value[x$site==d$site])-mean(x$value[x$site %in% rr])
    actual <- filter(r, endpoint==d$endpoint,MPA==d$MPA,site==d$site,window==d$window,.data$metric==.env$metric)
    stopifnot(nrow(actual)==1,abs(actual$difference-expected)<1e-10)
  }
}
cat("PASS: original-density reconstruction, fixed site/zone rosters, depth restrictions, matched years, direct contrasts, refits and input hashes\n")
