# Matched targeted/non-targeted fish diagnostic

Specified 23 October 2026 before estimating fish responses. This optional module
does not change primary inputs, fits or manuscript figures.

Use the 15 retained PISCO sheephead-density MPA comparisons and their recorded
analysis years. Use the declared, nonexcluded reference sites and mapped inside
sites. Freeze the roster of sites observed at least once in those years; do not
choose sites based on fish abundance. Require every roster site in a retained
year. Export missing sites and all coverage exclusions. This strict diagnostic
may retain fewer years than the primary analysis's changing site roster.

Fixed panel: targeted sheephead and kelp bass; non-targeted senorita, blacksmith,
rock wrasse, garibaldi, kelp surfperch, painted greenling and giant kelpfish.
Classifications are checked against the local PISCO attribute table and Hamilton
et al. (2010) Figure 4. The first four non-targeted species and both targeted
species have warm affinity. Report warm-only non-targeted comparisons as the
main comparison and all seven non-targeted species as sensitivity.

Count density, not biomass: retain observed counts and reconstruct absences on
recorded standard PISCO/VRG fish transects (LOOKEDFOR=yes), not absent visits.
Use UCSB/VRG only and exclude the distinct CRANE method. Primary sampling frame
requires BOT and MID observations in the same site/year/zone, averages transects
within each visit/level, visits within year/level, then sums BOT and MID means.
Average zones and sites equally. Exclude CAN and CNMD to avoid changing vertical
coverage with canopy presence; bottom-only is a sensitivity on exactly the same
sampled cells. This is a standardized diagnostic, not a reconstruction of the
paper's all-noncanopy sheephead index. Record the zone roster by site/year.

Use base species codes, excluding distinct YOY codes. Primary counts include all
sizes under these codes. Also remove fish at/below the recorded YOY size (or the
juvenile cutoff when no YOY cutoff exists). Call this a size-filtered sensitivity,
not a verified adult or universal YOY classification. Drop whole affected survey
cells from both versions if any panel species has missing count or positive count
with unknown size, ensuring identical denominators and years. Include offsets of
0.1, 0.5 (main) and 1 fish per nominal 60 m2 footprint before logging.

Within each pair/year calculate species log inside/reference contrasts. Average
species equally within the two targeted fish and within the prespecified
non-targeted panels. Fit their difference directly, retaining group and individual
species estimates. These are mean species log contrasts, not total community
biomass or identical ecological guilds.

Omit designation year. Before/after step comparisons require >=3 paired years
strictly before and >=5 strictly after. Label short baselines explicitly; do not
extrapolate fitted pretrends. Detrended pre-designation tracking requires >=6
years and nonconstant residuals. Separately fit post-designation calendar-year
linear trends with >=8 years and >=8 years between first and last observations.
Report calendar-gap HAC (lag 2) intervals and lag 1/3 sensitivity. No model
selection. Temporal fits are descriptive, not causal protection effects.

Equal-MPA pooled effects use the existing delete-shared-site-group jackknife
with t(G-1) intervals; connect pairs sharing any physical site (site_new as well
as raw site ID). Suppress pooled intervals below three groups. Do not treat
annual samples or species as independent MPA replicates. Report species-level
results and offset, size, sampling-frame and leave-group-out sensitivity. No
claim of reef equivalence follows from nonsignificant results. Non-targeted
fishes can respond to kelp, temperature, predation and protection itself.
