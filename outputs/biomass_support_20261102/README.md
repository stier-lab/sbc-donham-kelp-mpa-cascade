# Biomass and size-sample support — November 2, 2026

Audit covers all94 retained biomass histories:39 lobster/urchin histories,38 direct-fish or kelp histories, and17 satellite canopy histories with no diver size conversion. It preserves researcher exclusions and the current annual review observations. Paper and primary results unchanged.

## What is supported

All1191 positive protected/reference-year values in the39 animal histories have measured sizes under the current rules. There are1222 total role-years including zeros. The old15-history size inventory agrees exactly on overlapping measured-animal counts. For38 other fish/kelp histories, every positive station-year has source measurement rows. Direct measurements and resampled size pools are different evidence: one measured fish among one fish observed is not the same problem as using one fish to represent many unsized animals.

No new unsupported KFM predator biomass has been created. The two LTER lobster histories retain the five existing same-site/date missing-length corrections. These are estimates with explicit donor support, not recovered measurements. Satellite canopy has separate mask/imagery limitations, not animal-size uncertainty.

## Sparse and shared pools

23 positive role-year pools contain one measured animal:6 lobster,15 red-urchin and2 purple-urchin pools.91 have fewer than10 measured animals. These are descriptive flags, not adopted exclusion thresholds. Some pools use the same original record across histories, so these are not23 independent sampling failures.

Lobster singleton pools: Abalone2013 reference and2016/2017/2023 protected; Matlahuayl2023 protected; Point Vicente2012 protected. Campus Point red urchins are the next priority: LTER has7 singleton role-years and PISCO has5; the protected2016/2017/2018/2022/2023 sizes are shared across the two programs. Naples LTER red-urchin reference2010 is another singleton, and Matlahuayl PISCO red urchins have two. The purple-urchin singletons are LTER Campus2008 and Naples2001 protected.

571 positive role-years reuse size records in another retained history.325 have no donor animals from their own monitoring program;329 have no donor animals from the actual count stations in that role/year. That follows the existing MPA/status/year/species pooling rules; it is not newly introduced borrowing. PISCO urchins use diameters>=25mm, whereas KFM/LTER use unfiltered sizes. Cross-program donors may provide useful information, but their representativeness must be justified and they do not provide independent evidence across programs. Same-year pool matching does not establish same visit, depth or habitat. Missing KFM dates in the combined delivery prevent precise date matching; legacy source-date recovery remains a separate task.

663 positive role-years have one recorded survey-group identifier. Missing dates/transect identifiers can collapse multiple visits, so this must NOT be described as663 confirmed single-dive samples. A single-animal pool cannot estimate within-pool size variability by resampling that animal.

## Missing stipes cannot be zero biomass

The source again identifies two kelp plants with missing stipe measurements at UNDERWATER_ARCH in2023 (Point Vicente reference). The current station biomass is0, so checking only positive biomass would miss it. This is unknown biomass, not a demonstrated absence of kelp. Density can retain the measured plants.

An isolated conservative sensitivity withholds that entire paired-year instead of silently dropping a reference station and changing weights. The Point Vicente PISCO kelp recovery slope changes0.2425 [−0.2093,0.6943] over10 post-years to0.3420 [−0.1499,0.8338] over9. Both intervals include zero. This is a missing-data treatment, not a recovered biomass estimate or primary change. Full station-aware modeling could retain other measured stations with an explicit prediction population; that is not adopted here.

## Kelp and direct measurements

PISCO and LTER fish use lengths measured with counts; no positive retained station-year has a missing-length flag in the screened source events. LTER kelp uses direct frond measurements; PISCO kelp uses counted plants times their stipe measurements, with the exception above. These checks do not validate allometric-coefficient uncertainty, field detection or every historical aggregation decision.

KFM kelp converts counted plants using separately sampled plants from the same site/year. One positive station-year has only one sampled plant: Anacapa site a-k-33 in2007. Same-site/year support is verified; matching the actual count and size-survey dates and propagating plant-sampling uncertainty remain necessary. Equal-valued size records have not been deduplicated.

## Recommended handling

1. Keep density where counts are valid. Treat the Underwater Arch stipes as unknown, with the paired-year omission shown for review; do not substitute zero biomass.
2. Retain the existing direct-size calculations and documented LTER lobster corrections. A small number of directly counted animals is not automatically a failed size pool.
3. Prioritize uncertainty propagation and stress checks for the23 singleton pools, including shared Campus Point red-urchin animals. Preserve each original donor's identity across every use. Conditional bootstrap precision is not evidence of representativeness.
4. Test the consequences of cross-program/shared sizes on matched years before deciding biomass eligibility. A same-program-only rule would leave325 currently positive role-years without donors; do not silently replace them with zero, extrapolated sizes or density mislabeled as biomass. Keep these biomass histories provisional where representativeness remains unverified.
5. Extend the old lobster/purple-urchin uncertainty analysis to the current roster and red urchins; check KFM kelp plant sampling separately. No final biomass admission rule is chosen from outcome direction.

## Verification and limits

104125 donor occurrences independently checked against original source-row year/site/count/size, frequency sums and PISCO size threshold; sharing independently reconstructed. All earlier overlapping pool counts match.487 protected manuscript/result files remain unchanged. The screening covers source rows within retained site-years and stated protocol filters; it does not certify exact same-visit matching, new donor representativeness or allometric uncertainty. Raw donor membership stays in the ignored private directory, with public aggregate tables linked here. No emails or push.
