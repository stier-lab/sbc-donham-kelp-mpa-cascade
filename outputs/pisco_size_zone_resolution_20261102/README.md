# PISCO size-zone flags resolved to seasonal filtering

All 23 flagged species/year/role groups arise because the counts were filtered to May–October but their size pools were not. All 43 underlying station/year/species/zone groups have raw count surveys in the same zone on the size-measurement dates. Those dates are in March, April, November or December. No flag requires a zone relabeling or a missing count record to explain it. This supersedes the earlier unresolved explanation of the 23 flags in A61; it does not establish that every count record in the delivery is complete.

For example, Abalone Cove's 2014 MID-zone sizes and counts were both measured on December 22; only the count records were removed. Point Vicente West's 2013 shallow-zone sizes and counts were both measured on April 24, leaving only June DEEP-zone counts after filtering. Three Palms East's 2021 shallow-zone surveys were on November 2, twelve days after the retained October DEEP survey.

## Same-observation checks

I tested five changes separately: removing only the flagged zone measurements; requiring same-site sizes; requiring same-site/zone sizes; applying the existing May–October rule to the size pools; and requiring site, zone and season together. These are diagnostic restrictions, not adopted exclusions. The count data, offsets, researcher exclusions and models remain fixed. Each revised fit has a current-input fit on exactly the same retained years. Annual mean masses are reweighted from source frequencies and used to scale current biomass, so these remain conditional sensitivities rather than a complete production rebuild.

Applying May–October to the size pools and requiring site/zone/season gives identical recovery estimates for these six histories. All six retain matching-season support in most years, but three species-specific paired biomass years would be lost: Abalone Cove lobster2014, Abalone Cove purpleurchin2018 and Point Vicente redurchin2013. A whole paired biomass year is withheld if a positive role lacks eligible sizes; density remains supported. Zero counts remain zero.

The largest changes in these comparisons come from losing years, rather than from revising average animal size within the same years. Abalone Cove lobster changes from0.232 to0.190 simply by comparing the ten supported years, then to0.184 when sizes are restricted. Its temporal interval already includes zero after the year is removed, before changing sizes. Point Vicente redurchin changes from0.016 to−0.061 when2013 is removed, then to−0.063 with the size restriction; intervals include zero throughout. On matched years, the largest absolute slope adjustment from seasonal size matching is0.0091 log-ratio units/year. The other four histories show small changes in slope, with the same signs.

## Recommendation

Use a consistent survey-selection rule for counts and their size measurements. These findings do not independently justify May–October as the right biological season: some excluded visits fall only days outside its boundary. Before freezing the primary analysis, choose the count-survey timing policy and apply it to the corresponding size support. Do not exclude these sites because of the earlier zone flags, or borrow out-of-season sizes silently to preserve a biomass year. Retain valid density when biomass lacks support. The existing user hold on selecting a new seasonal cutoff remains in force.

The figure separates the current full history (A), current inputs on the years retained by the seasonal size check (B), and season-matched sizes on those same years (C). Lines are conditional temporal-model95% intervals, not total biomass uncertainty. These are post-protection recovery slopes using calendar years, not estimates of before–after protection effects. No pooled or downstream pathway model has been changed.

## Verification

All23original flags and43underlying station/species/zone groups were traced to raw files. Every flagged size date has a matching raw swath zone/date; all those count dates are outside May–October. Independent R checks reproduce66slopes and temporal standard errors. All487protected paper/result files are unchanged. Source hashes, pool counts, exact survey dates, coverage changes and matched-year fits are saved here. Code: `code/tools/audit_pisco_size_zones.py`, `code/tools/verify_pisco_size_zones.R`, `code/tools/report_pisco_size_zones.py`.
