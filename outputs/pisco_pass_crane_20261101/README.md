# PISCO single-pass records and CRANE2004

November1,2026. Review only. No paper, canonical data, researcher exclusions or primary fits changed.

## Recommendation

Keep supported CRANE2004 records at retained sites. For the combined bottom-plus-upper-water response, use completed BOT+(MID or CNMD) survey units and preserve measured single-pass data separately. Do not interpret missing passes as zero fish, combine passes separated by weeks/months into one visit, or exclude an entire site because one unit is incomplete. Keep the repeat-effort normalization proposed in the preceding audit. These recommendations follow survey documentation and measurement completeness, not slope direction.

## CRANE evidence

The2015 South Coast baseline report states that its approach uses PISCO protocols also used by CRANE, describes30m by2m fish belts and simultaneous bottom/midwater passes, and integrates CRANE2004 historical surveys. The public upstream processing script retains CRANE and explicitly removes the duplicate VRG copy of SBI2004. The delivered source has no VRG SBI2004 records remaining.

The raw file contains2,803 CRANE records across19 sites and288 survey units. Within the retained paired comparison windows,144 units from9 stations contribute to Anacapa Island SMR2003, Santa Barbara Island SMR and South Point SMR. Every one of those144 units has BOT+MID or BOT+CNMD. The four incomplete CRANE units elsewhere in the full delivery are outside these retained comparison windows. The method label therefore does not justify discarding these three2004 paired observations. Keep the label and existing researcher exclusions; this does not reopen excluded sites.

Sources: [South Coast baseline report, printed pp15–16 and18–19](https://caseagrant.ucsd.edu/sites/default/files/SCMPA-27-Final-Report_0.pdf); [pinned upstream merge code, lines111–113 and323–328](https://github.com/zoekitchel/PISCO-CRANE-processing/blob/39abf7e97e9b2720e4cb5110700312c39e22e611/scripts/from_PISCO/MLPA_data_merge_fish_2023.R). Downloaded report and code are retained with hashes. Shared methods support comparability; they do not prove identical habitat or every historical field condition.

## Single-pass evidence

The previous216-unit candidate-source screen reduces to67 unique single-pass units within retained comparison windows across9 MPAs (134 response-linked rows because each unit is linked to density and biomass). Three units at Cluster Point South2006 explicitly lack bottom surveys because of poor visibility. The other64 have no decisive explanation in the checked files. No opposite pass was recovered on the same date under another method or a whitespace/case variant of the level code.

Across the broader216-unit screen,33 units have an opposite pass on another date12–253 days away. These are not recovered simultaneous surveys. Some belong to different provider campaigns within a calendar year. Do not stitch them together merely because the transect number matches.

The full delivery has18,092 explicit NO_ORG records, all count0 with LOOKEDFOR=yes. This confirms a mechanism for recorded zero-fish passes; it does not establish perfect use of that code in every historical survey. Missing-pass records therefore remain unknown rather than automatically becoming zero. Notes such as no canopy or one side done do not, by themselves, prove whether the missing bottom or midwater pass was completed.

The SBI2005 notes also report overlap between Arch Point central transect4 and the neighboring Arch Point transect1. That year is already unavailable in the complete bottom+upper comparison. Preserve this overlap warning before any bottom-only reuse; do not silently delete equal-valued fish rows.

## Coauthor-email check

Targeted searches covered CRANE, the literal method code, CNMD, single-pass, canopy, midwater and bottom-only language, including Jenn, Emily and Scott Hamilton. No direct coauthor instruction to exclude CRANE or interpret a missing pass as zero was found. Recent Jenn/Emily replies were reread: Jenn explained that KFM/PISCO observations are separate and that the supplied raw fish sizes need explicit preprocessing; Emily emphasized same-site paired years and historical ecological exclusions based on Scott's advice. Neither resolves the missing-pass meaning. Search absence does not prove no private conversation occurred. Existing dated correspondence records remain authoritative; no email was sent and no private email bodies were copied here.

## What changes

Keeping CRANE and using complete passes retains417 of420 original history-year rows across34 sheephead histories at17 retained PISCO comparisons. The three unavailable response rows are SBI2005 density and biomass, plus South Point2005 biomass from the earlier missing-length problem. This is not three excluded sites, and it is not417 independent samples. Anacapa, SBI and South Point2004 remain available.

The incomplete-pass comparison is fitted on identical years in both arms, after repeat normalization. Five recovery comparisons change slightly: Abalone Cove, Gull Island, Point Vicente, Scorpion and South Point. Point Vicente's before–after log-ratio change increases from0.4167 to0.4714 for biomass and0.0894 to0.1496 for density. The original missing-year consequence is disclosed separately rather than included in this same-year effect comparison.

The CRANE contrast deliberately compares complete-pass histories with and without2004. Keeping2004 changes biomass recovery slopes from−0.0218 to−0.0023 at Anacapa,0.0351 to0.0501 at SBI, and−0.0191 to−0.0188 at South Point. These are calendar-year log-ratio slopes, not mean protected-site advantages or the manuscript's mixed effect-size definitions. All-program pooled recovery intervals still span zero. No new significance claim is made. Full density results, intervals and both contrasts are in effect_changes.csv and pooled.csv.

## Remaining evidence and verification

We can make a consistent complete-measurement decision now without guessing missing observations. To reinstate the64 unexplained single-pass units in a combined response, the useful question for Jenn/the data custodian is whether each absent opposite pass was unperformed, lost in export, or represented by another documented protocol, and whether original survey sheets or completion flags survive. No further broad CRANE inclusion question is needed unless contrary source evidence appears. Missing depth/season/reference-policy questions remain separate.

Reproduce: audit_pisco_pass_crane.py, compare_pisco_pass_policy.py, refit_pisco_pass_policy.R, report_pisco_pass_policy.py, verify_pisco_pass_policy.py from code/tools. Individual estimates and calendar-gap HAC3 SEs independently checked for224 fits;32 pooled models fitted with2 recorded BFGS fallbacks. Weakly identified pooled before–after biomass does not get a robust significance claim.487 protected artifacts and all35 artifacts from the preceding fish-effort audit unchanged.
