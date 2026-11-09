# Consolidated audit — November 2, 2026

Use audit_master.csv as the consolidated comparison record (206 candidates;175 retained review histories). audit_status_by_comparison.csv tracks all nine audit domains. The viewer now includes Audit status, Station checks and Decisions tabs. Latest annual coverage is labeled audited_; original inventory coverage remains separately recorded. Station/reference findings are integrated, but depth/timing, expanded size support and final model integration remain open. No manuscript changes.

[Latest station/reference findings and figures](../all_program_station_refits_20261102/index.html). Historical November1 source and email checks below retain their original dates; no new email search is claimed.

# Dataset, site and comparison audit

November 1, 2026. Review only. No manuscript, supplement, canonical data, models or source records changed. No messages sent. No inclusion policy was adopted.

## Summary

The audit inventories all four programs used in the paper: PISCO, KFM, SBC-LTER and Landsat. It reads nine biological source CSVs plus the annual satellite delivery (8,865,211 source rows; the large integrated files are narrowed to focal taxa after reading). Supporting coefficient, size-source and satellite reconstruction manifests are indexed separately. These files are the local delivered releases, not a claim that no newer release exists.

The tables retain all 337 mapped candidate station assignments, expanded to 1,856 source/response/method rows, all 206 reconstructed candidate comparisons, and all 143 current manuscript histories. Some stations occur in more than one MPA assignment; assignment counts retain cross-program rows from the inherited master map and are not counts of PISCO survey sites. Rows are not independent replicates. There are 1,839 candidate raw same-year pair records, before ecological exclusions or depth/date qualification.

| program | current_paper_histories | reconstructed_candidates | researcher_excluded | retained_review_candidates | mapped_station_assignments | observed_station_assignments | retained_with_three_prior_years | retained_passing_proposed_recovery_duration |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| KFM | 29 | 49 | 14 | 35 | 43 | 43 | 19 | 35 |
| LTER | 11 | 18 | 0 | 18 | 22 | 22 | 13 | 18 |
| Landsat | 14 | 18 | 1 | 17 | 0 | 0 | 17 | 17 |
| PISCO | 89 | 121 | 16 | 105 | 272 | 208 | 16 | 78 |

The 175 retained candidates are the previous researcher/ecological scope roster, not a newly approved primary cohort. Three prior years and five recovery years are coverage screens only. Neither establishes baseline stability, habitat equivalence, identifiable recovery or sufficient size measurements. Designation-year observations are counted separately; existing model routing is not replaced by this table.

## Read the tables

- [Dataset-by-dataset measurement audit](DATASET_AUDIT.md), [screen results](measurement_design_scan.csv), and [additional checks before freezing](additional_analysis_checks.csv). These distinguish verified source flags from untested numerical consequences.
- [Dataset inventory](dataset_inventory.csv): source/version/hash, sampling design, measured quantity, units, missingness conventions and data-use status.
- [Resolved site groups](site_group_resolution.csv) and [group-to-station links](site_group_links.csv) distinguish group records from sampled stations and preserve exclusions.
- [Site and response coverage](site_response_coverage.csv): every mapped station assignment, source method/response, actual years, internal gaps, months, depth fields, coordinates and changes in sampling coverage. Raw coverage includes observations outside historical seasonal filters; it is not final biomass-valid coverage.
- [Comparison eligibility](comparison_eligibility.csv): all reconstructed candidates, current-paper membership, researcher exclusions, baseline/recovery coverage, reference evidence, size limitations and issue IDs.
- [Compact supplement table preview](supplement_table_preview.csv): one row per comparison, including excluded candidates. This is a review artifact, not inserted into the supplement.
- [Candidate pair coverage](candidate_pair_coverage.csv): separate protected/reference identities with same-year coverage and nearest-date diagnostics. These are raw opportunities, not approved pairs. Broad all-years depth-range overlap is not a matched-stratum test.
- [Full satellite delivery](satellite_full_delivery.csv) includes all 24 supplied MPAs and all 72 polygon slots, with 14 empty masks explicitly marked unavailable.
- [Satellite candidate coverage](satellite_candidate_coverage.csv) covers all delivered polygons, including those without reconstructed comparisons.
- [MPA metadata](mpa_metadata.csv): dates, current-law snapshots, historical uncertainty and take-by-taxon fields. A reference row inherits the comparison MPA date for analysis context, not a claim the reference was protected.
- [Unresolved issues](unresolved_issues.csv): specific missing evidence, analysis choices and next actions.
- [Supporting sources](supporting_source_inventory.csv), [unmapped source groups](unmapped_source_records.csv) and [annual sampling audit](raw_station_year_audit.csv) retain the evidence behind coverage flags.

## Main findings

1. **The programs measure different things.** Diver density, size-derived standing biomass and satellite surface-canopy biomass are kept distinct. Fish, benthic surveys and separate size collections do not inherit one another's coverage.
2. **Source-map namespaces matter.** Of 64 assignments without a direct source match, 43 are KFM/LTER rows embedded in the PISCO master map; they are not evidence of missing PISCO surveys. The subsequent hierarchy audit resolves 19 of the remaining 21 rows as grouped-site labels linked to 40 existing stations; only Trailer and Willows remain unassigned. Group labels carry no additional observations; identity and exclusion are separate. [Station exceptions](station_link_exceptions.csv) and [resolved hierarchy](site_group_resolution.csv); the earlier [legacy-name-only check](alias_evidence.csv) is historical distinguish these cases. Provider campaign years also differ from calendar years in 130 distinct metadata cells; [the date ledger](survey_calendar_campaign_audit.csv) preserves both.
3. **Long time series can conceal changing sampling.** Raw annual station metadata now identifies changes in months, zones and unit labels. Changes are flags, not automatic exclusions: non-fixed transects can have legitimate changes in identifiers. Depth summaries are record-weighted descriptions, not area-weighted habitat estimates.
4. **A year with no record is not a biological zero.** Counts, explicit zeros, negative sentinels and missing measurements remain separate. Incomplete full dates retain a recorded year/month. Repeated selected-field values are never silently removed.
5. **Some reference rationales remain incomplete.** Scott supports specific PISCO habitat exclusions; Jenn supports the Naples reference choice but reports no close match; Emily confirms LTER consultation but not each station's rationale. This does not validate every inherited pairing. The user's newer exclusion of Gull KFM supersedes older documents that retained it; Gull PISCO remains a separate candidate.
6. **Size support is uneven.** Existing detailed lobster/purple-urchin inventories cover 15 canonical histories, with sparse/shared pools explicitly linked. This cannot certify expanded candidate years. KFM counts do not become a full biomass history simply because a limited separate size survey exists.
7. **Satellite mask discrepancies are resolved in the corrected delivery.** Existing verification covers 24 MPA and 34 nonempty reference arrays. This audit reads the current annual delivery and links that verification rather than pretending to repeat pixel reconstruction. Known empty reference slots remain unavailable, not zeros. Surface pixels are not assigned dive-depth fields.
8. **Five current protection labels were resolved.** The audit now distinguishes Point Vicente and Casino Point no-take SMCAs from the limited-take rules at Abalone Cove, Arrow Point–Lion Head Point and Crystal Cove. [Dated source links](protection_metadata_resolutions.csv) accompany each correction. These metadata corrections do not change the cohort.
9. **Historical metadata and proposed eligibility remain distinct.** Current protection labels do not reconstruct historical boundaries. The depth/date/five-year recommendation is still a proposal, not a filter retrospectively selected to preserve an effect.

## Primary-analysis freeze

See [freeze readiness](FREEZE_READINESS.md) for required decisions and checks, conditional comparison-level limits, and items that can remain disclosed limitations. Resolving a site name does not adopt new analysis rules.

## Questions still open

| issue_id | scope | issue | type | who_or_where |
| --- | --- | --- | --- | --- |
| U01 | KFM | Permission and redistribution | Contributor confirmation | Scott Gabara / park data manager |
| U02 | KFM and LTER | Physical depth/habitat match | Missing metadata | Program investigators / field metadata |
| U03 | LTER | Individual reference rationale | Contributor confirmation | Emily / relevant LTER investigators |
| U04 | PISCO Naples | Community matching evidence | Missing supporting analysis | Jenn / Spiecker technical-report materials |
| U05 | PISCO | Three missing source measurements | Missing measurements | Original field records / program data managers |
| U06 | PISCO | Season/depth matching and reference dependence | Analysis decision | Adrian / analysis review |
| U07 | KFM | Predator biomass support | Missing measurements / population match | Scott Gabara |
| U08 | PISCO/KFM/LTER | Sparse and shared size pools | Analysis work | Analysis pipeline |
| U09 | All | Historical protection and sampled boundaries | Missing historical metadata | Official archived designations / program records |
| U10 | PISCO | LECHUZA coordinates | Conflicting metadata | Source mapping / program data manager |
| U11 | Landsat | Reference design and interpretation | Design limitation | Analysis review / Tom if historical detail needed |
| U12 | All | Unmapped source records | Scope/crosswalk review | Analysis pipeline then contributors if needed |
| U13 | Size surveys | Incomplete exact dates | Metadata limitation | Source size metadata |
| U14 | All | Repeated records | Record identity screen | Source schema / field keys |

Detailed findings, next actions and evidence are in unresolved_issues.csv. Do not send new questions already answered: Li specified SBC127.5 Table2 Model1; Jenn explained the Naples matching rationale; Scott explained Harris/Cuyler and Yellowbanks; Emily explained intended station continuity. The remaining needs are narrower: missing measurements, unrecovered supporting records, historical boundaries, data permission and final analysis choices.

## Email and source checks

Targeted Gmail searches on November 1 covered recent replies from Emily, Jenn, Scott Hamilton, Scott Gabara, Li and Tom. The latest relevant reference replies found were October 29; no newer response to these specific questions appeared in the searched results. Scott's full permission reply was reread: park permission and redistribution terms remain to be confirmed. Li's full reply confirms the fixed kelp coefficients. Earlier reference replies were read in the preceding depth audit; its private dated evidence record is retained in the manuscript repository. Search coverage is not proof that no other private conversation occurred. No private email bodies are copied into this analysis report.

Historical inclusion labels in the station table are explicitly named historical; they do not override the current researcher-scope decision column. The older 145-history coverage snapshot is not used as the current primary roster; current membership is rebuilt from the live 143-history result keys. Exact raw sampling dates, source program and method are kept separate from the provider campaign year.

## Reproduction and boundaries

Run audit_all_dataset_sources.py, build_dataset_inventory.py, extend_dataset_audit.py and report_dataset_inventory.py from code/tools, then verify_dataset_inventory.py. The scan is read-only and writes only this audit directory. Original source hashes and table input hashes are recorded. Auxiliary files are indexed, not automatically certified by their filenames. No raw individual animal records or emails are copied into the deliverable.

The audit is a systematic inventory and comparability screen across the paper's datasets. It does not certify every historical field convention, reconstruct unprovided measurements, validate every site pairing or refit every biological result. Unresolved items are visible instead of being treated as passed checks. The next analysis should use this inventory to freeze a final cohort and verify response-specific survey and size support before any manuscript update.
