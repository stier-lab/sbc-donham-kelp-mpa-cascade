# Table definitions and interpretation

- **dataset_id:** stable readable source identifier. One integrated delivery can contain multiple monitoring programs; methods stay separate.
- **site_record_id:** hash of program, assigned MPA, source station, role, dataset, taxon and method. It is a join key, not a biological replicate.
- **comparison_id / series_id:** stable MPA–program–taxon–response identity. Current-paper membership is independent of researcher-scope eligibility.
- **candidate_pair_coverage:** raw candidate protected/reference pairing with the same response and method. These are not newly selected reference assignments. Excluded sites remain visible.
- **n_years:** unique recorded calendar years in the specified row. Station rows include raw recorded measurements and survey opportunities before seasonal/validity exclusions. Comparison rows describe the reconstructed annual analysis input. Never substitute one for the other.
- **n_missing / missing_years:** absent years between that row's first and last observations; years before monitoring began and after it ended are not counted.
- **elapsed_years:** last minus first year. Inclusive calendar span is elapsed_years + 1.
- **longest_internal_gap:** longest consecutive sequence of absent years, not elapsed distance between observed years.
- **n_strict_before / n_designation_year / n_strict_after:** years strictly before, equal to, or strictly after the recorded protection/network year. These describe coverage, not adopted model routing; within-year protection dates can complicate the designation year.
- **three/four_strict_prior_years:** numerical baseline coverage screens. They do not establish stable trends, continuity or pBACIPS eligibility.
- **proposed_recovery_duration_pass:** at least five strictly post-designation years spanning five elapsed years. This remains a proposal; date/depth/season identifiability is separate.
- **survey_dates:** valid complete source dates. Missing full dates are not invented from a year/month or campaign label.
- **months:** known calendar months; preserved even when the day is missing. Separate size records cannot be asserted to match an exact count visit merely because the year matches.
- **years_nearest_dates_within60/90:** minimum date separation between the two raw source records in a calendar year. No campaign/zone/effort validation is implied. Use the dedicated depth policy audit for the stronger five-MPA rule.
- **depth_min/max / depth_record_median:** recorded positive depths. Medians in the annual audit are record-weighted, not spatial estimates; no biomass weighting. No single nominal depth substitutes for missing measurements.
- **month/zone/unit_set_changes:** descriptive differences between annual source records. Not an error flag by itself; non-fixed transects and focal-species-only records can vary legitimately.
- **zero_count_rows:** explicit source zeros. **missing_count_rows** are blank/non-numeric values. **negative_count_rows** preserve sentinels for review. None is silently recoded here.
- **positive_count_missing_nonpositive_size:** broad raw screen only. Size fields differ among taxa and kelp stages; this is not a count of unusable biomass estimates.
- **selected_field_duplicate_rows:** identical values in selected scientific fields, omitting some observer/taxon-code/replicate metadata. They are not certified duplicate records and are not removed.
- **size_sparse/singleton/shared:** existing detailed canonical lobster/purple-urchin inventory. Sample counts are animals represented by frequencies, not independent surveys; scope is explicitly limited.
- **historical_*:** October31 decision/metadata snapshot retained as evidence. Current researcher exclusions and newer explicit user decisions take precedence.
- **protection_class / take_restrictions:** dated classification evidence. Whole-MPA area, depth and network date are not substituted for surveyed habitat, station depth or first historical protection.
- **issue_ids:** links to unresolved_issues.csv. A broad program-level issue signals a check still needed, not proof that every comparison has that defect.
- **unmapped_source_records:** other regions/methods/identifiers in the delivered releases. Presence alone does not justify adding them to the study; absence from a map is not absence of a survey.

The source scan includes all taxa as survey opportunities in PISCO and the combined size file, but focal taxa only in large integrated MBON and dedicated fish files. Invertebrate raw codes are preserved in taxon_code_inventory.csv. Scientific-code aliases do not merge independent sampling methods or convert juveniles into adults. SNI records in an integrated delivery remain outside the paper's four-program comparison inventory.

## Resolved site hierarchy update

- **inventory_record_type:** distinguishes observed source/response records, grouped-site labels, cross-program map entries and the two unassigned map labels. Group records carry no additional observations.
- **group_names / linked_constituent_stations:** exact source-map hierarchy, not fuzzy name substitution. Stations remain separate physical sampling units.
- **resolved_group_scope:** identity is resolved independently of inclusion; five of the19 mapped assignments have separate exclusion reasons.
- **inherited_pisco_reference_exclusion / inherited_whole_mpa_exclusion / inherited_kfm_nonfish_exclusion:** separate source rule flags. Do not transfer a KFM non-fish rule to PISCO or silently override fish exceptions or newer researcher-scope decisions.
- **mapped_protected_sites / mapped_reference_sites:** now list observed candidate station IDs only, including explicitly flagged exclusions; these are not fitted membership. Original broader map labels are preserved in historical_* columns.
- **freeze_readiness:** separates decisions before rule freeze, checks before result freeze, conditional issues, disclosed limitations and data-use requirements. No proposed threshold is adopted by updating this inventory.

## Measurement/design audit

`dataset_audit.csv`: one interpreted review per dataset (D01–D10), separating findings, status, next action and evidence limits. `measurement_design_scan.csv`: full-delivery source screens before cohort filters; missing physical keys depend on the listed fields and do not mean unusable observations. `measurement_protocols.csv`: source/program/method/taxon-code strata; size is not the same quantity in every dataset. `repeated_unit_dates.csv`: reused physical-key/year groups with multiple distinct complete dates, not repeated animals. `pisco_fish_level_coverage.csv`: standard UCSB/VRG surveyed units from all recorded taxa with LOOKEDFOR=yes and non-canopy levels, so effort is not conditioned on sheephead presence. `pisco_fish_revisited_units.csv`: physical unit-years with more than one distinct date. Site-table raw flags apply only to observed PISCO-fish rows and do not alter eligibility. `additional_analysis_checks.csv`: six unresolved design checks, not six verified numerical errors. Area and length screens do not certify units or habitat equivalence. Exact full-row repeats may be equal-valued individuals and must not be automatically removed.
