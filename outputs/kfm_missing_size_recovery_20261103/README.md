# Five KFM comparison-years with missing sheephead sizes

November 3, 2026. No additional measurements recovered from the checked files. The existing 52-year restricted biomass review, density analyses, researcher exclusions and paper are unchanged. The alternative below is a coverage calculation, not a new fitted analysis.

## What the sources show

| Comparison-year | Missing size support | Source evidence | What a matched-visit subset would retain |
|---|---|---|---|
| Anacapa 2007 | Admiral’s Reef, August 8 | Sizes exist July 13, 26 days earlier. Scott’s export and the original archive agree. Separate August roving counts recorded sheephead but provide no lengths. | All four sites, 6 of 7 visits; 24 of 28 count records and 22 of 28 observed fish. Omits the August visit and retains July sizes with July counts. |
| Harris Point 2010 | Wyckoff Ledge, July 29 | Both size files contain other fish that day but no sheephead lengths. A separate roving observer recorded a male sheephead. No same-year size pool for this site. | Not recoverable while keeping all three sites. Dropping Wyckoff would change the reference composition. |
| Scorpion 2007 | Scorpion Anchorage, September 5 | Both size files contain other species that day. Sheephead lengths appear October 4, 29 days later. A September roving observer recorded a female sheephead. | All seven sites, 13 of 14 visits; 52 of 56 count records and all 12 observed fish. Removes a zero-count visit; this is still a sampling change. |
| Scorpion 2022 | Cavern Point, Little Scorpion and Pedro Reef, May 18–19 | No fish-size records for those May visits in Scott’s export. Sheephead sizes match later count visits on August 22 or September 13, 95–118 days later. | All seven sites, 7 of 10 visits; 28 of 40 count records and 43 of 59 observed fish. This omits spring visits. |
| Scorpion 2023 | Scorpion Anchorage, May 9 | Other fish have size records that day, but no sheephead lengths anywhere in that site-year. | Not recoverable while keeping all five sites observed that year. Dropping Scorpion Anchorage would remove a protected site. |

A count record is a source survey replicate, not an independent site. Counts of fish across visits describe the records retained, not unique individuals.

The complete Scott export was checked across all species, stages, physical-site aliases and available dates. The 2007 and 2010 checks also used the original Kushner archive’s individual fish-size and roving-count files; extracted files are byte-identical to their ZIP members. Separate roving counts corroborate visits and fish presence but cannot substitute for visual-transect density or provide missing lengths. An absence from the size file is not proof that sheephead were absent.

Emily’s project data and GitHub-folder copies of the integrated fish and combined size files are each byte-identical to the current copies. The combined size file contains 265,095 KFM benthic/kelp rows and no KFM sheephead sizes. It cannot fill these gaps. Recent missing visits have no independent 2022–2023 individual-size source in the checked project files beyond Scott’s export; this audit does not establish that no additional records exist elsewhere.

## Keep or lose

- **Current rule: 52 comparison-years.** Every count visit has same-site size support within seven days. Retain this result as the current restricted biomass review. All five omitted years still support density.
- **Matched-visit sensitivity: 55 comparison-years.** Apply the same visit filter across the full cohort, retain every original site in each comparison-year, and renormalize the original count-record weights within each site. This could add Anacapa 2007 and Scorpion 2007 and 2022. It changes the observations and annual averages; call it a matched-visit analysis, not recovered annual biomass.
- **Harris Point 2010 and Scorpion 2023 remain unsupported.** Sizes would have to be borrowed from another site or year, or sites removed. Neither is adopted here.

The matched-visit alternative is not automatically a better comparison. At Scorpion 2022, retained protected visits range from June 23 to August 22, while references are September 12–14. At Scorpion 2007, removing a zero-count September visit may also condition the sample on where sheephead sizes were obtained. These are reasons to show a sensitivity rather than silently add the years. Anacapa’s retained reference visits span July–September around an August protected visit; dates remain explicit.

I recommend keeping the 52-year review unchanged and, if additional coverage is useful, testing the 55-year matched-visit alternative separately on density and biomass. Do not widen the seven-day rule to 26–118 days, relabel dates, or treat zero counts as zero uncertainty to obtain more years. No missing measurements were reconstructed or imputed.

## Evidence and reproducibility

- [Seven missing visits and source evidence](missing_visit_evidence.csv)
- [Five comparison-years: records and fish retained](five_year_options.csv)
- [All years under both coverage rules](all_year_coverage_options.csv)
- [Candidate matched-visit responses and seasonal coverage](matched_visit_subset_responses.csv)
- [Candidate visit weights](matched_visit_subset_weights.csv)
- [Available size events](available_size_events.csv)
- [Drive-copy checks](drive_copies.csv)
- [Input hashes](input_manifest.json)

Run `python3 code/tools/audit_kfm_missing_size_recovery.py`, then `python3 code/tools/check_kfm_missing_size_recovery.py`. Candidate responses are point estimates only; no uncertainty interval, model refit or manuscript change is implied.
