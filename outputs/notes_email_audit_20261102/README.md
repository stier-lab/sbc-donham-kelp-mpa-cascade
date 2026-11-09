# Emily’s notes and project correspondence audit

November 2, 2026. This checks historical instructions and contributor clarifications against the current guide and selected processing code. It does not update the paper, adopt new exclusions or certify an end-to-end model rerun.

The [finding ledger](findings.csv) contains 32 topics, each with evidence, current status, a code/table check and the next action. The [KFM station key](kfm_reference_followup.csv) maps Scott’s latest concerns to physical stations and the existing response-specific exclusion scope. The [current guide](../../docs/ANALYSIS_GUIDE.md) and [questions table](../dataset_inventory_20261101/unresolved_issues.csv) remain the operating documents.

## What needs attention

- **KFM reference suitability:** Scott’s November 2 reply flags Webster’s Arch, Rodes Reef and the north/south San Miguel comparison. Webster’s Arch, Rodes Reef and Miracle Mile already have benthic exclusions, while fish are handled separately. Wyckoff remains a reference contributor in at least one primary history. Reconstruct the retained contributors by response and test the consequences before deciding further exclusions. The station key records historical mapping, not proof that every mapped site contributes to every fit.
- **KFM predator biomass:** Scott confirms that KFM does not collect lobster carapace lengths. Fish sizes were collected by one diver during a 30-minute roving survey. His new fish-size attachment needs reconciliation with count sampling; individual measured fish do not provide independent survey replication. The existing density comparisons remain distinct from biomass feasibility.
- **LTER kelp:** Li has finished updating SBC50 annual biomass. The source registry still identifies the March 2024 annual kelp delivery, and our code calculates biomass from fronds/stipes. Compare matched records and formulas against the updated package before deciding whether any values need changing. SBC127.5 Table2 Model1 remains Li’s specified fixed conversion source.
- **Baseline convention:** Emily’s latest reply clarifies three observations at the same stations, assessed separately for fish and benthic surveys, with the implementation year counted as before. Our revised calendar convention must be identified as a deliberate departure. Do not pool short records from different physical stations to manufacture the minimum baseline.
- **Integrated biomass uncertainty:** shared references, shared size donors, roving samples outside count transects, changing survey effort and KFM plant sampling remain substantive concerns. Source reconciliation and separate sensitivity runs do not finish the coordinated rebuild or uncertainty model.
- **KFM permission:** the latest inspected permission thread still requires park confirmation of the permit and redistribution terms. A new data attachment does not resolve this.

## What the older record clarifies

Emily’s notes explicitly identified several issues now visible in the guide: KFM sizes entering PISCO biomass estimates, different fish survey levels, changing KFM quadrat areas, and shared references. Current KFM processing uses recorded row-level area. Current PISCO production code sums non-canopy fish counts over a seafloor footprint; the approved repeat-survey and complete-pass changes still need coordinated integration. Area versus volume is a response definition, not a substitute for matching survey effort.

The current figure code uses a common protected/reference maximum, as Emily requested. The old baseline-trend and Cook’s-distance exclusions are historical choices; they do not override Adrian’s later direction to retain ecological scope and evaluate model eligibility and influence separately.

Jenn clarified Naples’ change from Campus Point to Arroyo Quemado and confirmed that KFM and PISCO observations are independently collected. Her underlying Naples ordination has not been independently reproduced. Hamilton’s reason for avoiding Crook Point and Tyler Bight is habitat comparability, not poor data quality. Emily’s LTER reference choices followed consultation, although station-specific reasoning remains incompletely documented.

Tom’s corrected masks and annual tables supersede older mask files. Earlier emails discuss centroid distances while later descriptions refer to boundary bands; the saved masks and verified code define the extraction. Diver and satellite observations still measure different footprints and kelp components. Do not silently re-create reference polygons from an approximate verbal rule.

The historical `Sites_Years.xlsx` has 120 rows and 45 columns. Its year cells store no values; coverage is encoded in fill colors. A value-only import would therefore misread the workbook. The current raw-record coverage tables remain the authority; recovering the old color key is needed before using that workbook as an independent numerical check.

## Scope and limits

- Searched project-subject correspondence back through February 2021, before this project’s 2022 initiation, and exhausted the Emily-specific project-keyword query (246 hits). The broader subject search covered 1,248 messages; these counts overlap and are not unique project messages.
- Retrieved 54 candidate threads containing 312 messages, then a targeted Naples follow-up thread containing three messages. These include sent replies, notifications, unrelated kelp projects and quoted duplicates. Relevant original messages and inline answers were read; unrelated scientific projects were excluded from this paper’s decision ledger.
- Inventoried 124 entries in Emily’s project folder. Extracted seven Word documents, retrieved eleven native Google Docs and all 220 returned comment threads plus replies. Read the substantive standalone notes and comment texts, and checked relevant historical methods/code passages. Old manuscript versions provide historical context; this is not a new sentence-by-sentence certification of every draft or an exhaustive rereview of every attached coauthor manuscript.
- Examined the historical coverage workbook, original code snapshot and selected current processing/figure code. Review comments are evidence to assess, not instructions to adopt automatically; exploratory proposals and obsolete journal limits do not create new analysis requirements.
- Search coverage cannot recover deleted mail, conversations outside Adrian’s mailbox, undocumented meetings, or unprovided raw records. Emily explicitly said she no longer had access to some old correspondence. Some exact report/code replication and habitat rationale therefore remain unresolved.

Full private correspondence and extracted notes are retained in the ignored `local_archive/project_notes_email_audit_20261102/`. Public-facing artifacts contain only methodological paraphrases and source identifiers. No messages were sent, no new data release was substituted and no numerical or manuscript output was changed.
