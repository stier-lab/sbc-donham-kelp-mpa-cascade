# CA MPA kelp-forest analysis — assistant guide

Updated 2026-11-09. Short operating guide for coding assistants. Humans start with `README.md`.
The full session log that used to live here (entries A44–A100, Oct 30 – Nov 4) is preserved verbatim in
`docs/archive/2026-11-09-cleanup/CLAUDE_session_log_through_2026-11-04.md`. Dated findings are indexed in
`docs/AUDIT_FINDINGS_REGISTER.md`. Do not append new session logs here; update the register, the guide and
`.ai-continuity/CURRENT.md` instead.

## Where things stand

- **The paper uses the frozen pipeline.** `code/R/run_all.R` → `data/sumstats_final.csv`, `tables/`, `plots/`.
  The live draft (`MPA_Kelp_MS_V7_working.docx`, private manuscript repo) reports 143 histories from 21 MPAs.
- **The current analysis build is the integrated review**, not yet adopted for the paper:
  `python3 code/integrated_review/run.py` (policy: `config/integrated_review_policy.json`) and the wrapper
  `python3 code/tools/run_whole_dataset_biomass_review.py`. It has 175 retained candidate histories → 162 fitted
  (62 before–after, 100 post-only; 13 with too little post-protection coverage). All 31 researcher/ecological
  exclusions are kept. Results: `outputs/integrated_review_20261103/`, `outputs/whole_dataset_biomass_update_20261105/`.
- **What is current:** `docs/ANALYSIS_GUIDE.md` (rules and reasons), `outputs/dataset_inventory_20261101/index.html`
  (decision tables), `outputs/methods_results_reference_20261104/` (response- and site-level decision rows for writing methods),
  `outputs/README.md` (status of every output folder), `code/README.md` (script index).

## Holds

- Do not promote integrated-review results to the paper until Adrian reviews the figures and settles reference-site
  choices. Then rerun pathway/resilience modules, regenerate paper figures/tables and the concordance gate (module 24).
- Do not choose site omissions, tolerances or cohorts by the direction or significance of results.
- KFM raw-data redistribution needs NPS permission. Do not publish KFM raw records (`data/source_archive/` in the data repo).
- No emails, pushes, submissions or publishing without an explicit request in the current conversation.

## Decisions that stay in force

- Researcher-advised and ecological exclusions are retained; baseline-trend/result-direction exclusions are removed (Nov 1).
- Calendar time: designation year counts as after; three matched strictly earlier years for before–after eligibility
  (four as sensitivity); all observed follow-up years; post-protection trends reported separately (Nov 3 policy).
- Shared recorded depth zones; all valid survey months with repeat visits averaged; density retained where biomass lacks sizes.
- KFM urchin/kelp sizes: same site, within seven days (exact-date sensitivity). PISCO size borrowing is restricted to
  the reviewed reference-only cases; report direct-size-only results alongside (R15).
- Accepted sources: Emily's combined VRG size file; original 140 Pelican Bay 1995 red urchins (135 as sensitivity);
  original LTER SBC 19.33 sizes (group means labelled as such).
- Frozen-pipeline inputs: SBC 127.5 model-1 monthly kelp calibration; Hamilton & Caselle (2015) sheephead
  `W(g) = 0.0144 * TL(cm)^3.04`; Naples LTER sheephead biomass retained (Oct 12). Resample with
  `a[sample.int(length(a), n, replace = TRUE)]`, never `sample(a, n)` on a possible singleton.

## Repository boundary

Analysis code (`code/`), tracked analysis-ready inputs (`data/harmonized/`, `data/cache/`), small comparator inputs
(`data/external/`), generated tables, figures and analysis documentation. This repository is **public**. Keep Word
manuscripts, cover letters, correspondence, publisher PDFs and restricted raw records out. Raw contributor data are
processed in `sbc-donham-kelp-mpa-data`; the manuscript lives in the private `sbc-donham-kelp-mpa-manuscript`.

## Frozen pipeline map

`run_all.R` runs the core pBACIPS/meta-analysis workflow, in-pipeline resilience modules and the concordance gate
`24_resilience_pipeline_check.R`. `run_figures_only.R` uses `data/cache/figures_snapshot.rds`. `run_resilience.R`
runs the broader suite (17 and 20 need optional raw PISCO inputs and skip cleanly without them); membership is in
`code/R/resilience_modules.R`.

- `00_libraries.R`, `00b_color_palette.R`, `00c_analysis_constants.R`, `01_utils.R`, `02_pBACIPS_function.R`
- `03_load_harmonized_data.R`; checks `33_kelp_calibration_check.R`, `35_sheephead_allometry_check.R`
- `08_effect_sizes.R`, `09_meta_analysis.R`, `10_temporal_analysis.R`; uncertainty `36`, `37`
- `11_figures.R`, `12_results_summary.R`, `13_additional_analyses.R`
- Resilience: `14` AR1 heatwave, `15` Kumagai comparison, `16` environmental moderators, `17` Eisaguirre
  reproduction, `18` effectiveness prediction, `19` two-heatwave repeatability, `20` SSWD, `21` stability,
  `22` resistance/recovery (main Figure 5), `23` ecological memory
- `53_kfm_site_continuity.R`, `57` size sensitivity, `58` calendar-window validation

Prose-to-output maps: `docs/RESILIENCE_SYNTHESIS.md`, `docs/RESILIENCE_PIPELINE_INTEGRATION.md`, `docs/CONCORDANCE_AUDIT.md`.
If a manuscript claim changes, update the table/figure, the prose-facing docs and module 24 together.

## Hygiene

- Use `here::here()` or repo-relative paths; optional inputs via `KUMAGAI_SUBTIDAL_CSV`, `KUMAGAI_CS_GRID`,
  `KUMAGAI_HUMAN_GRAVITY`, `DONHAM_DATA_REPO_PISCO`, `SBC_KELP_ENV_NETCDF`. No user-specific absolute paths.
- Dated report folders keep their paths; code and hash checks read them. Record status in `outputs/README.md`
  instead of moving them.
- New one-off scripts go in `code/tools/` and get a line in `code/README.md`. Scratch work goes in `tmp/` (ignored).

<!-- shared-ai-workspace:begin -->
## Shared Codex / Claude continuity

Updated: 2026-10-30.

- Read both `AGENTS.md` and `CLAUDE.md` in the current workspace once per session, when present. They are complementary instructions for the same working files; do not repeatedly follow cross-links, and read a shared symlink target only once.
- Read `/Users/adrianstier/AI-Workspace/WORKFLOW.md`. Resolve this workspace in `/Users/adrianstier/AI-Workspace/projects.json`: prefer an explicit project named by Adrian, then an exact `canonical_path` match, then companion `canonical_roots`. Read only the matching handoff and relevant source documents automatically; do not ask Adrian to paste context already available locally.
- Preserve repo-local conventions and existing authoritative source pointers. Newer explicit user decisions and dated primary evidence supersede stale summaries; neither filename has blanket precedence. If a substantive conflict remains unresolved, surface it before dependent work.
- Before editing, inspect current changes and relevant status. Keep one active editor per working tree. For multi-folder projects, verify access and read applicable instructions in each folder before working there; do not confuse folder access with reviewed content.
- Automatically follow `/Users/adrianstier/AI-Workspace/ACTIVITY-WORKFLOW.md`: at session start run `python3 /Users/adrianstier/AI-Workspace/activity.py start --project "EXACT PROJECT NAME" --tool codex` (Claude uses `--tool claude`), retain the session ID, and read each mapped root's `.ai-continuity/CURRENT.md` and recent activity. Before finishing or switching, update authoritative status/handoff, then run `finish` with truthful summary, next steps, tests, holds and source pointers. Do this without asking Adrian to remind you. If interrupted or only mechanically testing, use `abandon` without claiming review. Respect explicit read-only requests that prohibit record writes; report that limitation. Do not bypass an overlapping active session.
- After substantive work or before switching tools, update the existing project status and concise shared handoff with dated verified changes, tests, holds, blockers and next steps. Update the authoritative instruction source when a durable rule changes; avoid maintaining duplicate project facts in both files. Do not overwrite unrelated instructions or create recursive imports.
- Shared files provide continuity; chat history, sign-ins, schedules and app permissions do not synchronize automatically. Bypass permissions is Adrian's requested execution preference, not authorization to send messages, publish, submit, or execute a backlog without a request. Preserve all existing sending and privacy holds.
<!-- shared-ai-workspace:end -->
