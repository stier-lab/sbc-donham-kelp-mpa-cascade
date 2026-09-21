# Orientation — Donham-Stier-CA-MPA-2026 (analysis repo)

_Created: 2026-09-21 · Last updated: 2026-09-21 (rev 2) · branch `main` · last analysis commit 2026-09-04 (`f3b5878`)_

## What this is

The public analysis/code repo for the Donham & Stier California MPA kelp-forest paper.
It uses **progressive-change BACI-Pairs (pBACIPS)** + multilevel meta-analysis on 34
southern-California MPA/reference pairs to ask whether the CA MPA network restored the
kelp-forest trophic cascade. Headline as the live draft states it: MPAs held **higher giant-kelp
biomass and stronger kelp resistance/recovery through the 2014–2016 marine heatwave**, but
trophic-cascade recovery is only **partial** — sheephead biomass (+219%) and lobster density
(+196%) are up, purple-urchin density is down (−77%), yet most predator→urchin and urchin→kelp
cross-taxa links are weak (all six Table 3 meta-regressions n.s.).

**This is one of three sibling repos** (see *Related*). Code lives here; raw/harmonized data
in `../Donham-Stier-CA-MPA-Data-2026`; the Word manuscript + SI + submission machinery in
`../Donham-Stier-CA-MPA-MS-2026`. `CLAUDE.md` defines the boundary: **no manuscripts, cover
letters, or publisher PDFs in this repo.**

## Data

| Dataset | Location | Range / sites | Notes |
|---|---|---|---|
| `harmonized_response_ratios.csv` | `data/harmonized/` (tracked, 535K) | 1982–2023, 5 focal taxa | 3,341 lnRR rows; PISCO + KFM/MBON + SBC LTER |
| `harmonized_raw_responses.csv` | `data/harmonized/` (596K) | 1982–2023 | 6,280 rows, in/out density & biomass |
| `harmonized_landsat_rr.csv` | `data/harmonized/` (119K) | 1984–2021 | Landsat kelp canopy (Bell et al.) |
| `harmonized_site_metadata.csv` | `data/harmonized/` | 34 S. CA MPAs | Channel Is., SB, Catalina, LA/Orange, San Diego |
| Bootstrap caches | `data/cache/*.rds` (14 tracked, ~8 MB) | — | Lets `run_figures_only.R` work from a clean clone |
| Kumagai 2024 comparators | `data/external/kumagai2024/` | — | Vendored; used by scripts 15/16/18/24 |
| Raw monitoring data (~1.3 GB) | **Not here** — sibling data repo, symlinked to Google Drive | — | Gitignored: `data/{PISCO,MBON,LTER,LANDSAT,MPA}` |

Verified 2026-09-21: the four harmonized CSVs are **byte-identical** to the data repo's
`output/harmonized/` and match `checksums.md5`. The two repos are in sync.

Taxa: *Macrocystis pyrifera*, *Strongylocentrotus purpuratus*, *Mesocentrotus franciscanus*,
*Panulirus interruptus*, *Semicossyphus pulcher*.

## Pipeline (run order)

`00_libraries` → `00b_color_palette` → `00c_analysis_constants` → `01_utils` →
`02_pBACIPS_function` → `03_load_harmonized_data` → `08_effect_sizes` → `09_meta_analysis` →
`10_temporal_analysis` → `11_figures` → `13_additional_analyses` →
resilience modules `14/15/16/19/21/22/23` → `24_resilience_pipeline_check` → `12_results_summary`

- **Entry point:** `source(here::here("code","R","run_all.R"))` — ~3 min, logs to `logs/`
- **Figures only:** `run_figures_only.R` (~1 min, uses `data/cache/figures_snapshot.rds`)
- **Full resilience suite:** `run_resilience.R` — adds `17` (Eisaguirre reproduction) and
  `20` (SSWD/sunflower star); both need optional raw PISCO and skip cleanly without it
- **Headline outputs:** `tables/table_02_meta_analysis.csv`, `plots/fig_03_mean_effects.pdf`,
  and Figure 5 = `plots/fig_kelp_resilience.pdf` (script `22`)
- **Gate:** `24_resilience_pipeline_check.R` fails the pipeline if manuscript/SI resilience
  values drift from generated output. Module membership is defined in `resilience_modules.R`.

## Current state

- **Done.** Pipeline runs clean end-to-end (last full run 2026-09-04 06:48, 3.02 min, all
  stages OK). Reviewer-readiness rerun verified from a fresh public clone. Repo boundary
  cleanup complete (`docs/REPO_READINESS_AUDIT.md`). Claim ledger written
  (`docs/CONCORDANCE_AUDIT.md`). Table 3 corrected to Knapp-Hartung throughout.
- **In progress (in the MS repo).** `MPA_Kelp_MS_V7_working.docx` is the live draft,
  pre-submission to *J. Applied Ecology*, 233/350-word Summary, 5,603/7,000 words, 5 figs,
  3 tables, 11 authors (Donham first, Stier last; Emily corresponding).
  `make jae-submit-check` does **not** yet pass under strict mode.
- **Not blocked on Emily** (corrected 2026-09-21 from the actual mail thread; the COO said
  otherwise for three weeks). She answered all four 2026-08-31 checks the same day —
  South La Jolla and Carrington Pt are **in**, Fig 1 shows only reported sites, the public
  archive is **derived data only**, and the lobster-coefficient citation reading is right.
- **The live question is authorship.** On 2026-09-17 Emily returned an edited draft and
  said she is unsure she wants to be first author. Unanswered. Nothing should go to
  coauthors until it is settled. Full record (project content only) is in the MS repo at
  `correspondence/EMILY_DONHAM_CORRESPONDENCE.md`.
- **Emily's 2026-09-17 draft is now in the MS repo** at
  `manuscript/from_emily/2026-09-17_MPA_Kelp_MS_V7_EMDedits.docx` (saved 2026-09-21), but
  is not yet reconciled into the working copy.

## Open threads (resume here)

- [ ] **Answer Emily's 2026-09-17 authorship email.** This is the only thing genuinely
      holding the paper, and the decision is hers.
- [ ] **Reconcile Emily's 2026-09-17 draft** into `MPA_Kelp_MS_V7_working.docx`. Saved
      2026-09-21 to `manuscript/from_emily/`; delta is edit-level, not structural.
- [ ] **`edonham@ucsb.edu` is dead** (since 2026-09-02) but is still the corresponding
      address in V7, and `front_matter_audit.md` marks it PASS. Use
      `emily.m.donham@gmail.com`, and make that gate catch dead addresses.
- [x] ~~Stale README front matter~~ — fixed 2026-09-21: title and 11-author list match V7.
- [x] ~~`docs/ANALYSIS_REVISIONS.md` superseded numbers~~ — banner strengthened 2026-09-21 to
      name the drift (memo +185%/+79%, current +219%/+109%) and point at
      `docs/RESULTS_SUMMARY.md` as the only current source.
- [ ] **Declarations unfilled** (MS repo): Funding, CRediT, COI, ethics/permits, Data
      Availability; cover-letter addressee; 3–5 suggested referees.
- [ ] **Archive decision.** No license, no DOI, no Dryad upload. `docs/REPO_READINESS_AUDIT.md`
      flags the residual risk: public git *history* still contains removed manuscript/PDF
      objects — current-head cleanup does not erase them.
- [ ] **Data-repo loose ends:** KFM urchin size-frequency not ingested (bootstrap borrows PISCO
      sizes); two "CORRECTION NEEDED" items in `METHODS_REFERENCES.md`; `MACRO_AVE_SLOPE`
      cross-check against Rassweiler 2018 still open.
- [x] ~~COO reconciliation~~ — settled 2026-09-21: `donham-mpa` and `mpa-trophic-cascade`
      are **two distinct papers**, not duplicates. This one is the statewide 34-pair pBACIPS
      meta-analysis (11 authors, JAE); the other is SBC-LTE-only paired SEM + thermal
      sensitivity, Adrian-led, no venue. `donham-mpa` is now P0/internal-review and the
      spurious `blocked_on: emily-donham` is cleared from both records.

## Stack & reproduction

R ≥ 4.0 with `renv.lock` (542 KB, pinned). No `targets`. From a clean clone:
`renv::restore()` then `source(here::here("code","R","run_all.R"))` — no raw data needed,
harmonized CSVs and caches are tracked. Optional external inputs are overridden by env vars
only (`KUMAGAI_SUBTIDAL_CSV`, `KUMAGAI_CS_GRID`, `KUMAGAI_HUMAN_GRAVITY`,
`DONHAM_DATA_REPO_PISCO`, `SBC_KELP_ENV_NETCDF`) — never hard-code absolute paths.

## Related

- Data repo: `../Donham-Stier-CA-MPA-Data-2026` → `stier-lab/Donham-Stier-CA-MPA-Data-2026`
  (see its `DATA_SOURCES.md` and `output/harmonized/DATA_DICTIONARY.md`)
- Manuscript repo: `../Donham-Stier-CA-MPA-MS-2026` (private; `make` targets `coauthor-check`,
  `jae-readiness-audit`, `jae-submit-check`)
- GitHub (this repo): `stier-lab/sbc-2026-donham-kelp-mpa-cascade`
- Catalog: `~/repo-catalog/repos/stier-lab__Donham-Stier-CA-MPA-Data-2026.md`; COO: `coo show donham-mpa`
- **Not the same project:** `~/mpa-trophic-cascade-sbc` is SBC-LTE-only paired SEM (reconciled 2026-09-21).
- Correspondence: `../Donham-Stier-CA-MPA-MS-2026/correspondence/EMILY_DONHAM_CORRESPONDENCE.md`
- In-repo docs: `docs/RESULTS_SUMMARY.md`, `docs/RESILIENCE_SYNTHESIS.md`,
  `docs/CONCORDANCE_AUDIT.md`, `docs/REPO_READINESS_AUDIT.md`
