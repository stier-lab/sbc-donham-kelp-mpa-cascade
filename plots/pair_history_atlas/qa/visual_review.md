# Pair-history atlas QA

Reviewed 2026-10-09. This is a structured self-review, not an independent second-model critique.

## Coverage and provenance

- All 144 final-primary MPA/abbreviated-taxon/source/response keys matched; all 21 MPAs are represented. No unmatched final series or duplicate series-year records were found.
- The atlas contains 2805 paired annual observations, spanning 1982-2023. Every series' paired count matches its final effect-size N.
- The source/response counts are KFM biomass 12, KFM density 17, LTER biomass 6, LTER density 4, PISCO biomass 44, PISCO density 45, and Landsat biomass 16.
- The 51 harmonized keys absent from the final primary set were excluded. No rejected series were added.
- The manifest has 144 rows, four per page, with an explicit MPA column, source, taxon, response, dataset and series ranges, and n_before/n_after.
- Only the two harmonized annual tables, harmonized site metadata and sumstats_final.csv were read. Their checksums were unchanged during rendering and verification.

## Visual review

All 36 final pages were inspected through six ordered contact sheets. Full-resolution representative pages covered every source/response combination, long MPA labels and gap-heavy histories. Page 027 was also inspected after rasterizing the combined vector PDF at 300 dpi.

| Criterion | Score / 5 | Evidence |
| --- | ---: | --- |
| Proportions | 5 | Four stacked facets at 170 x 200 mm, with consistent panel dimensions. |
| Margins | 5 | Strip labels, axis labels and right-hand direct labels remain inside the export. |
| Alignment | 5 | Shared calendar-year and normalized-response scales across all pages. |
| Labels | 4 | Full MPA/source/taxon/response identification; borderless white direct-label backing keeps leaders from obscuring text. |
| Readability | 5 | No clipping or overlapping labels remained in the final review. |
| Spacing | 5 | Four readable facets per page without figure titles or subtitles. |
| Accessibility | 4 | Blue solid/circles and grey dashed/triangles remain distinguishable in deutan, protan and greyscale checks. |
| Data ink | 4 | Light event shading and designation lines provide context without masking trajectories. |

The final score is 37/40. The compact direct-label gutter is the principal visual compromise. Endpoint leaders are annotations, not continuations of annual observations.

## Reproduction and export checks

- Standalone Rscript rendering completed successfully.
- Sourcing from outside the repository preserved sentinel pipeline objects, the working directory, attached packages, options and original RNG state, including an initially absent .Random.seed.
- The combined PDF has 36 pages; all 36 individual PDFs and PNGs exist. Arial regular and bold fonts are embedded in the combined PDF.
- All page PNGs are 2007 x 2362 pixels at 300 dpi, with an opaque 8-bit RGB background.
- Missing calendar years break lines; isolated observations remain markers. All selected MPA/reference values are finite.
- Designation lines use actual metadata MPA_Start. Shading covers the complete 1997-1998 and 2014-2016 annual windows.

## Limitations and integration

The supplied annual responses are normalized proportions, not physical biomass or density units. Values are unchanged and are not averaged across sources. The common 0-1.25 vertical scale retains supplied values above 1. Annual uncertainty estimates are unavailable and were not inferred.

There are 97 histories with internal missing years, totaling 256 absent calendar years, and 52 histories without Before-period observations. Neither gaps nor baselines were imputed.

For Figure S15, include each distinct manifest png once in ascending page order. The six contact sheets, accessibility variants and PDF rasterization under qa/ are review artifacts, not supplementary figure pages. The legend and input hashes are in docs/pair_history_atlas_legend.md.

Changed paths are limited to code/R/29_pair_history_atlas.R, plots/pair_history_atlas/*, plots/fig_s_pair_history_atlas.pdf, tables/table_s_pair_history_atlas_index.csv and docs/pair_history_atlas_legend.md. No manuscript, SI, orchestrator, existing script or data was edited; no commit or push was made.
