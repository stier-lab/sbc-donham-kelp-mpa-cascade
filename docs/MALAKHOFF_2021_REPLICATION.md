# Malakhoff and Miller replication status

Updated 16 October 2026. **Partial reconstruction, not a completed replication
with the newest data. None of these outputs is cleared for the manuscript.**

The best-available fixed-cohort extension is now complete as an exploratory
analysis: 75 additional fits compare 2019 with 2022 and partial 2023 using
three reserves and eight controls with measured before-after support. It
includes covariance, period-boundary and reserve-influence sensitivities,
before-period trends, observed trajectories and exported diagnostics. See
`docs/MALAKHOFF_2021_AVAILABLE_DATA_RESULTS.md` for estimates and interpretation.
Exact replication and the kelp endpoints remain incomplete for the reasons below.

The subsequent controlled comparison with the manuscript is documented in
`docs/MALAKHOFF_2021_METHOD_BRIDGE.md`. Its key finding is that the current
purple-density synthesis mixes nine fitted CI trajectory changes, three CI
mean contrasts and one BACIPS step; it is not directly comparable to a single
designation-focused BACI interaction. The bridge runs within this standalone
workflow and leaves the manuscript pipeline unchanged.

## Sources read

- Malakhoff and Miller (2021), *After 15 years, no evidence for trophic cascades
  in marine protected areas*, Proceedings B 288:20203061,
  https://doi.org/10.1098/rspb.2020.3061. All main-text sections and the
  supplementary captions and model table were read.
- Supplement: https://doi.org/10.6084/m9.figshare.13709928.v1. Downloaded file
  MD5 matches the public record: `e082519c91d2976381166760a720b127`.
- Author code: https://github.com/kmalakhoff/urchin-baci at
  `f892f7ace06e680c60175b1428ee8fc2463f4c71`. The repository contains
  `Data Prep.R` and `BACI paper script.R`, but no analysis CSVs.
- New NPS package: https://irma.nps.gov/DataStore/Reference/Profile/2318436,
  DOI `10.57830/2318436`, covers 1982-2025. The retrieved public profile marks
  files `Downloadability: Internal`, `HasCUI: true`, and
  `CanDownloadFiles: false`. Restricted files were not requested or downloaded.
  Only the public profile was saved.

## Design to preserve

The BACI cohort comprises KFM sites 1-11 and 14-16, with sites 2, 6, 9 and 14
as impacts and ten controls. Sites 12 and 13 belong to the older Anacapa reserve
and must be analysed separately with site 11. This is not the manuscript's
MPA-specific paired-reference design, and uses neither the pooled log-response
contrasts nor the current resilience cohort.

The published study window is 1991-2019 for BACI, with the older Anacapa series
discussed separately. Preserve wet urchin biomass and kelp stipe density as
distinct responses; do not substitute dry kelp biomass or plant density.

## Available input audit

The reconstruction reads the full local MBON count files before any manuscript
site filtering and uses only KFM sizes from `ALL_sizefreq_2024.csv`. It never
borrows PISCO or LTER sizes to fill missing KFM records. File hashes are recorded
in `outputs/malakhoff2021/input_manifest.csv`.

- Counts are observed counts divided by the row's sampled area. Missing values
  are not converted to zero. The quadrat protocol is filtered separately from
  predator bands and visual-fish transects.
- Size-frequency counts weight wet mass at each observed diameter, using
  0.00059 x diameter^2.917 for red and 0.00059 x diameter^2.870 for purple
  urchins. Diameters are in mm and biomass is wet grams per square metre.
- Sampling-row count keys and annual join keys are checked for duplicates.
  Repeated size bins are pooled by their recorded individual counts, not
  dropped by row order.
- Observed zero densities yield zero biomass in the corrected reconstruction.
  Positive densities without a same-site/year/species size distribution remain
  missing. No interpolation or cross-site size substitution is performed.
- Rodes Reef (site 5) and Fry's Harbor (site 7) have **no red size distributions**
  in the local file, despite positive counts. Their red biomass cannot be
  reconstructed. Rare observed purple zeros at these sites are not evidence of
  complete size coverage.
- Southeast Sea Lion Rookery (site 14), one of four impact sites, lacks sizes
  for 1991-2004. Its available post-period biomass does not provide a measured
  before-after contrast for that reserve.
- No 2020 counts occur in the selected local export. Only original sites 9,
  14 and 15 have 2023 counts, and site 14 lacks matching 2023 sizes. Treat 2023
  as partial coverage, not a network-wide update. Runs ending in 2022 are
  exported separately from those ending in 2023.
- The local kelp raw file supplies sampled plants' stipe counts. A confirmed
  sampled-area/effort table and the original annual stipe-density CSV are not
  available; a per-plant mean cannot be silently turned into stipes per area.

These limitations prevent both verification of published numerical results and
a complete newest-data update. In particular, the red models have only twelve
sites with any reconstructed biomass, not the published fourteen-site cohort.

## Code and paper differences to resolve

1. **Period boundary:** code labels years <=2004 as before, i.e. the first after
   year is 2005. The paper's before/after wording suggests 2004 is after.
   Both boundaries are tested explicitly without changing the manuscript.
2. **Biomass calculation:** older sizes are transformed individually and then
   averaged. The 2017-2019 prep code instead transforms mean diameter.
   Nonlinear allometry makes these unequal. Both constructions are retained;
   the code-era construction extends that later algorithm beyond 2019.
3. **Covariance:** several archived fits nest sites within year and use
   `corAR1()` without a time variable. A single-species annual year/site cell
   contains one observation, so this does not model between-year dependence
   within a site. The supplement's ST2 describes by-site random effects, adding
   another specification ambiguity. Four original-formula species fits were
   executed but tagged `not_temporal_AR1`; their P values are not approved.
4. **Small urchins:** code includes 35-mm individuals in the small class, while
   the paper specifies <35 mm. Code allocates biomass using count proportions,
   which is not a size-class-specific mass estimate. Its grouped denominator
   also needs verification: dropping the last species grouping under current
   dplyr semantics divides across species within a size-class/site/year cell.
   Code-cutoff, paper-cutoff, literal grouped-denominator and physically
   size-specific mass constructions are exported separately. These differences
   are not proof of what generated the published results without the original
   CSVs and software environment.
5. **Kelp labelling:** main Methods and supplement S3 use stipe density, whereas
   the S9 caption and a code axis call it biomass. The intended response must
   be resolved from the actual input file, not a caption alone.

## What ran

`Rscript code/replications/malakhoff2021/run.R` produces a standalone package
under `outputs/malakhoff2021/`, leaving all manuscript tables, figures, data and
models unchanged.

| Analysis | Current status |
|---|---|
| Red and purple BACI | Partial source reconstruction; original-formula runs and explicitly separate site/calendar AR1 sensitivities |
| Both-species BACI | Summed wet-biomass sensitivity only; not the paper's pooled-species random-slope model |
| Small-urchin BACI | Separate construction/boundary sensitivities, not verified reproduction of the published selected models |
| Size distributions | Annual proportions and mass constructions audited; all four-class published displays not yet reproduced |
| Kelp BACI and kelp-urchin associations | Blocked by missing verified stipe-density/effort inputs |
| Predator histories | Annual densities plotted descriptively from band and visual-fish counts; no new predator BACI inference |
| Anacapa | Separate exploratory ANCOVA for sites 11-13; available size support starts in 1985, not the paper's stated 1981 |
| Fishery-catch explorations | Not reproduced; original block-level catch inputs unavailable |

There are 78 site/calendar AR1 sensitivity fits and four original year/site
species fits. Every result has `eligible_for_manuscript = FALSE`, including
fits that converge. The site/calendar fits change the covariance structure and
are **not** exact original analyses. They omit the original year random effect
and do not resolve common year shocks, incomplete coverage, selection or
fishing-effort confounding. Raw nominal P values must not be cherry-picked from
the sensitivity grid. Residual/QQ diagnostic plots and residual lag-one
correlations are exported, but do not constitute a passed publication QA gate.

Available-data trajectories were rendered and inspected for clipping and
missing assets. They are internal coverage-dependent displays, not updated
versions of the published figures. With only one contributing site in a group,
SE is undefined and no ribbon is drawn.

Run `Rscript code/replications/malakhoff2021/test.R` after rebuilding. Tests
passed for the site cohort, weighted allometry, zero/missing distinctions,
size-class boundaries, join uniqueness and quarantine flags.

## Remaining requirements

1. Obtain the original `Biomass_1991_2019.csv`,
   `KFM_Raw_Biomass_Edited2.csv`, `Annual_Kelp_Mean_Stipe_Density.csv` and
   associated raw input/effort tables. Compare annual values before fitting.
2. Obtain authorised access to the full 1982-2025 KFM release, including counts,
   natural-habitat urchin sizes, kelp stipe surveys and effort metadata.
3. Resolve the original model choices against ST2, then match the published
   estimates, SEs, P values and AICs before extending those exact specifications.
4. Produce old-data/old-method, refreshed-old-window/old-method and
   updated-window/old-method comparisons. Keep statistical corrections in a
   fourth, clearly labelled sensitivity layer.

A question requesting local data or permission to contact Scott was presented
to Adrian. No request email, manuscript edit, Drive overwrite or push has been
made for this replication.
