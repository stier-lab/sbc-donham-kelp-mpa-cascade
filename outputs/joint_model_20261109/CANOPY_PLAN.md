# Landsat kelp canopy in the joint model: design, inference limits and pass rules

Written 2026-11-12, before any full canopy fit was run on GRIT. Review only; nothing is adopted. The paper, its figures
and the frozen pipeline are untouched. Seen before these rules were written: one short fit of the spring model to the real
data (2 chains, not converged; E10 +0.49 [−0.24, +1.15]) and one simulated data set. Neither is evidence for the rules.

## 1. Why the canopy is its own response

Reed, Rassweiler, Nelson and Kui (2026, *Front. Mar. Sci.* 13:1845447; monthly dives at three SBC LTER sites, 2003–2025,
5–10 m) show what a surface canopy can stand for, even when it is measured by divers without sensor error:

- the surface canopy is 12–38% of standing biomass;
- it explains 66–81% of the variance in subsurface biomass on average, 88–95% in February–May and 29–45% in
  September–November;
- it explains 19–35% of the variance in plant density on average, and 0–14% from June to November;
- much kelp can stand under little or no canopy; a large canopy implies much kelp below, not the reverse;
- their plots stop at 10 m and forests reach 20 m; they suggest deeper plants spend less time in the canopy.

Our diver kelp data are plants (KFM, PISCO) and stipes (LTER), mostly surveyed June–October. Landsat sees the canopy.
These are different quantities in different seasons, so the canopy is never pooled with diver kelp as one measurement.

## 2. Data

`prep_canopy.py` rebuilds quarterly canopy biomass for each MPA and reference polygon from the pixel file
(`CAkelpCanopyEnv_2021.nc`, 332,640 pixels × 152 quarters, 1984–2021) with Tom Bell's corrected masks (delivery of
2026-10-21). With his coverage rules it reproduces the delivered annual-maximum table exactly for all 58 non-empty
polygons (largest difference 0). The cohort is the 18 comparisons frozen on 2026-11-01
(`outputs/satellite_primary_20261101/CRITERIA.md`), with the designation years in its cohort file. Nothing is added or dropped.

Facts from the rebuilt table (all 58 polygons): canopy is absent in 16–23% of covered polygon-quarters; the annual
maximum falls in January–March 16% of the time, April–June 41%, July–September 25%, October–December 18%. So the
annual maximum mixes seasons from year to year. Median pixel depth per polygon is about 10 m (range 5–17 m).

One observation is one polygon in one quarter with at least 75% of its pixels observed: canopy density = biomass ÷
(observed pixels × 900 m²), kg wet per m².

## 3. Model (tag `can` in `run_multi.py`)

Latent log canopy density for each polygon and year, with the same structure as the density model (polygon baseline,
siting difference of the protected polygon, region × year effects with shared and regional trends, polygon AR(1),
protection curve by years since designation with MPA-specific rates). Observation model:

- quarter offsets (sum to zero);
- canopy detected ~ Bernoulli, with probability rising with the latent level and with the number of observed pixels;
- log density when detected ~ Student-t (5 d.f.) around the latent level.

A quarter without canopy is "no canopy detected". It lowers the latent level by an amount the data decide and is never a
biomass of zero. The effect reported is E10, the log change in canopy density ten years after designation. MPAs
designated in 2012 have nine years of canopy data after designation, so for them E10 leans on the curve's last step.

Main specification, fixed here: spring quarters (January–June), all depths, designation year counted as after (the joint
model's convention; in spring the 2012 MPAs were already in force), flat baseline, random-walk curve prior with the
independent-increment prior always reported beside it.

Regions are the proposed three (northern islands, southern islands, mainland). Point Conception is placed with the Santa
Barbara mainland. Both are proposals for Adrian to confirm.

## 4. Pass rules

| | Check | Rule |
|---|---|---|
| C1 | Independent arithmetic | `test_multi.py` check 4 on the real table: observation term equals a separate SciPy calculation (relative difference < 1e-9) for three season definitions |
| C2 | Convergence | main fits meet the same thresholds as M1 |
| C3 | Calibration | 20 simulated data sets on the real design at real-data noise levels: 90% intervals cover the true E10 in at least 16; 20 with no effect: at most 4 intervals exclude zero. Mean error reported |
| C4 | Shuffled labels | 40 shuffles per prior, protected label dealt at random among each comparison's polygons. R1a and R1b as in `VALIDATION_PLAN_v2.md`; randomization p reported |
| C5 | Time placebo | years before designation only, fake designation in 1990, 1993, 1996, 1999: the 90% interval covers zero in at least 3 of 4, for each prior |
| C6 | Seeds | two seeds agree within 0.05 |
| C7 | Leave one MPA out | no single MPA reverses the sign of an effect whose main interval excludes zero; shifts above 0.25 listed |

Reported side by side, never chosen by result (C8): season (spring, annual maximum, July–December, all quarters); designation
year as after or before; pixels to 10 m or deeper than 10 m; flat baseline or pre-trend; both curve priors.

Not built yet, so not gated: a predictive check of the canopy observation model (share of quarters without canopy by
role and period, spread of log density). Shared pixels between neighbouring comparisons' reference areas are not
modelled; they make those comparisons less than independent.

## 5. What the canopy result may be used to say

- It is a statement about surface canopy biomass in spring over whole reserves and their reference areas.
- It is not a statement about the number of kelp plants, their size structure, or kelp in autumn, and it is weaker for
  reefs deeper than about 10 m. The depth split is reported for that reason.
- A canopy effect and a diver kelp effect that disagree are not a contradiction; they measure different things at
  different scales (a transect against a polygon of hundreds to thousands of 30 m pixels).
- The earlier satellite analysis (annual maximum, log ratio with an additive offset for zeros, designation year as
  before) is kept as it is. The comparison between it and this model on the same 18 comparisons is part of the report.

## 6. Linking canopy to diver kelp (second step, after C1–C7)

A link is added only after the canopy response stands on its own. Order of work:

1. Describe first: for each comparison and role, the fitted diver kelp series (stipe density and biomass, not plant
   counts) against the canopy series, by season, with the number of dive stations and their depths. No model.
2. If the series overlap enough (Adrian to judge from that table), add an edge canopy ← diver kelp to the joint fit with
   separate coefficients for year-to-year fluctuation and for the protection function, a slope that varies by
   comparison, and no assumption that low canopy means low kelp (the detection term already carries that asymmetry).
3. Report the canopy effect with and without the link.
