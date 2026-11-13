# From the paper's Table 2 to the joint model: a step-by-step comparison

Written 2026-11-12 (evening) before any of the numbers below were computed. Descriptive; there are no pass rules and
nothing is adopted. The paper, its tables and the frozen pipeline are read, never changed.

## Why

The paper and the joint model give different numbers for the same species (density, log scale: lobster +1.05 against
+1.75; purple urchin −1.56 against −0.91; sheephead −0.25 against +0.20). They differ in three ways at once: what is
being estimated, which surveys go in, and how the estimate is made. This comparison changes one of those at a time.

## What the paper's number is

Table 2 pools one effect per history (a reserve, a species, a response, a survey program). The histories do not all
estimate the same thing (`data/sumstats_final.csv`, 143 histories):

- before-and-after histories (40): a fitted step, line or curve; the effect is a change relative to the before years;
- post-only histories fitted with a mean (61 "Mean" models in all): the effect is the average log ratio of protected to
  reference after designation, which includes any difference that was there before protection (siting);
- post-only histories fitted with a line: the change in the log ratio from the first to the last survey in the window.

The joint model's number is one thing for every reserve: the change in the protected sites relative to their own
baseline and their reference sites, ten years after designation, with siting differences estimated separately.

## Steps

All pooled values are per species and response. Each step reports the number of histories, the plain average, the
median, and a random-effects average (DerSimonian-Laird, with each history's own standard error). The random-effects
average is a simple stand-in for the paper's multilevel model; step 1 shows how close it comes.

| Step | Data | What is estimated | Tells us |
|---|---|---|---|
| 0 | Paper | Table 2 as published in the repository | starting point |
| 1 | Paper's per-history effects | the same effects, pooled simply; also split by model type (mean, line, step) and by program | whether simple pooling reproduces Table 2; how much each kind of history contributes |
| 2a | Paper's annual log ratios | inside against outside after designation (mean log ratio, after years) | the "siting included" quantity, made uniform |
| 2b | Paper's annual log ratios | after minus before (histories with at least 3 before years) | the before-and-after change, made uniform |
| 2c | Paper's annual log ratios | slope over the after years × 10 | the post-designation trend, made uniform |
| 3a–c | Joint model's survey data, turned into annual log ratios | the same three quantities | effect of the different survey inclusion (stations, years, months, programs) with the method held fixed |
| 4 | Joint model | siting difference; effect at ten years; siting + average effect over the after years | the model's version of 2a/3a (last item), 2b/3b (effect) |
| 5 | Joint model, before-and-after reserves only and post-only reserves only | effect at ten years | already run (round 2); included for completeness |

For step 3 the joint model's unit counts are turned into ratios in the plainest way: density per station and year
(count ÷ area searched), the unweighted mean over protected stations and over reference stations of a reserve in a year,
and log((protected + c) ÷ (reference + c)) with c = half the smallest positive station-year density of that species.
Years need both roles observed. Histories are kept separate by program, as in the paper. Sensitivity: c doubled and halved.

Density first. Biomass is compared at steps 0–2 from the paper's own ratios; at step 4 the joint model's biomass is
density plus mean mass. A biomass version of step 3 needs unit-level masses and is left for later.

Kelp in the paper pools diver biomass and Landsat canopy (14 of 31 histories are Landsat). They are shown separately.

## What will not be concluded

No step is "right". The comparison says which change moves each number. It does not say the paper's numbers are wrong,
and it does not make the joint model's numbers the new result.

## Change after the first run (2026-11-12, 21:30 PST)

Step 4's third item, "siting + average effect over the after years", was dropped. The model's siting term is an average
of logs across stations, and steps 2a/3a are logs of averages; with many low-count stations the two are not on the same
scale, so the composite invited a comparison it could not support. The siting term is still listed, labelled as such.
A proper model-based counterpart of 2a/3a needs posterior draws of expected densities per reserve and role; not built.
