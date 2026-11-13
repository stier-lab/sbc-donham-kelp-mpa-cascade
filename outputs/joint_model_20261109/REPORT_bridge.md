# From the paper's Table 2 to the joint model, one change at a time

Plan and limits: `BRIDGE_PLAN.md`, written before these numbers were computed. Descriptive; nothing is adopted and the paper is unchanged. All values are log ratios or log changes. "Pooled" is the paper's multilevel estimate at step 0, a simple random-effects average at steps 1-3, and the model's posterior mean at steps 4-5.

## Reading of the tables (written 2026-11-12 after the first run; interpretation, not a rule)

**1. The two data sets agree.** When the same plain quantity is computed from the paper's annual ratios (step 2) and from
the joint model's surveys (step 3), the answers are close for every species: lobster after-minus-before +1.54 against
+1.44; purple urchin −0.88 against −0.82; red urchin slope −0.70 against −0.87; sheephead inside-against-outside +0.44
against +0.55. The different survey rules (stations, months, years) are not what separates the paper from the joint model.

**2. What is being estimated is what separates them.** The paper's pooled number mixes three quantities, in different
proportions for each species (step 1). They differ a lot within a species:

| Species (density) | Inside against outside, after | After minus before | Slope over after years × 10 | Paper Table 2 | Joint model at ten years |
|---|---|---|---|---|---|
| Spiny lobster | +1.13 | +1.54 | +0.75 | +1.05 | +1.75 / +1.82 |
| Sheephead | +0.44 | +0.37 | −0.10 | −0.25 | +0.20 / +0.15 |
| Purple urchin | −0.95 | −0.88 | −1.19 | −1.56 | −0.91 / −0.28 |
| Red urchin | −0.02 | −0.24 | −0.70 | −0.46 | −0.63 / −0.27 |

(Columns 2–4 from the paper's own annual ratios; joint model under the two curve assumptions.)

- **Lobster.** The paper's +1.05 is mostly the inside-against-outside level (11 of 17 histories), which mixes in sites
  where little has changed. The before-and-after histories give +1.5 in both data sets. The joint model's +1.75 is an
  effect at year ten rather than an average over the first decade, so it sits above +1.5. The gap is the estimand.
- **Purple urchin.** The paper's −1.56 is mostly the fall of the ratio over the after years (9 of 13 histories, pooled
  −1.63). That trend is in both data sets. The joint model shrinks it to −0.91, or −0.28 when each step must earn its
  change, because the reserves with before surveys show about half as much.
- **Red urchin.** For density the paper and the joint model agree in sign and roughly in size. For biomass the paper is
  +0.57 and the joint model −0.22. Nine of the paper's 14 biomass histories are inside-against-outside levels (+0.62):
  protected sites hold larger red urchins. The trend after designation is downward (−0.39). The joint model reports the
  change, which is fewer animals (−0.63) partly offset by larger ones (+0.41).
- **Lobster biomass.** The paper has 5 histories with sizes (+0.42; PISCO −0.03, LTER +0.98). The joint model takes the
  density change from all 13 reserves (+1.75) and adds the size change from the sites that have sizes (+0.38).
- **Kelp.** The paper's 31 biomass histories are 14 Landsat (+0.75) and 17 diver (PISCO +1.00, KFM +2.14, LTER +2.08).
  The joint model's diver estimate (+0.96) and its separate canopy estimate (+0.49) are both lower than the paper's
  pooled +1.14, which is pulled up by the before-and-after KFM and LTER histories.

**3. One number in Table 2 does not follow from its own histories.** Sheephead density is −0.25 in Table 2. Of the 20
histories behind it, 15 are positive; their plain average is +0.39, their median +0.27, and both a fixed and a
random-effects average are positive (+0.23, +0.34). The other species' Table 2 values are within about 0.2 of the simple
pooling. The paper's model shares one random effect per reserve and one per program across all species and has no
separate term for differences among histories, so a few very precise histories and the other species at the same
reserves can pull a species' average. That is the likely cause; I have not rerun the frozen pipeline to confirm it.
This is about Table 2 as it stands in the repository and is worth a look whatever happens with the joint model.

## Spiny lobster, density

| Step | What | Histories or reserves | Plain average | Median | Pooled | SE |
|---|---|---|---|---|---|---|
| 0 | Paper Table 2 (multilevel model) | 17 |  |  | +1.05 | 0.46 |
| 1 | Paper's histories, pooled simply | 17 | +1.21 | +0.89 | +1.23 | 0.29 |
| 1 |   of which: mean after (siting included) | 11 | +0.99 | +0.61 | +1.01 | 0.39 |
| 1 |   of which: line over after years | 3 | +1.79 | +2.16 | +1.73 | 0.74 |
| 1 |   of which: before-and-after (step, line or curve) | 3 | +1.41 | +1.82 | +1.81 | 0.43 |
| 1 |   of which: KFM | 5 | +1.70 | +1.82 | +1.84 | 0.50 |
| 1 |   of which: LTER | 2 | +0.25 | +0.25 | +0.22 | 0.31 |
| 1 |   of which: PISCO | 10 | +1.15 | +1.36 | +1.15 | 0.39 |
| 2a | Paper's annual ratios: inside against outside, after years | 17 | +1.11 | +0.77 | +1.13 | 0.31 |
| 2b | Paper's annual ratios: after minus before | 6 | +1.49 | +1.59 | +1.54 | 0.41 |
| 2c | Paper's annual ratios: slope over after years x 10 | 17 | +0.70 | +0.63 | +0.75 | 0.24 |
| 3a | Joint model's surveys as ratios: inside against outside, after years | 18 | +1.36 | +0.75 | +1.37 | 0.33 |
| 3b | Joint model's surveys as ratios: after minus before | 6 | +1.37 | +1.61 | +1.44 | 0.67 |
| 3c | Joint model's surveys as ratios: slope over after years x 10 | 18 | +1.13 | +1.11 | +0.91 | 0.26 |
| 3a | Joint model's surveys as ratios: inside against outside, after years (offset doubled) | 18 | +1.22 | +0.64 | +1.23 | 0.29 |
| 3b | Joint model's surveys as ratios: after minus before (offset doubled) | 6 | +1.20 | +1.33 | +1.26 | 0.59 |
| 3c | Joint model's surveys as ratios: slope over after years x 10 (offset doubled) | 18 | +1.00 | +1.01 | +0.82 | 0.22 |
| 3a | Joint model's surveys as ratios: inside against outside, after years (offset halved) | 18 | +1.50 | +0.95 | +1.51 | 0.37 |
| 3b | Joint model's surveys as ratios: after minus before (offset halved) | 6 | +1.53 | +1.90 | +1.62 | 0.75 |
| 3c | Joint model's surveys as ratios: slope over after years x 10 (offset halved) | 18 | +1.26 | +1.20 | +0.99 | 0.29 |
| 4 | Joint model (trend carries forward): siting difference at designation (average over reserves, log scale per station; not a ratio of means) | 13 |  |  | -0.32 | 0.56 |
| 4 | Joint model (trend carries forward): effect at ten years | 13 |  |  | +1.75 | 0.42 |
| 5 | Joint model (trend carries forward): effect at ten years, before-and-after reserves only | 6 |  |  | +1.72 | 0.48 |
| 5 | Joint model (trend carries forward): effect at ten years, post-only reserves only | 7 |  |  | +1.11 | 0.61 |
| 4 | Joint model (each step near zero): siting difference at designation (average over reserves, log scale per station; not a ratio of means) | 13 |  |  | -0.40 | 0.56 |
| 4 | Joint model (each step near zero): effect at ten years | 13 |  |  | +1.82 | 0.44 |
| 5 | Joint model (each step near zero): effect at ten years, before-and-after reserves only | 6 |  |  | +1.82 | 0.53 |
| 5 | Joint model (each step near zero): effect at ten years, post-only reserves only | 7 |  |  | +1.26 | 0.67 |

## Spiny lobster, biomass

| Step | What | Histories or reserves | Plain average | Median | Pooled | SE |
|---|---|---|---|---|---|---|
| 0 | Paper Table 2 (multilevel model) | 5 |  |  | +0.42 | 0.31 |
| 1 | Paper's histories, pooled simply | 5 | +0.46 | +0.56 | +0.42 | 0.39 |
| 1 |   of which: mean after (siting included) | 4 | +0.19 | +0.48 | +0.22 | 0.42 |
| 1 |   of which: line over after years | 1 | +1.51 | +1.51 |  |  |
| 1 |   of which: LTER | 2 | +1.18 | +1.18 | +0.98 | 0.28 |
| 1 |   of which: PISCO | 3 | -0.03 | +0.39 | -0.03 | 0.53 |
| 2a | Paper's annual ratios: inside against outside, after years | 5 | +0.24 | +0.33 | +0.32 | 0.30 |
| 2c | Paper's annual ratios: slope over after years x 10 | 5 | +1.50 | +1.31 | +0.87 | 0.39 |
| 4 | Joint model (trend carries forward): effect at ten years | 13 |  |  | +2.13 | 0.44 |
| 4 | Joint model (each step near zero): effect at ten years | 13 |  |  | +2.05 | 0.47 |
| 4 | Joint model (trend carries forward): effect at ten years, mean mass | 13 |  |  | +0.38 | 0.15 |
| 4 | Joint model (each step near zero): effect at ten years, mean mass | 13 |  |  | +0.23 | 0.16 |

## Sheephead, density

| Step | What | Histories or reserves | Plain average | Median | Pooled | SE |
|---|---|---|---|---|---|---|
| 0 | Paper Table 2 (multilevel model) | 20 |  |  | -0.25 | 0.46 |
| 1 | Paper's histories, pooled simply | 20 | +0.39 | +0.27 | +0.34 | 0.11 |
| 1 |   of which: mean after (siting included) | 9 | +0.57 | +0.23 | +0.39 | 0.13 |
| 1 |   of which: line over after years | 5 | +0.21 | +0.31 | +0.18 | 0.38 |
| 1 |   of which: before-and-after (step, line or curve) | 6 | +0.27 | +0.40 | +0.31 | 0.36 |
| 1 |   of which: KFM | 4 | +0.50 | +0.74 | +0.52 | 0.46 |
| 1 |   of which: LTER | 1 | -0.42 | -0.42 |  |  |
| 1 |   of which: PISCO | 15 | +0.42 | +0.31 | +0.32 | 0.12 |
| 2a | Paper's annual ratios: inside against outside, after years | 20 | +0.42 | +0.16 | +0.44 | 0.14 |
| 2b | Paper's annual ratios: after minus before | 7 | +0.36 | +0.20 | +0.37 | 0.26 |
| 2c | Paper's annual ratios: slope over after years x 10 | 20 | -0.47 | -0.24 | -0.10 | 0.13 |
| 3a | Joint model's surveys as ratios: inside against outside, after years | 23 | +0.56 | +0.19 | +0.55 | 0.16 |
| 3b | Joint model's surveys as ratios: after minus before | 8 | +0.14 | +0.03 | +0.20 | 0.29 |
| 3c | Joint model's surveys as ratios: slope over after years x 10 | 23 | -0.83 | -0.42 | -0.21 | 0.16 |
| 3a | Joint model's surveys as ratios: inside against outside, after years (offset doubled) | 23 | +0.52 | +0.20 | +0.52 | 0.16 |
| 3b | Joint model's surveys as ratios: after minus before (offset doubled) | 8 | +0.08 | -0.00 | +0.13 | 0.22 |
| 3c | Joint model's surveys as ratios: slope over after years x 10 (offset doubled) | 23 | -0.62 | -0.32 | -0.16 | 0.14 |
| 3a | Joint model's surveys as ratios: inside against outside, after years (offset halved) | 23 | +0.58 | +0.19 | +0.57 | 0.16 |
| 3b | Joint model's surveys as ratios: after minus before (offset halved) | 8 | +0.21 | +0.07 | +0.28 | 0.34 |
| 3c | Joint model's surveys as ratios: slope over after years x 10 (offset halved) | 23 | -1.06 | -0.56 | -0.25 | 0.17 |
| 4 | Joint model (trend carries forward): siting difference at designation (average over reserves, log scale per station; not a ratio of means) | 17 |  |  | +0.12 | 0.37 |
| 4 | Joint model (trend carries forward): effect at ten years | 17 |  |  | +0.20 | 0.18 |
| 5 | Joint model (trend carries forward): effect at ten years, before-and-after reserves only | 7 |  |  | +0.15 | 0.26 |
| 5 | Joint model (trend carries forward): effect at ten years, post-only reserves only | 10 |  |  | +0.19 | 0.28 |
| 4 | Joint model (each step near zero): siting difference at designation (average over reserves, log scale per station; not a ratio of means) | 17 |  |  | +0.14 | 0.37 |
| 4 | Joint model (each step near zero): effect at ten years | 17 |  |  | +0.15 | 0.15 |
| 5 | Joint model (each step near zero): effect at ten years, before-and-after reserves only | 7 |  |  | +0.09 | 0.17 |
| 5 | Joint model (each step near zero): effect at ten years, post-only reserves only | 10 |  |  | +0.10 | 0.19 |

## Sheephead, biomass

| Step | What | Histories or reserves | Plain average | Median | Pooled | SE |
|---|---|---|---|---|---|---|
| 0 | Paper Table 2 (multilevel model) | 17 |  |  | +0.93 | 0.26 |
| 1 | Paper's histories, pooled simply | 17 | +0.73 | +0.71 | +0.78 | 0.12 |
| 1 |   of which: mean after (siting included) | 9 | +1.28 | +0.87 | +0.90 | 0.16 |
| 1 |   of which: line over after years | 3 | +1.08 | +0.64 | +0.59 | 0.21 |
| 1 |   of which: before-and-after (step, line or curve) | 5 | -0.44 | +0.71 | +0.09 | 0.59 |
| 1 |   of which: LTER | 2 | -2.53 | -2.53 | -2.57 | 0.93 |
| 1 |   of which: PISCO | 15 | +1.17 | +0.87 | +0.81 | 0.11 |
| 2a | Paper's annual ratios: inside against outside, after years | 17 | +1.02 | +0.82 | +1.03 | 0.20 |
| 2b | Paper's annual ratios: after minus before | 6 | -0.02 | +0.29 | +0.18 | 0.34 |
| 2c | Paper's annual ratios: slope over after years x 10 | 17 | -0.56 | -0.06 | +0.16 | 0.19 |
| 4 | Joint model (trend carries forward): effect at ten years | 17 |  |  | +0.53 | 0.23 |
| 4 | Joint model (each step near zero): effect at ten years | 17 |  |  | +0.41 | 0.21 |
| 4 | Joint model (trend carries forward): effect at ten years, mean mass | 17 |  |  | +0.33 | 0.14 |
| 4 | Joint model (each step near zero): effect at ten years, mean mass | 17 |  |  | +0.26 | 0.14 |

## Purple urchin, density

| Step | What | Histories or reserves | Plain average | Median | Pooled | SE |
|---|---|---|---|---|---|---|
| 0 | Paper Table 2 (multilevel model) | 13 |  |  | -1.56 | 0.47 |
| 1 | Paper's histories, pooled simply | 13 | -1.52 | -1.46 | -1.35 | 0.19 |
| 1 |   of which: mean after (siting included) | 3 | -0.44 | -0.39 | -0.50 | 0.12 |
| 1 |   of which: line over after years | 9 | -1.83 | -1.70 | -1.63 | 0.21 |
| 1 |   of which: before-and-after (step, line or curve) | 1 | -1.99 | -1.99 |  |  |
| 1 |   of which: KFM | 3 | -1.97 | -1.99 | -1.76 | 0.33 |
| 1 |   of which: PISCO | 10 | -1.39 | -1.44 | -1.21 | 0.21 |
| 2a | Paper's annual ratios: inside against outside, after years | 13 | -0.96 | -0.62 | -0.95 | 0.27 |
| 2b | Paper's annual ratios: after minus before | 4 | -0.88 | -0.79 | -0.88 | 0.57 |
| 2c | Paper's annual ratios: slope over after years x 10 | 13 | -1.28 | -1.33 | -1.19 | 0.29 |
| 3a | Joint model's surveys as ratios: inside against outside, after years | 18 | -0.54 | -0.26 | -0.54 | 0.25 |
| 3b | Joint model's surveys as ratios: after minus before | 8 | -0.86 | -1.03 | -0.82 | 0.52 |
| 3c | Joint model's surveys as ratios: slope over after years x 10 | 18 | -1.45 | -1.69 | -1.48 | 0.25 |
| 3a | Joint model's surveys as ratios: inside against outside, after years (offset doubled) | 18 | -0.55 | -0.26 | -0.53 | 0.25 |
| 3b | Joint model's surveys as ratios: after minus before (offset doubled) | 8 | -0.84 | -0.96 | -0.80 | 0.51 |
| 3c | Joint model's surveys as ratios: slope over after years x 10 (offset doubled) | 18 | -1.43 | -1.69 | -1.45 | 0.25 |
| 3a | Joint model's surveys as ratios: inside against outside, after years (offset halved) | 18 | -0.54 | -0.26 | -0.55 | 0.25 |
| 3b | Joint model's surveys as ratios: after minus before (offset halved) | 8 | -0.87 | -1.10 | -0.83 | 0.53 |
| 3c | Joint model's surveys as ratios: slope over after years x 10 (offset halved) | 18 | -1.47 | -1.70 | -1.49 | 0.26 |
| 4 | Joint model (trend carries forward): siting difference at designation (average over reserves, log scale per station; not a ratio of means) | 13 |  |  | +0.95 | 0.51 |
| 4 | Joint model (trend carries forward): effect at ten years | 13 |  |  | -0.91 | 0.55 |
| 5 | Joint model (trend carries forward): effect at ten years, before-and-after reserves only | 7 |  |  | -0.41 | 0.65 |
| 5 | Joint model (trend carries forward): effect at ten years, post-only reserves only | 6 |  |  | -1.51 | 0.85 |
| 4 | Joint model (each step near zero): siting difference at designation (average over reserves, log scale per station; not a ratio of means) | 13 |  |  | +0.84 | 0.49 |
| 4 | Joint model (each step near zero): effect at ten years | 13 |  |  | -0.28 | 0.41 |
| 5 | Joint model (each step near zero): effect at ten years, before-and-after reserves only | 7 |  |  | -0.07 | 0.32 |
| 5 | Joint model (each step near zero): effect at ten years, post-only reserves only | 6 |  |  | -0.60 | 0.72 |

## Purple urchin, biomass

| Step | What | Histories or reserves | Plain average | Median | Pooled | SE |
|---|---|---|---|---|---|---|
| 0 | Paper Table 2 (multilevel model) | 10 |  |  | -0.58 | 0.27 |
| 1 | Paper's histories, pooled simply | 10 | -1.01 | -1.00 | -0.82 | 0.22 |
| 1 |   of which: mean after (siting included) | 5 | -0.39 | -0.48 | -0.43 | 0.24 |
| 1 |   of which: line over after years | 5 | -1.64 | -1.42 | -1.38 | 0.24 |
| 1 |   of which: KFM | 2 | -1.40 | -1.40 | -1.27 | 1.50 |
| 1 |   of which: PISCO | 8 | -0.91 | -1.00 | -0.84 | 0.15 |
| 2a | Paper's annual ratios: inside against outside, after years | 10 | -0.65 | -0.44 | -0.59 | 0.24 |
| 2b | Paper's annual ratios: after minus before | 1 | -0.02 | -0.02 |  |  |
| 2c | Paper's annual ratios: slope over after years x 10 | 10 | -0.95 | -1.06 | -0.88 | 0.33 |
| 4 | Joint model (trend carries forward): effect at ten years | 13 |  |  | -0.66 | 0.58 |
| 4 | Joint model (each step near zero): effect at ten years | 13 |  |  | -0.19 | 0.43 |
| 4 | Joint model (trend carries forward): effect at ten years, mean mass | 13 |  |  | +0.25 | 0.18 |
| 4 | Joint model (each step near zero): effect at ten years, mean mass | 13 |  |  | +0.09 | 0.13 |

## Red urchin, density

| Step | What | Histories or reserves | Plain average | Median | Pooled | SE |
|---|---|---|---|---|---|---|
| 0 | Paper Table 2 (multilevel model) | 16 |  |  | -0.46 | 0.46 |
| 1 | Paper's histories, pooled simply | 16 | -0.58 | -0.47 | -0.51 | 0.22 |
| 1 |   of which: mean after (siting included) | 5 | +0.60 | +0.58 | +0.58 | 0.13 |
| 1 |   of which: line over after years | 7 | -1.54 | -1.23 | -1.50 | 0.33 |
| 1 |   of which: before-and-after (step, line or curve) | 4 | -0.37 | -0.45 | -0.33 | 0.42 |
| 1 |   of which: KFM | 5 | -0.48 | -0.43 | -0.41 | 0.45 |
| 1 |   of which: LTER | 1 | -1.19 | -1.19 |  |  |
| 1 |   of which: PISCO | 10 | -0.56 | -0.67 | -0.51 | 0.29 |
| 2a | Paper's annual ratios: inside against outside, after years | 16 | -0.07 | +0.36 | -0.02 | 0.20 |
| 2b | Paper's annual ratios: after minus before | 7 | -0.27 | -0.42 | -0.24 | 0.26 |
| 2c | Paper's annual ratios: slope over after years x 10 | 16 | -0.81 | -0.65 | -0.70 | 0.17 |
| 3a | Joint model's surveys as ratios: inside against outside, after years | 18 | -0.09 | +0.17 | -0.08 | 0.20 |
| 3b | Joint model's surveys as ratios: after minus before | 8 | +0.05 | -0.04 | -0.07 | 0.37 |
| 3c | Joint model's surveys as ratios: slope over after years x 10 | 18 | -0.68 | -0.78 | -0.87 | 0.17 |
| 3a | Joint model's surveys as ratios: inside against outside, after years (offset doubled) | 18 | -0.09 | +0.17 | -0.08 | 0.20 |
| 3b | Joint model's surveys as ratios: after minus before (offset doubled) | 8 | -0.01 | -0.03 | -0.07 | 0.36 |
| 3c | Joint model's surveys as ratios: slope over after years x 10 (offset doubled) | 18 | -0.68 | -0.77 | -0.91 | 0.17 |
| 3a | Joint model's surveys as ratios: inside against outside, after years (offset halved) | 18 | -0.08 | +0.17 | -0.07 | 0.20 |
| 3b | Joint model's surveys as ratios: after minus before (offset halved) | 8 | +0.11 | -0.06 | -0.07 | 0.38 |
| 3c | Joint model's surveys as ratios: slope over after years x 10 (offset halved) | 18 | -0.69 | -0.78 | -0.84 | 0.16 |
| 4 | Joint model (trend carries forward): siting difference at designation (average over reserves, log scale per station; not a ratio of means) | 13 |  |  | +0.98 | 0.35 |
| 4 | Joint model (trend carries forward): effect at ten years | 13 |  |  | -0.63 | 0.36 |
| 5 | Joint model (trend carries forward): effect at ten years, before-and-after reserves only | 7 |  |  | -0.16 | 0.35 |
| 5 | Joint model (trend carries forward): effect at ten years, post-only reserves only | 6 |  |  | -1.13 | 0.64 |
| 4 | Joint model (each step near zero): siting difference at designation (average over reserves, log scale per station; not a ratio of means) | 13 |  |  | +0.89 | 0.35 |
| 4 | Joint model (each step near zero): effect at ten years | 13 |  |  | -0.27 | 0.29 |
| 5 | Joint model (each step near zero): effect at ten years, before-and-after reserves only | 7 |  |  | -0.03 | 0.27 |
| 5 | Joint model (each step near zero): effect at ten years, post-only reserves only | 6 |  |  | -0.32 | 0.46 |

## Red urchin, biomass

| Step | What | Histories or reserves | Plain average | Median | Pooled | SE |
|---|---|---|---|---|---|---|
| 0 | Paper Table 2 (multilevel model) | 14 |  |  | +0.57 | 0.26 |
| 1 | Paper's histories, pooled simply | 14 | +0.32 | +0.40 | +0.36 | 0.22 |
| 1 |   of which: mean after (siting included) | 9 | +0.62 | +0.25 | +0.62 | 0.20 |
| 1 |   of which: line over after years | 2 | -2.00 | -2.00 | -2.02 | 0.53 |
| 1 |   of which: before-and-after (step, line or curve) | 3 | +0.97 | +1.10 | +1.19 | 0.21 |
| 1 |   of which: KFM | 5 | +0.32 | +1.10 | +0.32 | 0.66 |
| 1 |   of which: LTER | 1 | +1.09 | +1.09 |  |  |
| 1 |   of which: PISCO | 8 | +0.22 | +0.25 | +0.28 | 0.20 |
| 2a | Paper's annual ratios: inside against outside, after years | 14 | +0.52 | +0.48 | +0.54 | 0.19 |
| 2b | Paper's annual ratios: after minus before | 5 | +0.60 | +0.46 | +0.65 | 0.28 |
| 2c | Paper's annual ratios: slope over after years x 10 | 14 | -0.38 | -0.12 | -0.39 | 0.24 |
| 4 | Joint model (trend carries forward): effect at ten years | 13 |  |  | -0.22 | 0.41 |
| 4 | Joint model (each step near zero): effect at ten years | 13 |  |  | -0.06 | 0.36 |
| 4 | Joint model (trend carries forward): effect at ten years, mean mass | 13 |  |  | +0.41 | 0.20 |
| 4 | Joint model (each step near zero): effect at ten years, mean mass | 13 |  |  | +0.21 | 0.18 |

## Giant kelp, density

| Step | What | Histories or reserves | Plain average | Median | Pooled | SE |
|---|---|---|---|---|---|---|
| 3a | Joint model's surveys as ratios: inside against outside, after years | 18 | +0.62 | +0.55 | +0.62 | 0.23 |
| 3b | Joint model's surveys as ratios: after minus before | 8 | +1.11 | +1.36 | +0.97 | 0.42 |
| 3c | Joint model's surveys as ratios: slope over after years x 10 | 18 | +1.45 | +1.26 | +1.31 | 0.38 |
| 3a | Joint model's surveys as ratios: inside against outside, after years (offset doubled) | 18 | +0.58 | +0.51 | +0.58 | 0.22 |
| 3b | Joint model's surveys as ratios: after minus before (offset doubled) | 8 | +1.01 | +1.29 | +0.89 | 0.39 |
| 3c | Joint model's surveys as ratios: slope over after years x 10 (offset doubled) | 18 | +1.44 | +1.25 | +1.24 | 0.36 |
| 3a | Joint model's surveys as ratios: inside against outside, after years (offset halved) | 18 | +0.66 | +0.62 | +0.65 | 0.23 |
| 3b | Joint model's surveys as ratios: after minus before (offset halved) | 8 | +1.22 | +1.37 | +1.04 | 0.44 |
| 3c | Joint model's surveys as ratios: slope over after years x 10 (offset halved) | 18 | +1.45 | +1.26 | +1.38 | 0.39 |
| 4 | Joint model (trend carries forward): siting difference at designation (average over reserves, log scale per station; not a ratio of means) | 13 |  |  | -0.34 | 0.54 |
| 4 | Joint model (trend carries forward): effect at ten years | 13 |  |  | +0.96 | 0.58 |
| 5 | Joint model (trend carries forward): effect at ten years, before-and-after reserves only | 7 |  |  | +0.82 | 0.83 |
| 5 | Joint model (trend carries forward): effect at ten years, post-only reserves only | 6 |  |  | +1.10 | 0.62 |
| 4 | Joint model (each step near zero): siting difference at designation (average over reserves, log scale per station; not a ratio of means) | 13 |  |  | -0.30 | 0.55 |
| 4 | Joint model (each step near zero): effect at ten years | 13 |  |  | +0.60 | 0.54 |
| 5 | Joint model (each step near zero): effect at ten years, before-and-after reserves only | 7 |  |  | +0.32 | 0.63 |
| 5 | Joint model (each step near zero): effect at ten years, post-only reserves only | 6 |  |  | +0.56 | 0.58 |

## Giant kelp, biomass

| Step | What | Histories or reserves | Plain average | Median | Pooled | SE |
|---|---|---|---|---|---|---|
| 0 | Paper Table 2 (multilevel model) | 31 |  |  | +1.14 | 0.26 |
| 1 | Paper's histories, pooled simply | 31 | +1.22 | +1.12 | +1.08 | 0.19 |
| 1 |   of which: mean after (siting included) | 6 | +0.50 | +0.27 | +0.46 | 0.46 |
| 1 |   of which: line over after years | 7 | +1.74 | +1.67 | +1.69 | 0.23 |
| 1 |   of which: before-and-after (step, line or curve) | 18 | +1.26 | +0.91 | +1.04 | 0.25 |
| 1 |   of which: KFM | 5 | +2.34 | +2.44 | +2.14 | 0.43 |
| 1 |   of which: LTER | 2 | +2.28 | +2.28 | +2.08 | 0.74 |
| 1 |   of which: Landsat | 14 | +0.78 | +0.76 | +0.75 | 0.24 |
| 1 |   of which: PISCO | 10 | +1.06 | +1.30 | +1.00 | 0.36 |
| 2a | Paper's annual ratios: inside against outside, after years | 17 | +0.82 | +0.57 | +0.74 | 0.23 |
| 2b | Paper's annual ratios: after minus before | 8 | +1.29 | +1.48 | +1.13 | 0.40 |
| 2c | Paper's annual ratios: slope over after years x 10 | 17 | +1.66 | +1.64 | +1.55 | 0.25 |
| 4 | Joint model (trend carries forward): effect at ten years | 13 |  |  | +0.96 | 0.59 |
| 4 | Joint model (each step near zero): effect at ten years | 13 |  |  | +0.57 | 0.55 |
| 4 | Joint model (trend carries forward): effect at ten years, mean mass | 13 |  |  | +0.01 | 0.13 |
| 4 | Joint model (each step near zero): effect at ten years, mean mass | 13 |  |  | -0.03 | 0.11 |

