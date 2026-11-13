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
