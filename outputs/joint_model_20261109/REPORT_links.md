# Five species in one fit, with links between them

Fits `v_full_rw` and `v_full_pw` (density + mean mass for all five species, two-coefficient links, 4 chains × (3000 + 6000),
28.6 and 29.6 hours on GRIT, finished 2026-11-13). Review only; nothing adopted.

Each link has two coefficients:
- **year-to-year (beta):** when the source species is above its own usual level at a station one year, how far the target
  species is above or below its level the next year;
- **through protection (beta_prot):** how much of the target's protection effect is carried by the source's protection effect.

Only the second is a protection-driven cascade. The first is a correlation in time, open to shared local drivers.

## Convergence

| Fit | Worst R-hat | Smallest effective sample | Divergences |
|---|---|---|---|
| Trend carries forward (rw) | 1.014 (red urchin regional drift) | 402 | 0 |
| Each step near zero (pw) | 1.017 (red urchin regional drift) | 205 (below the 400 used for main fits) | 0 |

## Links (90% intervals)

| Link | Year-to-year, rw | Year-to-year, pw | Through protection, rw | Through protection, pw |
|---|---|---|---|---|
| Lobster → purple urchin | +1.52 [+1.17, +1.85] | +1.53 [+1.17, +1.87] | +0.23 [−0.08, +0.53] | +0.15 [−0.14, +0.44] |
| Sheephead → purple urchin | −0.42 [−0.71, −0.17] | −0.43 [−0.70, −0.15] | −0.05 [−0.75, +0.63] | −0.11 [−0.80, +0.59] |
| Lobster → red urchin | +1.04 [+0.80, +1.27] | +1.05 [+0.81, +1.28] | +0.06 [−0.14, +0.25] | −0.00 [−0.18, +0.18] |
| Sheephead → red urchin | −0.30 [−0.48, −0.11] | −0.30 [−0.49, −0.11] | +0.13 [−0.49, +0.69] | +0.08 [−0.52, +0.66] |
| Purple urchin → kelp | −0.45 [−0.61, −0.29] | −0.44 [−0.61, −0.28] | −0.44 [−0.95, +0.12] | −0.47 [−1.00, +0.03] |
| Red urchin → kelp | −0.22 [−0.51, +0.06] | −0.22 [−0.52, +0.06] | −0.06 [−0.83, +0.62] | −0.09 [−0.77, +0.64] |

The year-to-year coefficients are now well determined (in the short fit of 11 November the lobster links had R-hat 7).

## Effects at ten years, with the part carried by each link

| Species | Direct (rw) | Total (rw) | Carried through links (rw) | Total (pw) | Single-species fit (rw) |
|---|---|---|---|---|---|
| Spiny lobster | +1.79 | +1.79 [+1.15, +2.47] | — | +1.87 | +1.75 |
| Sheephead | +0.20 | +0.20 [−0.09, +0.48] | — | +0.14 | +0.20 |
| Purple urchin | −1.30 | −0.90 [−1.81, +0.00] | via lobster +0.41 [−0.09, +0.98]; via sheephead −0.01 | −0.25 | −0.91 |
| Red urchin | −0.70 | −0.58 [−1.13, −0.03] | via lobster +0.10; via sheephead +0.02 | −0.23 | −0.63 |
| Giant kelp (divers) | +0.55 | +0.91 [−0.00, +1.90] | via purple urchin +0.33 [−0.15, +0.94]; via red urchin +0.03 | +0.44 | +0.96 |

Biomass totals (rw): lobster +2.18, sheephead +0.52, purple urchin −0.64, red urchin −0.16, kelp +0.91. Mass effects match
the single-species fits within 0.01.

## Reading

1. **The joint fit reproduces the single-species fits.** Every total effect is within 0.06 of the species fitted alone.
   Linking the species does not change any headline number.
2. **The year-to-year links have the cascade sign for sheephead and for purple urchin → kelp, and the wrong sign for
   lobster.** A lobster-rich year at a station is followed by more urchins, not fewer. The 11 November simulations showed
   this is what a shared local driver (good habitat or recruitment for both) produces when no real predation link exists.
   These coefficients cannot be read as predation.
3. **No link carries a protection effect that is distinguishable from zero.** The closest is purple urchin → kelp
   (−0.44 [−0.95, +0.12]; pw −0.47 [−1.00, +0.03]). Under the main curve prior about a third of kelp's ten-year effect
   (+0.33 of +0.91) is attributed to the purple urchin decline, with wide uncertainty; under the other prior, +0.08 of +0.44.
4. **Lobster's "via" effect on urchins is positive,** opposite to a cascade, and spans zero. Protection's effect on urchins
   is carried by the direct term, not by lobsters.

What this supports saying: the data are consistent with an urchin-to-kelp step under protection but do not establish it,
and they give no evidence that lobster recovery is what reduced urchins inside reserves.
