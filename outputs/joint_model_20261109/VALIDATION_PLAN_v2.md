# Validation plan, round 2: randomization placebo, evidence split, pre-trend

Written 2026-11-11 (late evening) before any of the runs below were launched on GRIT. It adds to
`MASS_VALIDATION_PLAN.md` (rules M1–M11, unchanged). Nothing here is adopted for the paper.

## Why a second round

Round 1 exposed three weaknesses in how the model was being tested, not only in the model:

1. **The space placebo was not independent across seeds.** Each seed promoted one reference station per MPA from a short
   list, so two seeds shared most of their relabelled stations (seeds 2 and 10 shared 7 of 9 for lobster). Ten seeds
   were closer to two or three tests than to ten. The lobster lean (+0.46, all positive) could not be read either way.
2. **Pooled effects hid how thin the before data are.** At the six lobster MPAs with before-and-after surveys, protected
   stations counted 48 lobsters before designation, 40 of them at Scorpion. Kelp has the same shape at the three KFM
   island MPAs (5 plants counted before at protected stations across Harris Point, Santa Barbara Island and Scorpion).
3. **The sheephead time placebo fails for a real reason.** Sheephead at the three to-be-protected KFM stations fell
   relative to the four reference stations between 1985 and 2002 (0.09 vs 0.48 fish per 100 m² in 2000–2002). The model
   treats the before period as flat, so any fake designation year inside that decline returns a negative effect. Two
   structure options (MPA-cluster random walk, station slopes) did not change this and were not adopted.

Exploratory runs already seen before these rules were written (short, 2 chains, sandbox): the cluster and station-slope
options on the sheephead time placebo and the lobster space placebo; lobster density split by before-and-after vs
post-only MPAs; sheephead with `--pretrend` on the real data and on the 1993 and 1996 time placebos. Their results are in
`NOTES_20261111_evening.md` and are not evidence for or against the rules below.

## R1. Randomization (label-permutation) placebo — replaces the space placebo as the test of record

`--permute-seed k`: within each MPA and its home reference stations, the protected label is dealt at random to the same
number of stations. Every station and every survey is kept; designation years are unchanged. Each station has the same
chance of a protected label, so if the model is unbiased the permuted effects centre on zero whether or not protection
has a real effect (a real effect adds spread, not a shift).

- 40 permutations per taxon (seeds 1–40), density only, recommended structure, both curve priors, 4 chains × (1500 +
  3000). A permutation counts only if R-hat ≤ 1.05 and divergences ≤ 1% of draws; report how many are dropped.
- **R1a (bias).** The mean of the permuted E10 must lie within 2 standard errors (SD/√n) of zero. Fail = the model
  manufactures an effect from the design alone; report the size of the shift.
- **R1b (calibration).** The 90% interval of each permuted fit should cover zero in 80–97% of permutations (the range
  that 40 draws from a true 90% allow, roughly). Below 80% = intervals too narrow for this design.
- **R1c (the real effect against the design).** Report the randomization p-value: (1 + number of permutations with
  |E10| ≥ |real E10|) / (n + 1). This is reported for every taxon and is not a gate; it is the number that says how
  unusual the real estimate is among estimates the same data could give with the labels shuffled.
- The old space placebo runs (v_big tasks 0–49) are still reported, labelled as non-independent.

## R2. Evidence split — reported for every taxon, no gate

`evidence_table.py` writes, per taxon and MPA, the animals counted and survey years before and after designation at
protected and reference stations. For each taxon and curve prior, fit three versions (density + mass where mass exists):
all MPAs (the main fit), before-and-after MPAs only, post-only MPAs only. Report all three side by side with the counts.

Labels used in any summary (proposed; Adrian to confirm the cut-offs):
- "before-and-after supported": the before-and-after-only interval excludes zero **and** at least three of those MPAs
  each counted ≥ 20 animals at protected stations before designation.
- "rests on protected-vs-reference contrast after designation": anything else. For these, the estimate cannot separate
  protection from where reserves were sited.

By the counts alone, lobster cannot earn the first label (one MPA with ≥ 20 before). That is known before the runs.

## R3. Sheephead pre-trend — a sensitivity, reported next to the flat baseline

`--pretrend`: at MPAs with before data, protected stations get their own linear trend relative to references (per MPA,
partially pooled), extended through the after years; the protection function then measures departure from that line.
This is the model-based version of the Nov 3 policy ("flat-baseline recovery main and pre-trend sensitivity").

- Run for all five taxa, both curve priors, real data; and for the sheephead and urchin time placebos (the same fake
  years as round 1: 1990, 1993, 1996, 1999).
- **R3a.** With `--pretrend`, the sheephead time-placebo 90% interval should cover zero at ≥ 3 of the 4 fake years.
  If it does not, the decline is not linear enough for this fix and the sheephead before-and-after contrast stays flagged.
- **R3b.** Report the real-data E10 flat and pre-trend side by side for every taxon. Neither is chosen by its result.
  The flat baseline stays the main specification unless Adrian decides otherwise; where the two disagree in sign or the
  pre-trend interval is more than twice as wide, the summary says the effect depends on the baseline assumption.
- Simulation check before trusting it: 20 simulated datasets with a true pre-trend (mu 0.5 per decade, sd 0.3) and a
  known E10; coverage of E10 with `--pretrend` must be within 0.80–0.99 at nominal 90%, and without it is reported to
  show the size of the bias the option removes.

## What stays open and is not addressed by this round

- Kelp space placebo from stage 1 (R1 will retest it with independent labels).
- Urchin count misfit in the posterior checks; the dependence of urchin and kelp effects on the curve prior.
- The lobster → urchin fluctuation link is not identified (sign flips between chains).
- Landsat canopy and the controlled comparison against the paper's numbers are separate pieces of work.
