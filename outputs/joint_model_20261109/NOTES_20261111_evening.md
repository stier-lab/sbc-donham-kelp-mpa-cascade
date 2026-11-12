# Evening notes, 2026-11-11 (to commit when the Mac is reachable)

Sandbox tests of two structure options on placebos only (short fits, 2 chains):
- Sheephead time placebo (KFM before 2003, fake designation 1993), E10: current -2.43 [-4.39, -0.43];
  local trend per MPA cluster -2.23 [-3.97, -0.48]; station slopes -2.34 [-4.33, -0.28]; both -2.17 [-4.13, -0.39].
- Lobster space placebo seed 2, E10: current +1.07 [-0.02, +2.19]; local trend +1.00 [-0.10, +2.13];
  station slopes +1.09 [-0.05, +2.13]. Seed 10 current +0.72 [-0.33, +1.70].
Neither option changes a placebo. Not adopted; flags --cluster-drift and --station-trend stay off by default.

Raw aggregates behind the placebos:
- KFM sheephead, protected-to-be vs reference, fish per 100 m2: 1985-89 1.07 vs 1.22; 1990-92 0.57 vs 0.97;
  1993-96 0.95 vs 1.40; 1997-99 0.43 vs 1.07; 2000-02 0.09 vs 0.48 (6 fish in 6,800 m2 vs 46 in 9,600 m2).
  Log ratio: -0.13, -0.53, -0.39, -0.92, -1.69. Three protected and four reference stations.
- Lobster space placebo: seeds 2 and 10 relabel 7 of 9 identical stations. Relabelled stations with zero lobsters
  counted before 2012: Whites Point (0 -> 10.8 per 1000 m2; references 3.3 -> 17.1), Nicholas Canyon W (0 -> 16.7;
  references 1.5 -> 19.4), Escondido W (0 -> 0.9).
- Lobster before-designation counts at protected stations of the six before-after MPAs: 0, 0, 1, 6, 1, 40 (Scorpion).

Lobster density effect by evidence split (sandbox, 1500 + 3000, 2 chains): before-after MPAs only +1.72 [+0.93, +2.47];
post-only MPAs +1.11 [+0.10, +2.15] (R-hat 1.024, not converged); without Scorpion +1.73 [+0.95, +2.50].

## Late evening (sandbox only; Mac and GRIT unreachable)

**Evidence table** (`evidence_table.py`, totals only). Animals counted at protected stations before designation, at the
MPAs that have before-and-after surveys:

| Taxon | BA MPAs | Counted before (protected) | Largest MPA's share | BA MPAs with < 20 counted |
|---|---|---|---|---|
| Lobster | 6 of 13 | 48 | 83% (Scorpion) | 5 |
| Sheephead | 7 of 17 | 884 | 33% | 1 |
| Purple urchin | 7 of 13 | 118,944 | 32% | 0 |
| Red urchin | 7 of 13 | 21,819 | 35% | 0 |
| Kelp | 7 of 13 | 28,594 | 64% (Naples) | 3 (Harris Point 1, Santa Barbara Island 0, Scorpion 4 plants) |

Kelp shares lobster's weakness at the three KFM island MPAs: almost no plants at protected stations before 2003, so the
"before" level there is a near-zero count and the log change is set by the prior and the zero handling as much as by data.

**Pre-trend option on the sheephead time placebo** (short runs, 2 chains, rw prior):

| Fake year | Flat baseline | With `--pretrend` | Estimated pre-trend (per decade) |
|---|---|---|---|
| 1993 | −2.43 [−4.39, −0.43] | −1.85 [−4.15, +0.63] | −0.48 [−1.63, +0.66] |
| 1996 | (round 1, GRIT: −3.18) | −2.53 [−5.19, +0.50] | −0.48 [−1.38, +0.51] |

The option widens the interval to include zero and moves the estimate about a quarter of the way. It does not remove the
problem: the decline at the three to-be-protected KFM stations is not a straight line (it steepens after 1996), and three
MPAs cannot pin a trend down. This is a description of the before period, not a fix. Rules for the full runs are in
`VALIDATION_PLAN_v2.md` (written after these two short runs; said so there).

**New options, all off by default:** `--permute-seed`, `--evidence {ba,post}`, `--pretrend`, `--sim-pretrend`. Code checks
1–3 still pass with defaults. Task lists: `jobs/v2_perm.tsv` (400), `v2_split.tsv` (20), `v2_pre.tsv` (34), `v2_sim.tsv` (80).
