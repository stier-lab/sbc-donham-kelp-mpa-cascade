# Validation plan for the joint density, mean mass and biomass model

Written 2026-11-11, 18:45, before any of the runs below had results. Review only. A check that fails stays failed in
the report; no rule is loosened afterwards, and no specification, prior, site or cohort is chosen by the direction or
significance of a result. Effects are log changes at ten years after designation (density, mean mass, biomass);
intervals are 90%. "2 SE" is two binomial standard errors at the stated number of intervals.

| # | Check | Runs | Pass rule |
|---|---|---|---|
| M1 | Convergence | main fits: 5 taxa, random-walk prior (2 seeds) and independent-increment prior | R-hat < 1.01, effective sample size > 400, zero divergent transitions, every fit |
| M2 | Reproducibility | the two seeds | density, mass and biomass effects differ by < 0.05 |
| M3 | Calibration when the model is right | simulated data on the real designs, 20 data sets per taxon in each of: effect present (`std`), no effect on density or mass (`null`), no effect on mass plus a regional trend in size (`trend`) | coverage of the true mass effect and of the true biomass effect, pooled over taxa (100 intervals each): at least 0.90 − 2 SE = 0.84. In `null` and `trend`: share of mass-effect intervals that exclude zero at most 0.10 + 2 SE = 0.16. Mean error of the mass effect per taxon below 0.10 in absolute value |
| M4 | Robustness when an assumption is wrong | 20 data sets per taxon in each of: two size cohorts instead of a Gamma (`bimodal`); extra visit-to-visit noise in mean size (`visitnoise`); mean size falling with station density while the model assumes independence (`coupled`) | same coverage rule as M3. A failure here is reported as a limit of the model and named in the summary |
| M5 | Placebo in time | KFM surveys before the real designation, fake designation in 1992; purple urchin, red urchin, kelp; both priors (sheephead and lobster have no KFM sizes before 2003) | mass-effect and biomass-effect intervals include zero |
| M6 | Placebo in space | protected stations dropped, one reference station per MPA relabelled as protected; 10 random relabellings per taxon | per taxon and effect: at most 25% of intervals exclude zero, and not (same sign in all ten with mean beyond ±0.22) |
| M7 | Predictive check of the mass part | main random-walk fits, 200 posterior draws per chain | no protected-versus-reference contrast statistic outside its 90% replicate interval; at most 10% of all statistics outside |
| M8 | Prediction of data not shown to the model | 20% of interior station-years of size data held out, 2 random splits per taxon | per taxon: 90% predictive intervals cover 85% to 95% of held-out group means |
| M9 | Leave one MPA out | each MPA dropped in turn, every taxon | no single MPA reverses the sign of an effect whose main interval excludes zero; shifts above 0.25 are listed |
| M10 | Independent arithmetic | `code/joint_model/test_multi.py` on the real inputs | density part equals the single-taxon code; mass part equals a separate NumPy/SciPy calculation (relative difference < 1e-9) |
| M11 | Input integrity | `code/joint_model/prep_mass.py` | every mass group reproduces the integrated review's value from source rows |

Already known and not covered by this plan (density part, 2026-11-10): the sheephead time placebo and the kelp space
placebo fail; the count layer misfits urchin counts; calibration misses its rule narrowly; the two curve priors give
different effect sizes. The biomass effect inherits these through the density effect.
