# Joint hierarchical model: architecture, priors and assumptions (v2)

Status: design v2, 2026-11-09. No model has been fitted. Nothing here changes the frozen paper pipeline or the integrated review.
Evidence: `outputs/joint_model_design_20261109/` (data profile, `profile_joint_model_data.py`) and
`design_checks.json` (`code/tools/check_joint_model_design.py`). All checks use density histories from the integrated review.

## 0. What changed from v1, and why

| # | v1 | v2 | Evidence or reason |
|---|---|---|---|
| 1 | Step + linear trend on the log scale | A regional **response curve by years since designation** (piecewise-linear, 4-year knots) plus an MPA-specific rate deviation | Mean change in log ratio from ages 0–3: lobster +0.41, +0.83, +1.54, +1.18 (rises, then flattens after ~12 yr); purple urchin −0.11, −0.63, −0.94, −1.43 (keeps falling); red urchin +0.08, −0.54, −0.62, −0.83. One shape does not fit all taxa. v1's slope prior also allowed a ×358 change over 20 years. |
| 2 | MPA effects independent given the model | A **protected × year shock shared across MPAs** | Detrended post-protection residuals correlate between MPAs: mean r = 0.05 (red urchin), 0.07 (lobster), 0.10 (sheephead), −0.02 (purple urchin). With 13 MPAs, r = 0.05–0.10 inflates the variance of a pooled annual mean 1.6–2.2×. Ignoring it makes pooled intervals too narrow. |
| 3 | AR(1) with one ρ prior, N(0.3, 0.3), read off the data | **AR(1) process plus white noise, taxon-specific ρ**, program-specific noise | Lag-1 autocorrelation: purple urchin 0.49, red urchin 0.25, sheephead 0.07, lobster −0.12. By program: KFM 0.32, LTER 0.31, PISCO 0.02. PISCO histories are dominated by sampling noise (median 4 transects per visit). A single ρ is wrong, and the white-noise share differs by program. |
| 4 | Lognormal lengths, allometry priors "from published SEs" | **Mean mass modelled directly** (Gamma likelihood on individual mass, mean-parameterized); allometry fixed, exponent varied in a sensitivity | Biomass needs the arithmetic mean mass. A mean-parameterized likelihood estimates it consistently even when sizes are bimodal (KFM). The repository stores allometry coefficients without standard errors, so v1's prior had no source. Exponent uncertainty barely matters for contrasts (§4.2). |
| 5 | "Density" treated as one quantity | **Density is defined per program**, with an optional common ≥25 mm definition | PISCO counts and measures urchins ≥25 mm only; KFM and LTER count all sizes. Recruitment pulses enter KFM density and not PISCO density. |
| 6 | Visit effect "only identified with multiple visits" | Visit term is the white-noise component for every visit | It is identified from the autocovariance of long series (KFM has 18–21 before-years), not only from repeat visits. v1 was wrong on this. |
| 7 | Units aggregated to visits | Aggregated to **visit × zone** | Negative-binomial sums are exact only when units share a mean. Zones differ. |
| 8 | Siting difference λ centred on zero | **Mean siting difference estimated** from MPAs with before data; post-only MPAs get a free intercept | A non-zero mean is Bob Miller's siting concern. The model can estimate it directly instead of assuming it away. |
| 9 | Cascade test = correlation of MPA effects | That, plus a **dynamic link**: last year's predator density predicts this year's urchin change at the same station | 13 MPAs cannot estimate a correlation well. The station-year link uses every station and year. |
| 10 | Validation by simulation only | Adds **placebo tests** and a **frozen specification** before fitting real protection labels | Three island MPAs have 18–21 KFM years before 2003 for four taxa. A fake designation in the early 1990s should give no effect. |
| 11 | Priors tuned to the same data; prior check passed after I widened two scales | Generic weakly informative scales, stated per decade; prior check targets implausible outcomes | Tuning priors to the data and then showing they cover the data proves little. |
| 12 | No statement of what can be detected | **Detectable effect sizes stated** | Current pooled half-widths: purple urchin ±0.08/yr, red urchin ±0.10, lobster ±0.15, sheephead ±0.08. That is a ×2.2 to ×4.4 change per decade. No model removes this limit. |

Checked and dropped: a "take allowed" covariate. Where recorded, no retained MPA allows lobster or urchin take (`mpa_metadata.csv`). The record is blank for 14 of 22 MPAs; fill those before relying on it.

Not verified (the connection to Adrian's computer dropped): whether PISCO swath counts were subsampled in the field and scaled up. If they were, the negative binomial must use the area actually searched (§6, item 13).

## 1. Estimands

The model reports one quantity per taxon, MPA and age: **f(a) = log of (density at the protected stations a years after designation) ÷ (density at the same stations had they not been protected)**. The unprotected counterfactual comes from the stations' own before-years where they exist, and from reference stations through shared region-year effects.

Two regimes, never mixed silently:

- **Before–after MPAs** (6–8 of 13 for urchin and lobster density; 25 histories for kelp): f(a) is identified from zero at designation. Causal reading needs parallel trends between protected and reference stations (testable in before-years; §8).
- **Post-only MPAs:** only *change after the first survey* is identified. The level at first survey mixes siting and early response. These MPAs inform the shape and rate of the response, not its level.

Headline numbers: **E10 = g(10)**, the regional mean log change at 10 years; the full curve at 4, 8, 12, 16 years; and the mean annual rate over the first decade (comparable to current slopes). The same three for mean individual mass and for biomass.

Predicted directions, stated before fitting: sheephead +, lobster +, purple urchin −, kelp +. Red urchin has no prediction: it is fished (protection should raise it) and prey (predators should lower it).

**What can be detected.** With 13 MPAs and this noise, pooled effects smaller than about a doubling or halving per decade are unlikely to be resolved for urchins and sheephead, and smaller than about ×4 for lobster. The joint model should use the data more efficiently, but by how much is unknown until stage 1 runs.

## 2. Data facts that drive the design

| Fact | Consequence |
|---|---|
| Integer counts per unit with known area (PISCO 60 m² swaths; KFM 1–2 m² quadrats, 40–60 m² bands, 100–200 m² fish transects; LTER 1 m² quadrats, 300 m² lobster bands) | Count likelihood with area offset; no +offset log ratios |
| Zeros: lobster 78–88% of units, LTER red urchin 75%, sheephead 45–63% | Negative binomial; test for zero-inflation, don't assume it |
| Within-visit overdispersion 1/φ = 0.02–2.8 | Dispersion per program × taxon |
| **Median 1 protected station per history (max 3); references median 2 (max 9)** | "MPA effect" is mostly one station's effect. No within-MPA replication; pooling across MPAs carries the inference |
| Protected and reference annual means correlate 0.23–0.76 | Shared region × year effect |
| Before-years: 6–8 of 19 density histories have ≥3; **urchin biomass has none** | Urchin biomass supports post-only inference |
| PISCO size pools small (urchins median 7–10 animals; lobster 3; 12–34% singletons). KFM pools ~200 | Partial pooling of mean mass; open decision on how far to pool (§9) |
| KFM within-pool SD of log mass 1.1–1.4 (mixed cohorts) | No lognormal size assumption |
| 19–30% of PISCO/LTER station-years have >1 visit; protected and reference stations are sometimes surveyed weeks apart | Month effect and visit-level noise |
| Designation is 2003 (islands) or 2012 (mainland) for every MPA | Ages beyond 11 years come only from the 2003 island cohort |

## 3. Indices

*k* taxon; *p* program; *m* MPA; *c* cluster (an MPA with its reference stations; MPAs sharing references form one cluster); *r* region; *s* station (once); *t* year; *v* visit; *z* zone; *a* = *t* − *T_m* years since designation; *A* = 1 if *t* ≥ *T_m*; *P_s* = 1 for protected stations.

## 4. Architecture

### 4.1 Counts

    y[v,z] ~ NegBin(mean = area[v,z] · exp(η[s,t,k] + ν[v,k] + ζ[s,z,k] + m_k(month_v)),  n[v,z]·φ[p,k])

- Units are summed within visit × zone (exact for equal-area units).
- ν[v,k] ~ N(0, σ_ν[p,k]): **white noise per visit**. Covers within-year timing, conditions and patchiness beyond the negative binomial. Program-specific, because PISCO's is much larger.
- ζ[s,z,k] ~ N(ζ̄[z,k], σ_ζ): **station × zone offsets**. Main run keeps the approved shared-zone rule. Because zones are modelled, an all-zones run is a clean sensitivity that recovers years the rule drops.
- m_k(month): month effect per taxon (sum to zero). Removes bias when protected and reference stations are surveyed in different months.

### 4.2 Mean mass and biomass (fished taxa)

    mass_j = α_k · length_j^β_k                       (fixed published coefficients)
    mass_j ~ Gamma(shape κ[p,k], mean M[s,t,k])        (frequency-weighted)
    log M[s,t,k] = c[s,k] + γᴹ[r,t,k] + P_s·A·(dᴹ[m,k] + Eᴹ[m,k]·a/10) + e[s,t,k]
    biomass density[s,t,k] = exp(η[s,t,k]) · M[s,t,k]

- **Why Gamma on mass:** M is the arithmetic mean, the quantity biomass needs. The estimate of M stays consistent when sizes are skewed or bimodal; only its precision depends on the variance assumption, which posterior predictive checks test.
- **Additivity:** log biomass effect = log density effect + log mean-mass effect. All three are reported. Mean mass tests "bigger animals, not just more".
- **Allometry.** The intercept α cancels in every protected/reference contrast. The exponent (≈2.9 for urchins and lobster) matters only through size differences: a 0.1 error in β moves a log contrast by 0.01 when protected animals are 10% longer, and by 0.03 when 30% longer. Year-to-year noise is ~1. So allometry is fixed in the main fit and re-run at β ± 0.1. Standard errors from the LTER coefficient package should replace ± 0.1 if they exist.
- **LTER class means:** a class mean of n animals enters as Gamma(n·κ, mean M), exact under the Gamma model.
- **KFM sheephead:** roving-diver sizes, transect counts. Used as is in the main run; a sensitivity drops KFM sheephead biomass.
- **Station-years with counts and no sizes:** M is predicted from the station effect c, region-year γᴹ and e's prior. The reach of that prediction is decision 1 (§9). The model reports, per station-year, how much of M's posterior comes from animals measured there.
- **Common ≥25 mm density (sensitivity):** for KFM and LTER, the share of measured urchins ≥25 mm, q[s,t], is modelled (binomial, partially pooled) and density ≥25 mm = exp(η)·q. Biomass is already nearly comparable across programs because small urchins weigh little.

### 4.3 Latent state

    η[s,t,k] = b[s,k] + γ[r,t,k] + χ[c,t,k] + δ[s,t,k] + P_s · F[m,k](t)

- b[s,k]: station baseline (program-specific mean and spread).
- γ[r,t,k]: region × year effect. χ[c,t,k]: cluster × year effect (an MPA and its references move together beyond the region). Both are shared by protected and reference stations, so neither can create a protection effect.
- δ[s,t,k]: station AR(1), with ρ_k per taxon and ρ^gap across missing years. White noise sits in ν (§4.1), so ρ describes the population process, not sampling.

### 4.4 Protection function

    F[m,k](t) = λ[m,k] + A · ( g_k(a) + D[m,k]·a/10 + ω[k,t] )   [+ π[m,k]·(t − T_m): pre-trend sensitivity]

- **g_k(a): regional response curve.** Piecewise-linear with knots at 0, 4, 8, 12, 16, 20 years, g_k(0) = 0. Increments between knots follow a random walk, which shrinks the curve toward a straight line unless the data bend it. A step-like, linear, saturating or delayed response are all special cases; no per-history model selection (pBACIPS's AICc step) is needed.
- **D[m,k]: MPA-specific rate deviation** (log change per decade), D ~ N(0, τ_D[k]). One random effect per MPA: with about one protected station, a second (a level shift) is not separately estimable. A level-shift version is a sensitivity.
- **ω[k,t]: protected × year shock shared across MPAs**, ω ~ N(0, σ_ω[k]). It absorbs years when protected sites as a group diverge from references (lobster 2015 +0.67, 2018 +0.71; sheephead 2016 −0.58 in the residuals). Without it, those years count as independent evidence in every MPA. It also carries the heatwave-resistance signal the paper cares about: the 2014–2016 values of ω for kelp are a direct estimate.
- **λ[m,k]: siting difference.** For MPAs with before data, λ ~ N(μ_λ[k], σ_λ[k]); **μ_λ is reported** as the average protected-minus-reference difference before protection. For post-only MPAs, λ is a free intercept (N(0, 3)), so their level tells the model nothing about the response; only their change with age does. Sensitivity: draw post-only λ from the same distribution.
- Ages beyond 11 years rest on the nine 2003 island MPAs. The curve is reported with that label.

### 4.5 Cascade links (two tests)

1. **Between-MPA:** D[m,·] across taxa ~ MVN with correlation Ω. Asks whether MPAs with faster predator gains show faster urchin declines. Weak with 13 MPAs; report the whole posterior.
2. **Within-station dynamics (stage 2b):** the urchin process gains a term β_pred·(η_pred[s,t−1] − station mean), and kelp gains β_urch·(η_urch[s,t−1] − mean). It uses every station-year, protected or not, and latent states handle predator measurement error. Region-year and cluster-year effects absorb shared environmental drivers, but unmeasured local drivers can still fake a link. Report it as an association consistent with predation, with the mediated share of the protection effect (through predators vs direct) as a derived quantity.

### 4.6 Kelp

Diver kelp uses §4.1 on frond or stipe counts, with the adopted SBC 127.5 monthly conversion and its uncertainty (largely cancelling when both roles are surveyed in the same month; the month effect covers the rest). Landsat canopy keeps a separate observation model (hurdle-lognormal on polygon values, existing pixel-overlap kernel) with its own g and D, correlated with diver kelp. Stage 4.

## 5. Priors

Log scale; HN(s) half-normal. Scales are generic and weakly informative. They are no longer fitted to the profile, which is used only to check that they are not absurd.

| Parameter | Prior | Reasoning |
|---|---|---|
| g increments (per 4 years) | N(0, σ_g), σ_g ~ HN(0.4) | Implies a 95% prior range of ×3.5 at 10 years and ×5.8 at 20 years, either direction. Covers reported reserve effects of order ×2–5 without the ×358 tail of a linear slope prior. |
| τ_D, between-MPA SD of decade rate | HN(1.0) | Same scale as the regional effect: MPAs may differ as much as the mean. Most influential prior; run at 0.5 and 2.0. |
| μ_λ; σ_λ | N(0, 1); HN(1) | Siting differences of several-fold are possible. |
| σ_ω, shared protected × year shock | HN(0.5) | Residual year means reach ±0.7. |
| σ_γ, σ_χ (region-year, cluster-year) | HN(1), HN(0.5) | Order-of-magnitude swings (heatwave, wasting disease) are real. |
| ρ_k | Beta(2, 2) on (0, 0.95) | Populations of long-lived animals do not have negative carry-over; observed negatives are noise. Centre 0.47, wide. |
| σ_δ (AR innovation) | HN(1) | |
| σ_ν (visit white noise) | HN(1), by program × taxon | PISCO noise is large. |
| 1/√φ | HN(1) | |
| σ_b (station baselines) | HN(1.5) | Stations differ up to ~10×. |
| ζ̄, σ_ζ (zones); month effects | N(0, 1), HN(0.5); N(0, 0.5) | |
| Eᴹ, dᴹ (mean-mass decade rate, step); τs | N(0, 0.5), N(0, 0.3); HN(0.3) | A decade's change in mean mass beyond ×2.7 is implausible for slow-growing animals. |
| Gamma shape κ | lognormal(0, 1) | CV of individual mass 0.4–2.7. |
| σ_c, σ_e (mass: station, station-year) | HN(1), HN(0.5) | |
| π pre-trend (sensitivity) | N(0, σ_π), σ_π ~ HN(0.05)/yr | Unregularized pre-trends in the review explode (between-MPA SD up to 11). |
| β_pred, β_urch (dynamic links) | N(0, 0.3) | A unit log change in predators shifting next year's prey by more than ~×1.8 is implausible. |
| Ω | LKJ(2); run LKJ(1), LKJ(4) | |

**Prior checks to run before fitting:** (1) prior probability that any taxon's regional 20-year effect exceeds ×20: target <2%; (2) simulated station counts stay within physically possible densities; (3) prior-only fits return the prior (code check).

## 6. Assumptions

Testable in the model:
1. Negative binomial counts; no extra zero-inflation.
2. Gamma variance function for individual mass.
3. AR(1) + white noise per station; common ρ per taxon.
4. One protection response per MPA across programs (sensitivity: program offset on D).
5. Flat baseline before protection (sensitivity: pre-trend; placebo test).
6. The response curve is the same function of age for the 2003 and 2012 cohorts (testable only for ages 0–11).

Not testable here; stated and bounded:
7. **Parallel trends:** absent protection, protected stations would have followed their references' region- and cluster-year effects. Placebo tests (§8) probe it; they cannot prove it.
8. **Equal detectability** at protected and reference stations within program and year.
9. **Measured animals represent counted animals.** Weakest for KFM sheephead and pooled PISCO sizes.
10. **Missing at random** given station, year, region, month.
11. **Stable survey footprint within station × zone.**
12. **No spillover** to references (if present, effects are biased toward zero).
13. **Counts are complete counts of the stated area.** Unverified for PISCO swaths, where the field protocol may allow subsampling dense species. If the delivered counts were scaled up, use the searched area and raw count.
14. **Each protected station represents its MPA.** With one station per MPA this is an assumption, not a result.
15. **Allometry form** holds across sites and years.

## 7. Identifiability

- **Curve vs calendar.** g(a), γ, χ and ω all vary in time. g is identified only by protected-minus-reference divergence that is smooth in age and common to MPAs; ω takes year-specific common divergence; γ and χ take everything shared with references. Two designation cohorts separate age from year for ages 0–11 only.
- **One protected station.** D[m] and that station's AR path trade off. ρ and σ_δ are learned from all ~120 stations, mostly references.
- **AR vs white noise.** Separated by long series (KFM, LTER). PISCO's short, noisy histories borrow ρ from them; this assumes the same population process across programs.
- **Risk of over-building.** v2 has more variance terms than v1. Each of χ, ω, ζ and month effects is added one at a time in stage 1 and kept only if posterior predictive checks or leave-one-year-out prediction improve. A term whose SD posterior matches its prior is dropped and reported as not identified.

## 8. Validation and pre-specification

1. **Freeze first.** Commit this specification, the estimands, the predicted directions and the reporting rule (below) before any fit on real protection labels. We have already seen two-stage results, so the freeze limits further tuning; it does not make the analysis blind.
2. **Simulation.** Parameter recovery and interval coverage on the real station × year × MPA layout, at realistic truth (effects of 0, ×1.5 and ×2.5 per decade). This also gives the power of the design directly.
3. **Placebo in time.** KFM 1982–2002 at Harris Point, Santa Barbara Island and Scorpion (18–21 years, four taxa): assign a fake designation in 1992, fit the same model. Effects should centre on zero. Three MPAs is a small test, but it is the only direct check of parallel trends.
4. **Placebo in space.** Relabel a reference station as protected within each cluster (where ≥3 references exist). The pooled effect should be zero; the spread of placebo estimates calibrates the false-positive rate.
5. **Posterior predictive checks:** zero fractions, visit totals, dispersion, mass distributions, protected–reference trajectories, residual cross-MPA correlation (should vanish once ω is in).
6. **Agreement with the integrated review** on the same data, with a linear g for comparability. Explain every disagreement.
7. **Sensitivities:** τ_D, LKJ, σ_g scales; pre-trend; post-only λ rule; all zones; ≥25 mm density; β ± 0.1; reference omissions (Harris Point, Webster's Arch, Wyckoff); direct-size-only biomass; leave-one-MPA-out.

**Reporting rule (proposed).** Smallest effect of interest: a 25% change per decade (|E10| > 0.22). For each taxon and quantity: *supported* if Pr(E10 beyond 0.22 in the predicted direction) ≥ 0.90; *practically absent* if Pr(|E10| < 0.22) ≥ 0.90; otherwise *inconclusive*. This separates "no effect" from "can't tell", which significance testing cannot do. The thresholds are a judgment call for Adrian and coauthors.

**Computation.** Stan (cmdstanr or cmdstanpy) for the full model. Stage 1 is also fitted in glmmTMB or TMB (negative binomial, AR(1), random effects) as a fast independent implementation; the simulations in step 2 need hundreds of refits, which only the fast version makes practical. About 4,700 visit-level count observations for retained density data. Run time is unverified.

## 9. Decisions for Adrian

1. **Mean mass without measured animals:** pool from the same station's other years and the region-year (recommended, with the measured-information share reported), same station only, or measured station-years only.
2. **Regions and clusters:** proposed regions are Northern Channel Islands, Santa Barbara Island, Santa Barbara mainland, Malibu, Palos Verdes, Catalina, San Diego.
3. **Smallest effect of interest and the 0.90 thresholds** in the reporting rule.
4. **Flat baseline as main, pre-trend as sensitivity** (current policy).
5. **Dynamic cascade link (stage 2b):** in the paper, or supplementary only?
6. **Is density per program acceptable as the main definition,** with ≥25 mm as a sensitivity?

## 10. Build plan

| Stage | Scope | Gate |
|---|---|---|
| 0 | Freeze spec; data builder; verify PISCO counting protocol and LTER urchin size threshold; fill take-rule blanks; simulation harness | Simulated recovery and coverage acceptable |
| 1 | Purple urchin density. (a) linear g, compare with −9.7%/yr [−16.9, −1.7]; (b) response curve; (c) add χ, ω, zones, month one at a time | Placebo tests pass; each added term earns its place |
| 2 | All taxa, density, Ω. 2b: dynamic links | Prior sensitivity on τ_D, LKJ |
| 3 | Mean mass and biomass for urchins, lobster, sheephead | Decision 1; Gamma check |
| 4 | Kelp: diver and Landsat | Footprint comparability |

Code: `code/joint_model/`. Outputs: `outputs/joint_model_<date>/`. Station selection and exclusions come unchanged from the integrated review.
