# Fixed monthly kelp-coefficient refit — 29 October 2026

The current analysis uses SBC 127.5, Macrocystis biomass relationships, model 1 (Nelson et al. 2026; DOI 10.6073/pasta/3a6272c195d7de89eb412cf1f0c4d63e). The provider CSV MD5 is `8ddaf77458d21adfb6d0398cc1c51fda`. Units remain kg dry foliar mass per frond, multiplied by 1000 for grams. Each recorded survey month is converted before annual aggregation; KFM resamples dated plant records using the existing seeds and matching rules.

## Scope and checks

- Rebuilt all source processors into an isolated candidate folder before promotion.
- Only kelp biomass was promoted: 531 of 6280 raw-response rows and 354 of 3341 paired-response rows changed. Non-kelp rows were preserved byte-for-byte; density, satellite files and membership did not change.
- All 145 primary history identities remain. Seventeen kelp-biomass contrasts changed; no non-kelp primary contrast changed.
- Two PISCO histories changed selection from mean to linear: Naples (new slope p = 0.04961) and South Point (p = 0.04707). Both cross the existing 0.05 rule.
- The joint biomass model shares MPA and program variance among taxa. Thus pooled non-kelp estimates can change without changes to their observations or individual contrasts.

## Pooled estimates

Exponentiated contrasts combine different estimands; these are not uniform biomass ratios.

| Taxon | Response | Previous estimate | Current estimate | Current 95% CI | Current BH p |
|---|---|---:|---:|---|---:|
| S. purpuratus | Biomass | 0.668 | 0.498 | 0.307–0.806 | 0.011518 |
| M. franciscanus | Biomass | 2.057 | 1.690 | 1.080–2.644 | 0.039985 |
| M. pyrifera | Biomass | 2.160 | 3.096 | 1.940–4.941 | 0.000067 |
| P. interruptus | Biomass | 1.657 | 1.442 | 0.810–2.570 | 0.270392 |
| S. pulcher | Biomass | 3.216 | 2.612 | 1.650–4.134 | 0.000375 |
| S. purpuratus | Density | 0.137 | 0.137 | 0.044–0.426 | 0.002633 |
| M. franciscanus | Density | 0.543 | 0.543 | 0.176–1.676 | 0.318256 |
| P. interruptus | Density | 2.590 | 2.590 | 0.838–8.001 | 0.145195 |
| S. pulcher | Density | 0.727 | 0.727 | 0.236–2.239 | 0.572539 |

## Interpretation

All nine pooled directions persist. Purple-urchin biomass now passes the nine-test BH correction and all calendar-HAC specifications. That result is conditional on contrast selection. Substituting mean contrasts for the two near-threshold PISCO kelp trends keeps the kelp result positive (log estimate 0.7725, 95% CI 0.3687–1.1763), but leaves purple-urchin biomass uncertain (−0.4054, −0.8271–0.0162; BH p = 0.0946). This uses current measurements throughout. It is a targeted sensitivity, not a substitute primary analysis or an estimate of full selection uncertainty. Module 51 reproduces it in the main pipeline.

During the heatwave the kelp baseline-relative advantage is 4.24 (95% CI 1.37–13.12). Inside/reference median fold changes are 2.3/0.9, 2.3/0.8 and 2.2/0.6 across the three periods; all focal adjusted signed-rank p values remain 0.014. Historical satellite results are unchanged. These results do not establish a predator-mediated protection effect on resilience.

## Reproduction

`code/promote_fixed_kelp.py` validates and promotes only the scoped biomass inputs. `data/harmonized/kelp_fixed_calibration.json` extends the existing verified release chain. `code/R/run_all.R` runs the primary and supporting analyses, including the current-input model-selection check in `51_kelp_selection_sensitivity.R`. The optional resilience modules 17, 18 and 20 were also rerun. `code/audit_fixed_kelp_refit.R BASELINE_DIRECTORY` compares a frozen preceding release and reproduces the selection check. Results are recorded in `tables/table_s_fixed_kelp_*.csv`.

The preceding files are frozen locally under `local_archive/kelp1275_2026-10-29/before`; source bytes for the previous calibration remain in `data/external/sbc11215/`. Historical audit records describe their dated releases and should not be read as current results. Li Kui confirmed on 29 October 2026 that SBC127.5, Table 2, Model 1 is the intended fixed calibration for future SBC kelp biomass calculations. It is the exact source used here; no further numerical rerun is required by her clarification.
