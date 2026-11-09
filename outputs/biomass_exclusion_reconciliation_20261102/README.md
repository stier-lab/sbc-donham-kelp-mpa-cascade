# KFM donor exclusions and later-source check

The size-site discrepancies are present in the2015NPS release as well as the1985–2011legacy data. All four2015comparisons reproduce the complete published grouped frequencies and sample totals: a-k-01contains Wyckoff Ledge1 plus Miracle Mile21; a-k-14matches Southeast Reef37. The2015tables contain grouped sizes, so they are used to verify identity, not to invent individual measurements or nonlinear mean masses. Published Mean_size fields were not used for biomass conversion.

## Why this matters

Emily's original pBACIPS_PISCO_V10.R explicitly removes sites21and37 from SizeFreq.south. Current03_data_import.R/01_utils.R applies the same site exclusion list. The combined source labels conceal those original sites, so exclusion by delivered label does not enforce the intended physical-site exclusions. This is an input/implementation mismatch, not a new exclusion chosen from results. Its origin in the upstream harmonization remains unlocated: the checked original script reads the combined file; it does not create it.

Six retained biomass histories use these labels: Harris Point red/purpleurchins in KFM and PISCO, and Santa Barbara Island red/purpleurchins in KFM. Cross-program reuse carries the Harris Point discrepancy into PISCO. Suspect labels occur through2022, but physical identity is directly checked only where original records are available; do not assert every later record has been identified.

## Same observations, intended donor exclusions

For years through2011 with exact source matches, replace the grouped Wyckoff/Miracle Mile measurements with original Wyckoff measurements and the mislabeled Southeast Reef measurements with original Southeast Sea Lion measurements. Other size donors remain unchanged. PISCO retains its25mm size rule. Counts, retained years and previously verified offsets stay fixed. Do not split combined frequency rows arbitrarily.

Across the full observed post-protection histories, this partial reconciliation changes Harris Point red slopes from−0.150to−0.160(KFM) and−0.041to−0.049(PISCO); purple from0.023to0.003(KFM) and0.011to−0.016(PISCO). All four conditional intervals still include zero. Santa Barbara Island red changes−0.123to−0.006 and purple−0.068to−0.008, reproducing the previous check. Units are annual change in log protected/reference biomass ratio.

These are partial diagnostic estimates: only seven or eight historical years per history were corrected. Remaining years retain unresolved inputs. The effects table also includes fits restricted to reconciled years, keeping those distinct from full-history partial replacements. No final primary estimate or pooled conclusion is updated.

## Recent source located

NPS package2318436 (1982–2025; DOI10.57830/2318436) lists natural-habitat size, site information, survey schedule, quadrat and kelp-size tables. Its file access isInternal and CanDownloadfalse. No restricted file download was attempted. The associated public processing script2316006 was downloaded and inspected; it preserves original event/site identifiers and Size_mm but requires the underlying Access database to run. It was not executed. These are acquisition leads, not replacement inputs. Scott Gabara is a listed package author. Requesting the original size/site/survey tables through2022 would allow the existing analysis window to remain unchanged.

## What can now be resolved and what cannot

- The intended exclusion rule is clear in the original code and is unchanged.
- Historical allowed-site sizes are recoverable;24refits independently reproduce in R.
- The2015source confirms persistence but cannot supply exact individual size distributions for nonlinear conversion.
- Recent original site-specific size records and the script that assembled ALL_sizefreq_2024.csv are still needed to rebuild all affected years consistently.
- Density remains available. The six biomass histories remain provisional; do not silently zero or drop their uncertain years.

The main inventory and assumption ledger record this distinction. Raw data and manuscript/frozen results are unchanged. The old source documentation's claim that the combined file lacks KFM/LTER size records is stale; source documentation is corrected with a dated current-state note.
