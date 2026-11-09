# Biomass assumptions and source reconciliation

The assumptions/change ledger separates confirmed source facts, adopted corrections, diagnostic assumptions and unresolved decisions. It is maintained in docs/BIOMASS_ASSUMPTIONS_AND_CHANGES.csv. Earlier results remain archived; superseded assumptions are identified rather than silently erased.

## What the source checks resolved

The original KFM metadata confirms species11005 is red urchin,11006 is purple urchin, and legacy test diameters are in centimeters. Matching complete size-frequency distributions after conversion to millimeters gives873 exact matches among904 overlapping site/year/species groups. Dates can be recovered for58,369 frequency records representing159,963 animals. These are frequency records, not58,369 independent dives or newly collected measurements. In retained donor pools,26,015 unique original records benefit from the overlay. Raw source data remain unchanged.

Recovered dates do not prove same-visit sampling. All873 matched groups have count dates in MBON;401 have size dates that match a count date exactly. The largest nearest-date gap is219days. Among511 recovered donor-site/year/species uses in retained histories,193 have a nonzero gap. These are donor-site timing comparisons, not confirmation that borrowed sizes coincide with the target count sites. No cutoff or automatic exclusion was applied.

## Two site labels require special treatment

**Wyckoff Ledge:** all51 overlapping species/year groups labeled a-k-01 match the combined legacy distributions from Wyckoff Ledge1 and Miracle Mile21. Fifteen disagree with Wyckoff alone; the other36 predate/add no Miracle Mile measurements. Both sites are references for Harris Point. This supports a grouped source identity, not a license to assign every animal to Wyckoff or divide aggregated frequencies arbitrarily.

**Santa Barbara Island:** all14 red/purple species-year distributions labeled a-k-14 in2005–2011 exactly match Southeast Reef37. The corresponding MBON count densities match Southeast Sea Lion14 in all14checks. Thus the counts and sizes represent different sites. The previous audit correctly matched the delivered labels but its physical same-site interpretation was wrong. This finding supersedes that interpretation and qualifies the earlier count-site-only donor sensitivities. Later years have not been source-reconciled.

## What changes in an isolated same-observation test

Using original Southeast Sea Lion sizes for2005–2011, with existing offsets fixed, changes the17-year red-urchin recovery slope from−0.123 [−0.193,−0.052] to−0.006 [−0.041,0.028] log-ratio units per year. Purple changes from−0.068 [−0.202,0.066] to−0.008 [−0.160,0.145]. These intervals describe temporal model uncertainty conditional on the annual inputs, not total biomass uncertainty.

This is a partially reconciled diagnostic: only2005–2011 sizes were replaced; the remaining years retain current inputs. It is not a corrected final17-year estimate. The table also reports the seven-year matched comparison separately. The result shows that source-site identity matters and must be resolved before interpreting this biomass history. It does not justify choosing one source because of the direction of its slope. Density is unchanged.

## Protocol boundaries

The PISCO technical report supports a roughly2.5cm minimum for conspicuous benthic swath organisms. Its greater-than wording does not independently settle whether exactly25mm is included; the current>=25mm implementation is retained and its boundary remains explicit. KFM size protocols include small animals and sample across the permanent transect area. Changing only the size threshold would not standardize counts. LTER-specific count/size protocol reconciliation remains open.

## Next checks

- Resolve a-k-14 and other grouped site labels in later size releases and the original harmonization code.
- Recover recent KFM survey dates and assess count-size timing without selecting a cutoff by outcome.
- Recompute physical count-site support using verified donor identities; earlier label-based percentages are not final physical-coverage estimates.
- Complete joint count/size uncertainty and KFM kelp plant-sampling checks.

No raw data, manuscript or frozen primary result was changed in this pass. Eight new model fits were independently verified in R; the14count identities and14size identities are exact within numerical precision. The source date overlay remains private; public tables contain aggregated checks.

## Sources

[Kushner et al. KFM metadata](https://esapubs.org/archive/ecol/E094/245/metadata.php), [size variable definitions](https://esapubs.org/archive/ecol/E094/245/Table13B_invertebrate_size_frequency_variables.html), and [PISCO technical report](https://caseagrant.ucsd.edu/sites/default/files/R_MPA-43_Kelp%20Forest%20Technical%20Report%20Narrative.pdf). Downloaded source copies and record-level overlays are retained privately with checksums.
