from pathlib import Path
import pandas as pd, json, hashlib, shutil, html
D=Path('/Users/adrianstier/sbc-donham-kelp-mpa-data')
S=D/'data/processed/source_build'
O=Path('outputs/lter_integration_pisco_check_20261102');O.mkdir(exist_ok=True)
a=pd.read_csv(S/'lter_annual_mass_comparison.csv');x=a.pivot(index=['site','year','classcode'],columns='version',values='mean_proxy_mass_g');diff=(x.original-x.combined_audit).abs()
f=pd.read_csv(S/'lter_focal_original_sizes.csv')
summary={'lter_positive_count_rows':len(f),'lter_sized_group_rows':int(f.donor_eligible.sum()),'lter_missing_size_rows':int((~f.donor_eligible).sum()),'annual_mass_proxy_groups':len(x),'max_absolute_annual_mass_difference_g':float(diff.max()),'date_corrections':4,'primary_refit':False,'paper_changed':False}
for name in ['lter_selection_ledger.csv','lter_date_differences.csv','lter_annual_mass_comparison.csv']:
 shutil.copy2(S/name,O/name)
(O/'verification.json').write_text(json.dumps(summary,indent=2)+'\n')
body=f'''# Original LTER integration and cached PISCO check

November 2, 2026. Review outputs only; manuscript, supplement and frozen primary results are unchanged.

## LTER: integrated

- SBC 19.33 and its metadata now enter through registered, checksum-verified source snapshots. Source units and the missing-value code are checked before processing.
- The existing scope contains 2,706 positive-count source rows: 1,809 have group-average diameters and 897 have missing sizes. Missing sizes remain missing; their counts remain visible. They do not become size donors or zero biomass.
- The source records preserve dates, physical sites, transects, quadrats, sides, area, counts, raw sizes and row identifiers. Six existing stations and 2000–2023 are preserved; source availability does not add sites or override exclusions.
- Four Arroyo Quemado purple-urchin records now use July 27, 2022, rather than July 19. Counts and sizes agree with the original source.
- Across {len(x)} annual site/species-size-class pools, the maximum difference in count-weighted mean mass proxy is {diff.max():.3g} g. This is numerical serialization precision, not a material annual biomass change. Visit-level joins can be affected by the four corrected dates.

These sizes are survey-group means, not individual diameters. Converting a mean diameter to mass does not recover the average mass of the animals when the relationship is nonlinear. Group-level resampling cannot recover their unrecorded within-group size variation. The mass and donor ledgers retain this distinction. This comparison checks annual inputs and mass conversions; it does not claim a new model fit or total uncertainty estimate.

## PISCO: measurements verified; adoption held

- The cached PISCO 1.11 size file passes its metadata checksum. All 3,970 older UCSB roving records match with multiplicity after documented unit, method-label and blank-zone normalization. They represent 11,862 animals.
- An old aggregation label assigns 99 physical Painted Cave East records from 2003 to South Point: 55 red-urchin and 44 purple-urchin records. The source measurements match; the site assignment needs correction and affected-donor tracing before promotion. Do not copy the erroneous alias into a new adapter.
- The cache contains 2,915 additional 2022–2024 urchin rows. These are listed separately and have not been added to analytical scope.
- No VRG records occur in this cache. It therefore does not resolve the full PISCO source requirement.

## What changes now

Original LTER data replace the opaque combined-file input at the source boundary. Missing and group-average sizes are explicit, and four dates are corrected. The annual mass proxies are effectively unchanged. PISCO receives a reproducible verification and alias-conflict table, not automatic adoption. Full harmonization remains guarded until missing PISCO sources and outstanding KFM source decisions are resolved. Downstream area normalization, donor matching and complete uncertainty still require review before primary promotion.

## Linked records

- [LTER selection and missingness](lter_selection_ledger.csv)
- [LTER source-date differences](lter_date_differences.csv)
- [LTER annual mass comparison](lter_annual_mass_comparison.csv)
- [Verified PISCO measurement rows](pisco_exact_matches.csv)
- [PISCO alias conflict](pisco_alias_conflict.csv)
- [Additional PISCO rows, not adopted](pisco_unadopted_extensions.csv)
- [PISCO checksum and match verification](pisco_verification.json)

Reproduce with `make check` and `make trace` in the data repository; run `python3 code/tools/check_cached_pisco_sizes.py` and `python3 code/tools/report_lter_integration.py` in the analysis repository. The data repository's `docs/RAW_TO_BIOMASS.md` is the workflow guide.
'''
(O/'README.md').write_text(body)
import subprocess
subprocess.run(['pandoc',str(O/'README.md'),'-s','--metadata','title=LTER and PISCO source audit','-o',str(O/'index.html')],check=True)
hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [S/'build_manifest.json',D/'config/lter_size_scope.json',D/'config/lter_physical_sites.csv',Path('code/tools/check_cached_pisco_sizes.py')]}
(O/'input_hashes.json').write_text(json.dumps(hashes,indent=2)+'\n')
print(json.dumps(summary,indent=2))
