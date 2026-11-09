#!/usr/bin/env python3
"""Refresh isolated 16-pair inputs; never chooses sites by fitted response."""
from pathlib import Path
import pandas as pd
out=Path('outputs/satellite_baseline_comparison_20261101')
# Frozen scope copied from the pre-screen 145-history cohort, not current primary eligibility.
co=pd.read_csv(out/'cohort.csv');assert len(co)==16 and co.MPA.is_unique
d=pd.read_csv('data/harmonized/harmonized_landsat_rr.csv').rename(columns={'CA_MPA_Name_Short':'MPA'}).merge(co,on='MPA',validate='many_to_one')
assert not d.duplicated(['MPA','year']).any()
d.to_csv(out/'annual_input.csv',index=False)
# Freeze external metadata with this review analysis; source paths/hashes in source_manifest.json.
prot=pd.read_csv(out/'protection_metadata.csv').set_index('MPA')
rows=[]
for name,g in d.groupby('MPA'):
 b=g[g.BA=='Before'];a=g[g.BA=='After'];start=int(g.MPA_Start.iloc[0]);missing=sorted(set(range(g.year.min(),g.year.max()+1))-set(g.year))
 assert g.lnDiff.notna().all() and len(b)>=3 and len(a)>=3
 rows.append(dict(MPA=name,designation_year=start,protection_level=prot.loc[name,'Protection'],protection_exception=prot.loc[name,'exceptions'],protection_source=prot.loc[name,'source'],first_year=g.year.min(),last_year=g.year.max(),paired_years=len(g),missing_years=len(missing),missing_year_list=','.join(map(str,missing)),before_years=len(b),strict_before_years=int((g.year<start).sum()),after_years=len(a),designation_year_label='Before',reference_id='supplied satellite footprint: '+name,reference_matching='Geometric assignment documented; ecological equivalence not established',overlap='320 reference pixels shared with '+('Naples' if name=='Campus Point SMCA' else 'Campus Point') if name in ['Campus Point SMCA','Naples SMCA'] else 'See footprint overlap audit; absence of overlap not assumed',data_check='Unique finite annual paired log ratios',review_eligibility='Include in requested 16-pair comparison; not final ecological validation',current_primary='Excluded by baseline trend screen' if name in ['Santa Barbara Island SMR','Swamis SMCA'] else 'Included',evidence='docs/SOURCE_RECONCILIATION_2026-10-17.md; data/external/tom_bell_2026-10-20/; tables/table_s_reference_assignment_audit.csv'))
pd.DataFrame(rows).to_csv(out/'eligibility.csv',index=False)
