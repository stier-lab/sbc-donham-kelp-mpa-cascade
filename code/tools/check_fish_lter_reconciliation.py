"""Verify source hashes, unchanged primary outputs, and comparison invariants."""
from pathlib import Path
import json,hashlib
from collections import Counter
import pandas as pd,numpy as np
A=Path(__file__).resolve().parents[2];D=A.parent/'sbc-donham-kelp-mpa-data';O=A/'outputs/fish_lter_reconciliation_20261102';P=D/'data/source_archive/fish_lter_reconciliation_20261102'
for r in json.loads((O/'input_manifest.json').read_text()):assert hashlib.sha256((P/r['file']).read_bytes()).hexdigest()==r['sha256'],r['file']
protected=json.loads((A/'outputs/researcher_scope_20261101/protected_hashes.json').read_text())
for f,h in protected.items():assert hashlib.sha256(Path(f).read_bytes()).hexdigest()==h,f
reg={s['source_id']:s for s in json.loads((D/'config/source_manifest.json').read_text())['sources']}
k=['YEAR','MONTH','DATE','SITE','TRANSECT','QUAD','SIDE','AREA','FRONDS'];old=pd.read_csv(D/reg['lter_kelp']['path']);new=pd.read_csv(P/'Annual_Kelp_All_Years_20260923.csv');new=new[new.YEAR.le(2023)]
assert Counter(map(tuple,old[k].fillna('').to_numpy()))==Counter(map(tuple,new[k].fillna('').to_numpy()))
f=pd.read_csv(O/'fish_count_size_matches.csv');assert len(f)==423 and f.new_has_sizes.sum()==282 and f.all_size_dates_match.sum()==279
p=pd.read_csv(O/'fish_size_pools.csv');assert p.year.min()==2007 and p.year.max()==2025
x=pd.read_csv(O/'lter_transect_comparison.csv');assert not x[['provider_dry_g_m2','provider_density']].eq(-99999).any().any();assert x._merge.eq('both').sum()==968 and x._merge.eq('right_only').sum()==132
z=pd.read_csv(O/'lter_formula_decomposition.csv');finite=z.density_residual.dropna();assert finite.abs().max()<1e-10
r=pd.read_csv(O/'lter_role_year_comparison.csv');assert r.loc[r.missing>0,'provider'].isna().all()
q=pd.read_csv(A/'outputs/dataset_inventory_20261101/unresolved_issues.csv').set_index('issue_id');assert '282' in q.loc['U07','finding'] and 'SBC50.21' in q.loc['U16','finding']
summary={'input_files_verified':len(json.loads((O/'input_manifest.json').read_text())),'adult_raw_rows_unchanged':len(old),'protected_files_unchanged':len(protected),'fish_count_site_years':len(f),'fish_size_supported_site_years':int(f.new_has_sizes.sum()),'lter_transects_matched':968,'missing_values_preserved':True,'primary_or_paper_changed':False}
(O/'verification.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
