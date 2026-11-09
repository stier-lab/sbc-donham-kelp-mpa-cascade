from pathlib import Path
import pandas as pd,numpy as np,hashlib,json
A=Path(__file__).resolve().parents[2];O=A/'outputs/kfm_missing_size_recovery_20261103'
manifest=json.loads((O/'input_manifest.json').read_text())
for r in manifest:assert hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()==r['sha256'],r['path']
e=pd.read_csv(O/'missing_visit_evidence.csv');f=pd.read_csv(O/'five_year_options.csv');c=pd.read_csv(O/'all_year_coverage_options.csv');w=pd.read_csv(O/'matched_visit_subset_weights.csv');s=pd.read_csv(O/'matched_visit_subset_responses.csv')
assert len(e)==7 and len(f)==5 and e.recovered_measurements.sum()==0
assert c.all_visits_supported.sum()==52 and c.matched_visit_subset_possible.sum()==55
assert f.matched_visit_subset_possible.sum()==3
assert e.scott_sheephead_rows_on_date.eq(0).all()
assert e.loc[e.year.le(2011),'legacy_sheephead_rows_on_date'].eq(0).all()
assert np.allclose(w.groupby(['MPA','year','site_id']).new_weight.sum(),1)
assert pd.read_csv(O/'drive_copies.csv').matches_current.all()
assert (f.n_retained_count_rows<=f.n_original_count_rows).all()
protected=json.loads((A/'outputs/researcher_scope_20261101/protected_hashes.json').read_text())
for path,digest in protected.items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest,path
result=dict(source_hashes_verified=len(manifest),missing_visits=7,comparison_years=5,raw_sizes_recovered=0,current_review_years=52,candidate_matched_visit_years=55,protected_files_unchanged=len(protected),primary_changed=False)
(O/'verification.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
