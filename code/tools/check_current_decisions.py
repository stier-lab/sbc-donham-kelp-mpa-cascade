from pathlib import Path
import pandas as pd
import hashlib,json
O=Path('outputs/dataset_inventory_20261101')
c=pd.read_csv(O/'comparison_decisions.csv');m=pd.read_csv(O/'audit_master.csv');r=pd.read_csv(O/'current_rules.csv');h=pd.read_csv(O/'decision_history.csv');v=pd.read_csv(O/'measurement_visit_decisions.csv');b=pd.read_csv(O/'measurement_support_current.csv')
assert len(c)==206 and c.series_id.is_unique and c.exclude.sum()==31
pd.testing.assert_frame_equal(c[['series_id','exclude','reason']],m[['series_id','exclude','reason']])
assert set(r.rule_id)=={f'R{i:02d}' for i in range(1,16)} and r.rule_id.is_unique
assert all(Path(x).exists() for x in r.evidence_path)
assert all(set(x.split(';'))<=set(r.rule_id) for x in c.rule_ids)
assert len(b)==94 and b.series_id.is_unique
expected_history=51 if Path('outputs/kfm_scott_reconciliation_20261102/verification.json').exists() else 49 if Path('outputs/nonlinear_biomass_calibration_20261102/report_verified.json').exists() else 46 if Path('outputs/joint_borrowed_biomass_20261102/report_verified.json').exists() else 43
assert len(h)==expected_history and h.id.is_unique
assert h.set_index('id').loc['B33','lifecycle']=='Earlier stage; see current interpretation'
source=pd.read_csv('outputs/pisco_biomass_borrowing_20261102/visit_biomass_support.csv')
pd.testing.assert_frame_equal(v,source)
assert v.method.value_counts().to_dict()=={'measured_same_visit':624,'zero_count':285,'borrowed_other_site':131,'unsupported':38,'borrowed_other_visit':8}
assert not ((v.method.str.startswith('borrowed')) & v.status.ne('reference')).any()
s=pd.read_csv('outputs/pisco_biomass_borrowing_20261102/summary.csv').set_index('series_id')
bs=b.set_index('series_id').loc[s.index]
for col in ['original_post_years','local_only_post_years','borrowing_post_years']:
    assert (bs[col]==s[col]).all()
assert bs.borrowing_post_years.sum()==44 and bs.original_post_years.sum()==63
for path,sha in json.loads((O/'current_decision_inputs.json').read_text()).items():
    assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==sha,path
print(f'Current rules, 206 eligibility rows, 94 biomass histories, 1086 visit uses, {len(h)} history entries and input hashes verified.')

assert not c.loc[c.Source.eq('Landsat'),'rule_ids'].str.contains('R03').any()
assert c.loc[c.Source.eq('PISCO') & c.Taxa.eq('S. pulcher') & ~c.exclude,'rule_ids'].str.contains('R05').all()
six_den=set(s.index.str.replace(r'\|Bio$', '|Den',regex=True))
assert c.loc[c.series_id.isin(six_den),'rule_ids'].str.contains('R04').all()

next_steps=Path('outputs/biomass_next_steps_20261102')
if (next_steps/'completed.json').exists():
    density=pd.read_csv(next_steps/'current_density_annual.csv')
    assert len(density)==84 and (density.t>0).sum()==63
    assert density.loc[density.t>0,'ln_ratio'].notna().all()
    assert density.ln_ratio.isna().sum()==3 and density.full_zone_evidence_lnratio.notna().all()
    assert pd.read_csv(O/'donor_attempts_current.csv').shape[0]==1416
    assert pd.read_csv(O/'measurement_visit_support_detailed.csv').shape[0]==1086
    assert c.loc[c.series_id.isin(set(density.series_id)),'roster_basis'].str.contains('rebuilt').all()
