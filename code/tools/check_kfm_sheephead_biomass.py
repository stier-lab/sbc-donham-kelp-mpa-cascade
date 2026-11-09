from pathlib import Path
import json,hashlib
import pandas as pd,numpy as np
A=Path(__file__).resolve().parents[2];O=A/'outputs/kfm_sheephead_biomass_20261103'
for row in json.loads((O/'inputs.json').read_text()):assert hashlib.sha256(Path(row['path']).read_bytes()).hexdigest()==row['sha256']
r=pd.read_csv(A/'outputs/fish_lter_reconciliation_20261102/fish_count_size_matches.csv');v=pd.read_csv(O/'visit_support.csv');a=pd.read_csv(O/'annual_responses.csv');years=pd.read_csv(O/'year_eligibility.csv');t=pd.read_csv(O/'trends.csv');s=pd.read_csv(O/'size_record_use.csv');cal=pd.read_csv(O/'variance_calibration_membership.csv')
allowed=set(map(tuple,r[['site_id','year']].to_numpy()))
assert set(map(tuple,cal[['site_id','year']].to_numpy()))<=allowed
assert set(map(tuple,s[['site_id','year']].to_numpy()))<=allowed
assert not s.source_row.duplicated().any()
ledger=pd.read_csv(O/'site_year_support.csv')
for row in ledger.itertuples():
 used=s[s.site_id.eq(row.site_id)&s.year.eq(row.year)]
 assert row.n_size_rows==len(used) and row.n_fish==used.n_fish.sum()
 assert row.n_size_visits==used.date.nunique()
assert years.eligible.sum()==52
assert a.groupby('MPA').year.nunique().to_dict()=={'Anacapa Island SMR 2003':13,'Harris Point SMR':13,'Scorpion SMR':12,'South Point SMR':14}
for row in a.itertuples():
 sites=row.sites.split(';');c=r[r.MPA.eq(row.MPA)&r.year.eq(row.year)&r.status.eq(row.status)]
 assert set(sites)==set(c.site_id)
 assert np.isclose(c.count_density.mean(),row.density)
 z=v[v.site_id.isin(sites)&v.year.eq(row.year)]
 assert not z['mode'].eq('unsupported').any()
 assert (abs((pd.to_datetime(z.count_date)-pd.to_datetime(z.size_date)).dt.days)<=7).all()
 site_b=(z.biomass_g_m2*z.annual_weight).groupby(z.site_id).sum()
 assert np.isclose(site_b.mean(),row.biomass_g_m2)
 assert row.biomass_upper>=row.biomass_lower>=0
 assert row.biomass_upper>0 # no exact certainty at observed zero
for variant in ['empirical','pool_strength4','count_prior01','dependence_positive','dependence_negative']:
 x=pd.read_csv(O/variant/'annual_responses.csv');assert np.allclose(x[['density','biomass_g_m2']],a[['density','biomass_g_m2']])
 z=pd.read_csv(O/variant/'trends.csv');assert np.isfinite(z.slope).all()
assert np.isfinite(t.slope).all() and not t.upper_boundary.fillna(False).any()
f=pd.read_csv(O/'equal_year_slopes.csv');u=f[f.response.eq('biomass')].set_index('MPA').slope;w=f[f.response.eq('biomass_count_only')].set_index('MPA').slope;assert np.allclose(u,w)
# Existing researcher exclusions and production files stay unchanged.
protected=json.loads((A/'outputs/researcher_scope_20261101/protected_hashes.json').read_text())
for path,digest in protected.items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest,path
result=dict(protected_files_unchanged=len(protected),input_hashes_verified=True,eligible_comparison_years=52,visit_products_verified=True,calibration_restricted_to_roster=True,zero_count_uncertainty_preserved=True,variant_point_estimates_identical=True)
(O/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
