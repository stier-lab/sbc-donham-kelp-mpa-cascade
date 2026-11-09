"""Scientific safeguards for the frozen six-history review build."""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
O=Path('outputs/pisco_biomass_borrowing_20261102');v=pd.read_csv(O/'visit_biomass_support.csv');a=pd.read_csv(O/'biomass_annual.csv');d=pd.read_csv(O/'selected_density_annual.csv');old=pd.read_csv('outputs/pisco_all_months_20261102/annual_values.csv');old=old[(old.scenario=='allmonths_sharedzones_visit_sizes')&old.series_id.str.endswith('|Den')];pd.testing.assert_frame_equal(d.reset_index(drop=True),old.reset_index(drop=True))
assert v.loc[v.method=='zero_count','biomass_point'].eq(0).all();assert v.loc[v.method=='unsupported','biomass_point'].isna().all()
for r in a.itertuples():assert r.supported == (not ((v.series_id==r.series_id)&(v.year==r.year)&(v.method=='unsupported')).any())
c=pd.read_csv(O/'outer_site_calibration.csv');assert (c.outer_site!=c.cal_site).all();assert all(r.outer_site not in r.donor_sites.split(';') for r in c.itertuples())
x=pd.read_csv(O/'private/visit_draw_checks.csv');x=x.merge(v[['series_id','site','year','zone','date','predicted']],on=['series_id','site','year','zone','date'],validate='one_to_one');x=x[x.method.str.startswith('borrowed')]
for draw in range(3):
 x['factor']=x['draw_'+str(draw)]/x.density/x.predicted
 for _,z in x.groupby(['series_id','site','year','method']):assert np.allclose(z.factor,z.factor.iloc[0])
protected=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text());assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in protected.items())
print('Density preserved; zero/unknown guards, full reference membership, outer-site holdouts, shared error draws and protected hashes pass')
