"""Independent keyed checks of review scope, mass corrections and current guide."""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
O=Path('outputs/kfm_scott_reconciliation_20261102')
a=pd.read_csv('outputs/biomass_resolution_20261102/review_annual_input.csv');v=pd.read_csv(O/'review_annual_input.csv');r=pd.read_csv(O/'all_corrections.csv');c=pd.read_csv(O/'recent_corrections.csv');h=pd.read_csv(O/'kfm_current_holds.csv')
assert len(c)==54 and len(r)==99 and len(h)==100 and c.histogram_absolute_difference.eq(0).all()
assert c.allowed_source_site.isin(['Wyckoff Ledge','SE Sea Lion Rookery']).all()
assert c.max_source_count_gap_days.le(1).all() and c.allowed_replacement_records.gt(0).all()
assert (c.cross_program_borrowing.sum()==14) and c.sparse_lt10.sum()==1
assert c.loc[c.sparse_lt10,'allowed_replacement_records'].iloc[0]==7
pd.testing.assert_frame_equal(a.drop(columns='lnDiff'),v.drop(columns='lnDiff'))
keys=['series_id','year'];old=a.set_index(keys);new=v.set_index(keys);rc=r.set_index(keys)
assert np.allclose(new.loc[rc.index,'lnDiff'],rc.new_lnDiff)
assert np.allclose(old.loc[rc.index,'lnDiff'],rc.old_lnDiff)
untouched=~old.index.isin(rc.index);pd.testing.assert_frame_equal(old.loc[untouched],new.loc[untouched])
geo=pd.read_csv('data/MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv').set_index('site_id');cross=[]
for sid,name,allowed in [('a-k-01','Wyckoff Ledge',True),('a-k-21','Miracle Mile',False),('a-k-14','SE Sea Lion Rookery',True),('a-k-37','Southeast Reef',False)]:
 cross.append(dict(site_id=sid,scott_label=name,geolocation_label=geo.loc[sid,'geolocation'],allowed_donor=allowed,evidence='MBON geolocation plus exact source-frequency reconciliation; existing researcher exclusions'))
pd.DataFrame(cross).to_csv(O/'physical_site_crosswalk.csv',index=False)
files=['outputs/biomass_resolution_20261102/review_annual_input.csv','outputs/biomass_support_20261102/animal_pool_inventory.csv','outputs/biomass_support_20261102/private/pool_members.csv','outputs/all_program_station_refits_20261102/offsets.csv','outputs/biomass_exclusion_reconciliation_20261102/annual_donor_corrections.csv','outputs/biomass_next_steps_20261102/kfm_year_holds.csv','data/MBON/SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv','data/MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv']
record=dict(recent_exact_matches=54,annual_corrections=99,verified_zero=1,cross_program_recent_records=14,sparse_allowed_samples=1,scope_unchanged=True,input_hashes={f:hashlib.sha256(Path(f).read_bytes()).hexdigest() for f in files})
(O/'scope_checks.json').write_text(json.dumps(record,indent=2));print('54 source matches; 99 annual corrections; 100 affected records; unchanged counts/years/other histories verified.')
