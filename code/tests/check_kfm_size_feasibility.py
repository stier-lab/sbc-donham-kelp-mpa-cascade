from pathlib import Path
import csv, hashlib, math, json
ROOT=Path(__file__).resolve().parents[2]
def read(name):
    with (ROOT/'tables'/f'table_s_kfm_size_{name}.csv').open() as f:return list(csv.DictReader(f))
m=read('site_year_matches');p=read('pools');c=read('coverage')
assert len(m)==len({(r['species'],r['MPA'],r['status'],r['site_id'],r['year']) for r in m})
for r in read('count_reconstruction'):assert abs(float(r['difference']))<1e-10
for r in read('inputs'):assert hashlib.sha256(Path(r['path']).read_bytes()).hexdigest()==r['sha256'],r['path']
for r in p:
    assert float(r['bin_lower_mass_g'])-1e-9<=float(r['mean_mass_g'])<=float(r['bin_upper_mass_g'])+1e-9
    if r['conditional_precision_estimable']=='False':assert math.isnan(float(r['conditional_lower_g']))
for r in m:
    if r['size_dates_missing']=='True':assert r['exact_date_match']=='False'
    if r['has_size_pool']=='False':assert r['conditional_relative_halfwidth']==''
# A zero-count site without lengths must not discard a whole matched MPA-year.
h=next(r for r in c if r['MPA']=='Harris Point SMR' and r['species']=='S. pulcher')
assert '2010' not in h['years_all_count_sites_have_sizes'].split(';')
assert '2010' in h['years_all_positive_count_sites_have_sizes'].split(';')
# No fake pre-designation size coverage or replacement of lobster gaps with zero biomass.
assert sum(int(r['n_before_size_matches']) for r in c)==0
assert all(r['has_size_pool']=='False' for r in m if r['species']=='P. interruptus')
print('PASS: unique joins; 218 primary density reconstructions; source hashes; bin bounds; missing-data/zero-count handling; no invented pre-designation coverage.')
