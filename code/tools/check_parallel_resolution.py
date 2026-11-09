"""Coordinator checks of load-bearing parallel-audit claims; no scientific writes."""
from pathlib import Path
import pandas as pd,json,hashlib
A=Path('outputs/parallel_resolution_20261102')
B=Path('outputs/pisco_biomass_borrowing_20261102')
v=pd.read_csv(B/'visit_biomass_support.csv');u=v[v.method.eq('unsupported')];annual=pd.read_csv(B/'biomass_annual.csv')
x=u.merge(annual[['series_id','year','t']],on=['series_id','year'],validate='many_to_one');post=x[x.t>0]
assert len(u)==38 and len(u[['site','date']].drop_duplicates())==24
assert len(post)==26 and len(post[['series_id','year']].drop_duplicates())==19
assert set(map(tuple,post[['series_id','year']].drop_duplicates().values))==set(map(tuple,annual[(annual.t>0)&~annual.supported][['series_id','year']].values))
s=pd.read_csv('outputs/pisco_timing_donors_20261102/size_sample_inventory.csv');z=v[v.method.eq('measured_same_visit')].merge(s,on=['series_id','site','status','year','zone','date'],validate='one_to_one');assert len(z)==624 and z.n_sized.eq(1).sum()==84
val=pd.read_csv(B/'heldout_interval_validation.csv');t=val[val.available&val.kind.eq('other_visit')];assert len(t)==8 and t.site.nunique()==2 and t.year.nunique()==1 and not t.Taxa.eq('S. purpuratus').any()
st=pd.read_csv('outputs/all_program_station_refits_20261102/verified_station_values.csv');pisco=[]
for sid,g in st[st.Source.eq('PISCO')].groupby('series_id'):
    p=g[g.t>0];ys=p.year.unique();r=p[(p.status=='reference')&p.value.notna()];sets=[set(r.loc[r.year==y,'site']) for y in ys];fixed=set.intersection(*sets) if sets else set()
    pisco.append((len(ys),bool(fixed)))
assert len(pisco)==105 and sum(not f for n,f in pisco)==17
assert sum(n for n,f in pisco if not f)==187
assert sum(n>0 and not f for n,f in pisco)==15
roster=pd.read_csv(A/'coverage/density_scope_differences.csv').set_index('series_id')
for sid in pd.read_csv(B/'summary.csv').series_id:
    den=sid.replace('|Bio','|Den');d=st[st.series_id.eq(den)];b=st[st.series_id.eq(sid)];lost=d.merge(b[['site','year']],on=['site','year'],how='left',indicator=True);lost=lost[lost._merge.eq('left_only')]
    assert len(lost)==roster.loc[den,'omitted_density_station_years']
    assert (lost.t>0).sum()==0
assert roster.omitted_prior_station_years.sum()==33
h=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text());assert all(Path(p).exists() and hashlib.sha256(Path(p).read_bytes()).hexdigest()==sha for p,sha in h.items())
print('Coordinator checks passed: unsupported use/year identity, singletons, temporal validation replication, fixed-reference coverage, baseline roster losses; 487 protected files unchanged.')
