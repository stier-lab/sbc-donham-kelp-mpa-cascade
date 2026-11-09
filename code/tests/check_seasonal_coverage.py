#!/usr/bin/env python3
"""Verify the seasonal coverage export against retained years and source hashes."""
import csv
import hashlib
import os
from pathlib import Path
root=Path(__file__).resolve().parents[2]
data=Path(os.environ.get('DONHAM_DATA_REPO',root.parent/'sbc-donham-kelp-mpa-data'))
def rows(name):return list(csv.DictReader((root/name).open()))
def key(x):return tuple(x[k] for k in ('MPA','Taxa','Source','Resp'))
manifest=rows('tables/table_s_seasonal_inputs.csv')
for x in manifest:
    p=(data if x['repository']=='data' else root)/x['path']
    assert hashlib.md5(p.read_bytes()).hexdigest()==x['md5'],f'Stale seasonal input: {p}'
expected={(key(x),y) for x in rows('data/harmonized/reference_primary_keys.csv') if x['Source']!='Landsat' for y in x['observed_years'].split(',')}
pairs=rows('tables/table_s_seasonal_paired_years.csv')
assert len(pairs)==len(expected)==len({(key(x),x['year']) for x in pairs})
assert {(key(x),x['year']) for x in pairs}==expected
for x in pairs:
    a=set(x['protected_months'].split(';'))-{''};b=set(x['reference_months'].split(';'))-{''}
    assert a and b, 'Dates missing for retained count coverage'
    assert int(x['n_common_months'])==len(a&b)
    assert (x['same_month_set']=='TRUE')==(a==b)
    assert all(1<=int(m)<=12 for m in a|b)
for s in rows('tables/table_s_seasonal_summary.csv'):
    r=[x for x in pairs if all(x[k]==s[k]for k in ['Source','Taxa','Resp'])]
    assert int(s['paired_years'])==len(r)
    assert int(s['no_shared_month'])==sum(int(x['n_common_months'])==0 for x in r)
print(f'Seasonal coverage verified for {len(expected)} retained response-years; all source hashes current. No seasonal refit is implied.')
