"""Independent checks of retained years, source counts, overlap and reconciliation."""
import csv, hashlib, itertools, math, os
from pathlib import Path
from collections import defaultdict, Counter
root=Path(__file__).resolve().parents[2]
dr=Path(os.environ.get('DONHAM_DATA_REPO',root.parent/'sbc-donham-kelp-mpa-data'))
def read(p):
    with p.open() as f:return list(csv.DictReader(f))
def tab(n):return read(root/f'tables/table_s_primary_size_{n}.csv')
for x in tab('inputs'):
    p=(root if x['repository']=='analysis' else dr)/x['path']
    assert hashlib.md5(p.read_bytes()).hexdigest()==x['md5'],p
inv=tab('inventory'); summaries=tab('comparisons')
keys=[x for x in read(root/'data/harmonized/reference_primary_keys.csv') if x['Resp']=='Bio' and x['Taxa'] in ('P. interruptus','S. purpuratus')]
expected={(x['series_id'],y,s) for x in keys for y in x['observed_years'].split(',') for s in ('mpa','reference')}
key=lambda x:(x['series_id'],x['year'],x['status'])
assert {key(x) for x in inv}==expected and len(inv)==len(expected)
private=root/'local_archive/primary_size_inventory_2026-10-30'
use=read(private/'record_membership.csv');bykey=defaultdict(list)
sf=read(dr/'data/ALL_sizefreq_2024.csv');dl=read(dr/'data/LTER/Lobster_Abundance_All_Years_20230831.csv')
for x in use:
    original=(sf if x['source_file']=='ALL_sizefreq_2024.csv' else dl)[int(x['source_row'])-1]
    assert float(x['count'])==float(original['count' if x['source_file']=='ALL_sizefreq_2024.csv' else 'COUNT'])
    bykey[key(x)].append(x)
for x in inv:
    if float(x['value'])>0:
        assert float(x['n_sized'])==sum(float(z['count']) for z in bykey[key(x)])
for s in summaries:
    rows=[x for x in inv if x['series_id']==s['series_id']]
    positive=[x for x in rows if float(x['value'])>0]
    assert int(s['positive_side_years'])==len(positive)
    assert int(s['under10_side_years'])==sum(float(x['n_sized'])<10 for x in positive)
    assert int(s['singleton_side_years'])==sum(float(x['n_sized'])==1 for x in positive)
# Calculate shared records independently by original file/row, never by matching lengths.
byrecord=defaultdict(list)
for x in use:byrecord[x['source_file'],x['source_row']].append(x)
sharing=Counter()
for xs in byrecord.values():
    for a,b in itertools.combinations(sorted(xs,key=lambda x:x['series_id']),2):
        if a['series_id']!=b['series_id']:
            sharing[a['series_id'],b['series_id'],a['year'],a['status'],b['status']]+=float(a['count'])
observed={(x['series_id_a'],x['series_id_b'],x['year'],x['status_a'],x['status_b']):float(x['shared_measured_animals']) for x in tab('sharing')}
assert dict(sharing)==observed
# Existing independently reconciled primary lobster means use the same full-year pools.
lookup={(x['MPA'],x['year'],x['status']):x for x in inv if x['Taxa']=='P. interruptus' and x['Source']=='PISCO' and float(x['value'])>0}
n=0
for x in read(root/'tables/table_s_lobster_reconciliation_annual.csv'):
    y=lookup.get((x['CA_MPA_Name_Short'],x['year'],x['status']))
    if y and x['primary_pool_mass'] not in ('','NA'):
        assert math.isclose(float(y['mean_mass']),float(x['primary_pool_mass']),rel_tol=1e-12)
        n+=1
assert n>0
for x in tab('direct_validation'):assert abs(float(x['reconstructed'])-float(x['value']))<1e-10
missing=read(private/'missing_lengths.csv');assert sum(float(x['COUNT']) for x in missing)==5
assert all(float(x['SIZE_MM'])==-99999 for x in missing)
print(f'PASS: {len(keys)} histories, {len(inv)} exact retained side-years; source counts, shared records, {n} prior lobster pool means, direct biomass and missing-length flags verified.')
