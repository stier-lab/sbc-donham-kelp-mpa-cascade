#!/usr/bin/env python3
"""Audit observed coverage and delivered pixel overlap; never alter primary data."""
import csv
import hashlib
from pathlib import Path
from openpyxl import load_workbook

ROOT = Path(__file__).resolve().parents[2]

def read(path):
    return list(csv.DictReader(path.open()))

def key(row):
    return tuple(row[x] for x in ('MPA', 'Source', 'Taxa', 'Resp'))

def write(name, rows):
    with (ROOT / 'tables' / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        w.writeheader()
        w.writerows(rows)

files = [ROOT / 'data/sumstats_final.csv',
         ROOT / 'data/harmonized/reference_primary_keys.csv',
         ROOT / 'data/external/tom_bell_2026-10-20/03_MPA_Ref_Locations.xlsx']
hashes = {p: hashlib.md5(p.read_bytes()).hexdigest() for p in files}
primary = {key(r): r for r in read(files[0])}
records = read(files[1])
assert len(primary) == len(records) == 143
assert {key(r) for r in records} == set(primary)
rows = []
for r in records:
    year = sorted({int(x) for x in r['observed_years'].split(',')})
    start = int(float(primary[key(r)]['MPA_Start']))
    rows.append({**{x: r[x] for x in ('MPA', 'Source', 'Taxa', 'Resp')},
                 'designation_year': start, 'first_year': min(year), 'last_year': max(year),
                 'span_years': max(year)-min(year), 'paired_survey_years': len(year),
                 'strict_before_years': sum(y < start for y in year),
                 'baseline_years_through_designation': sum(y <= start for y in year),
                 'after_years': sum(y > start for y in year),
                 'observation_within_one_year': int(any(abs(y-start) <= 1 for y in year))})
assert all(r['span_years'] >= 10 and r['after_years'] >= 4 for r in rows)
write('table_s_selection_history_coverage.csv', rows)

sheet = list(load_workbook(files[2], read_only=True, data_only=True)['Pixel Locations'].values)
sets = {}
for j in range(1, len(sheet[0]), 2):
    if sheet[0][j] is None:
        continue
    sets[int(sheet[0][j]), sheet[1][j]] = {
        (round(float(r[j]), 6), round(float(r[j+1]), 6)) for r in sheet[2:]
        if j+1 < len(r) and isinstance(r[j], (int, float)) and isinstance(r[j+1], (int, float))}
overlap = []
for a in ('MPA', 'Ref_1', 'Ref_2'):
    for b in ('MPA', 'Ref_1', 'Ref_2'):
        x, y = sets[1, a], sets[17, b]
        overlap.append({'site_a': 1, 'label_a': 'San Miguel Island SC', 'group_a': a,
                        'site_b': 17, 'label_b': 'Judith Rk SMR', 'group_b': b,
                        'pixels_a': len(x), 'pixels_b': len(y), 'shared_pixels': len(x & y)})
write('table_s_selection_judith_pixel_overlap.csv', overlap)
assert hashes == {p: hashlib.md5(p.read_bytes()).hexdigest() for p in files}
write('table_s_selection_criteria_inputs.csv', [
    {'input': str(p.relative_to(ROOT)), 'md5': h} for p, h in hashes.items()])
exceptions = [r for r in rows if not r['observation_within_one_year']]
print(f'{len(rows)} histories meet span/post-year criteria; {len(exceptions)} lack a survey within one year.')
print('Exceptions:', sorted({r['MPA'] for r in exceptions}))
print('Judith/San Miguel shared MPA pixels:', overlap[0]['shared_pixels'])
