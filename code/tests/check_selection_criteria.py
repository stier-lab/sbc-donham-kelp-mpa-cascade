#!/usr/bin/env python3
"""Check source freshness and manuscript-critical selection claims."""
import csv
import hashlib
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[2]
data = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else root.parent/'sbc-donham-kelp-mpa-data'
def rows(p):
    return list(csv.DictReader(p.open()))
def key(r):
    return tuple(r[k] for k in ('MPA','Source','Taxa','Resp'))
for manifest in ('table_s_selection_source_inputs.csv','table_s_selection_criteria_inputs.csv'):
    for r in rows(root/'tables'/manifest):
        base = data if r.get('repository') == 'data' else root
        p = base/r['input']
        assert hashlib.md5(p.read_bytes()).hexdigest() == r['md5'], f'Stale selection audit: {p}'
coverage = rows(root/'tables/table_s_selection_history_coverage.csv')
routes = {key(r):r for r in rows(root/'tables/table_s_processing_routes.csv')}
assert len(coverage) == len(routes) == 143
assert {key(r) for r in coverage} == set(routes)
exceptions = [r for r in coverage if r['observation_within_one_year']=='0']
assert len(exceptions)==24
assert {r['MPA'] for r in exceptions} == {'Anacapa Island SMR 2003','South Point SMR'}
assert all(routes[key(r)]['file_label']=='CI' for r in exceptions), 'Near-date exceptions must not be described as before/after contrasts'
assert all(int(r['span_years'])>=10 and int(r['after_years'])>=4 for r in coverage)
branch = [r for r in rows(root/'tables/table_s_selection_branch_coverage.csv') if r['scope']=='all_retained_branch_responses']
assert all(int(r['pairs_three_through'])==0 for r in branch if r['branch']=='benthic')
assert {r['MPA'] for r in branch if r['branch']=='fish' and int(r['pairs_four_through'])>0} == {'Campus Point SMCA','Naples SMCA','Point Vicente SMCA'}
checked = [r for r in rows(root/'tables/table_s_selection_branch_coverage.csv') if r['scope']=='fish_count_length_checked']
assert {r['MPA'] for r in checked if int(r['pairs_four_through'])>0} == {'Campus Point SMCA','Naples SMCA','Point Vicente SMCA'}
overlap = rows(root/'tables/table_s_selection_judith_pixel_overlap.csv')
assert next(r for r in overlap if r['group_a']==r['group_b']=='MPA')['shared_pixels']=='364'
print('Selection audit inputs are current; 143 histories, 24 CI exceptions, branch coverage and pixel-overlap claims agree.')
