#!/usr/bin/env python3
"""Export retained processing routes without relabelling or refitting primary models."""
import csv
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
KEY=('MPA','Taxa','Source','Resp')

def key(row):return tuple(row[k]for k in KEY)

def build_routes(primary,coverage):
    indexed={key(r):r for r in coverage}
    if len(indexed)!=len(coverage) or len({key(r)for r in primary})!=len(primary):
        raise ValueError('Duplicate retained time-series identity')
    if {key(r)for r in primary}!=set(indexed):raise ValueError('Coverage cohort differs from primary')
    rows=[]
    for p in primary:
        c=indexed[key(p)]
        label=p['AnalysisType']
        if label!=c['primary_design']:raise ValueError('Conflicting processing label')
        before=int(c['n_before'])
        special=(p['MPA']=='Campus Point SMCA' and p['Taxa']=='P. interruptus' and p['Source']=='LTER' and p['Resp']=='Bio')
        route=('direct_CI_trend_with_pBACIPS_file_label' if special else
               'PISCO_direct_CI' if p['Source']=='PISCO' else 'retained_'+label)
        rows.append({**{k:p[k]for k in KEY},'file_label':label,'model':p['Model'],
            'paired_before_years':before,'paired_after_years':int(c['n_after']),
            'three_years_each_side':before>=3 and int(c['n_after'])>=3,
            'construction_route':route,
            'pretrend_failure_status':'not_established_by_processing_label'})
    return rows

def main():
    with (ROOT/'data/sumstats_final.csv').open()as f:primary=list(csv.DictReader(f))
    with (ROOT/'tables/table_s_design_matched_coverage.csv').open()as f:coverage=list(csv.DictReader(f))
    rows=build_routes(primary,coverage)
    assert len(rows)==145
    assert Counter(r['file_label']for r in rows)=={'CI':108,'pBACIPS':37}
    assert all(r['file_label']=='CI'for r in rows if r['Source']=='PISCO')
    assert sum(r['file_label']=='CI'and r['paired_before_years']>=3 for r in rows)==24
    assert sum(r['construction_route']=='direct_CI_trend_with_pBACIPS_file_label'for r in rows)==1
    path=ROOT/'tables/table_s_processing_routes.csv'
    with path.open('w',newline='')as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print('Verified 145 routes: 84 CI records have <3 before-years; 24 have >=3; one direct CI trend retains a pBACIPS file label.')

if __name__=='__main__':main()
