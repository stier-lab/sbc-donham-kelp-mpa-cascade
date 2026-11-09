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
        if special and label!='CI':raise ValueError('Campus Point lobster trend must be classified CI')
        route=('direct_CI_trend' if special else
               'PISCO_fixed_station_pBACIPS' if p['Source']=='PISCO' and label=='pBACIPS' else 'PISCO_direct_CI' if p['Source']=='PISCO' else 'retained_'+label)
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
    assert len(rows)==143
    assert Counter(r['file_label']for r in rows)=={'CI':103,'pBACIPS':40}
    assert sum(r['file_label']=='pBACIPS' for r in rows if r['Source']=='PISCO')==6
    assert sum(r['file_label']=='CI'and r['paired_before_years']>=3 for r in rows)==18
    assert sum(r['construction_route']=='direct_CI_trend'for r in rows)==1
    path=ROOT/'tables/table_s_processing_routes.csv'
    with path.open('w',newline='')as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");writer.writeheader();writer.writerows(rows)
    print('Verified 143 routes: 85 CI records have <3 before-years; 18 have >=3; one direct CI trend is consistently labeled CI.')

if __name__=='__main__':main()
