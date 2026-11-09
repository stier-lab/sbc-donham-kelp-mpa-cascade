#!/usr/bin/env python3
"""Audit actual KFM predator count-site membership; no primary outputs modified.

Uses the size-feasibility audit's independently reconstructed count-site means.
Presence includes observed zeros. A missing survey is never made into a zero.
"""
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
def read(p):
    with (ROOT/p).open() as f:
        return list(csv.DictReader(f))
def write(name, rows):
    with (ROOT/'tables'/('table_s_kfm_continuity_'+name+'.csv')).open('w') as f:
        w=csv.DictWriter(f,fieldnames=rows[0],lineterminator='\n'); w.writeheader();w.writerows(rows)

paths=['tables/table_s_kfm_size_site_year_matches.csv','data/harmonized/reference_primary_keys.csv',
       'data/harmonized/harmonized_response_ratios.csv','data/sumstats_final.csv']
rows=read(paths[0]);keys=read(paths[1]);primary=read(paths[3])
lookup={(r['MPA'],r['Taxa']):r for r in keys if r['Source']=='KFM' and r['Resp']=='Den'}
groups=defaultdict(list)
for r in rows:groups[(r['MPA'],r['species'],r['status'],r['site_id'])].append(r)
sites=[];persistent=set()
for (mpa,taxon,status,site),rr in sorted(groups.items()):
    k=lookup[(mpa,taxon)];before=set(k['before_years'].split(','));after=set(k['after_years'].split(','))
    years={r['year'] for r in rr};b=years&before;a=years&after
    if a and b:persistent.add((mpa,taxon,status,site))
    sites.append(dict(MPA=mpa,Taxa=taxon,AnalysisType=k['AnalysisType'],status=status,site_id=site,
        site_name=rr[0]['site_name'],n_before=len(b),n_after=len(a),
        before_years=';'.join(sorted(b)),after_years=';'.join(sorted(a)),
        membership='both_periods' if a and b else 'after_only' if a else 'before_only'))
write('sites',sites)
annual=[];dropped=[]
rr_lookup={(x['CA_MPA_Name_Short'],int(x['year'])):x for x in read(paths[2]) if x['source']=='KFM' and x['y']=='Semicossyphus pulcher' and x['resp']=='Den'}
for mpa,taxon in sorted({(s['MPA'],s['Taxa']) for s in sites if s['Taxa']=='S. pulcher' and s['AnalysisType']=='pBACIPS'}):
    rr=[r for r in rows if r['MPA']==mpa and r['species']==taxon]
    by=defaultdict(list)
    for r in rr:by[(r['year'],r['status'])].append(r)
    means={k:sum(float(r['count_density']) for r in a)/len(a) for k,a in by.items()}
    offset=min(v for v in means.values() if v>0)/2
    for year in sorted({r['year'] for r in rr}):
        a={'MPA':mpa,'year':int(year),'offset':offset}
        if any(not any((mpa,taxon,status,r['site_id']) in persistent for r in by[(year,status)]) for status in ['mpa','reference']):
            dropped.append({'MPA':mpa,'year':int(year),'reason':'No persistent reference or protected site observed; omitted from both comparison arms'})
            continue
        for status in ['mpa','reference']:
            group=by[(year,status)];core=[r for r in group if (mpa,taxon,status,r['site_id']) in persistent]
            assert core,(mpa,year,status,'No persistent site: cannot preserve annual cohort')
            a[status+'_original']=means[(year,status)]
            a[status+'_persistent']=sum(float(r['count_density']) for r in core)/len(core)
            a[status+'_n_original']=len(group);a[status+'_n_persistent']=len(core)
        for scenario in ['original','persistent']:
            a['lnRR_'+scenario]=math.log((a['mpa_'+scenario]+offset)/(a['reference_'+scenario]+offset))
        match=rr_lookup[(mpa,int(year))]
        assert abs(a['lnRR_original']-float(match['lnDiff']))<1e-10,(mpa,year,'Original lnRR mismatch')
        annual.append(a)
write('annual',annual)
manifest={'scope':'KFM primary lobster and sheephead density site membership; fixed-site refit inputs only for two pBACIPS sheephead histories',
 'inputs_sha256':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
 'site_histories':len(sites),'after_only_in_pbacips':[s for s in sites if s['AnalysisType']=='pBACIPS' and s['membership']=='after_only'],
 'annual_lnRR_original_verified':len(annual),'years_dropped_from_both_arms':dropped,
 'method':'Retain sites with observations in both inherited before/after periods; compare original versus persistent membership over the same available paired years, retaining zeros, original offset, selected model and clock. No threshold or clock change.',
 'limits':'Does not audit non-predator biomass contributions, habitat equivalence, or every source program. Persistent membership does not require sampling in every year.'}
(ROOT/'tables/kfm_continuity_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps({k:v for k,v in manifest.items() if k not in ['inputs_sha256','after_only_in_pbacips']},indent=2))
print('After-only sites in pBACIPS:',[(s['MPA'],s['site_name']) for s in manifest['after_only_in_pbacips']])
