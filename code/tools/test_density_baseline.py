#!/usr/bin/env python3
"""Independent csv/stdlib source reconstruction plus numerical edge tests.
Run from any directory: python test_density_baseline.py --output DIR --repo REPO.
Writes only DIR/tests.json; uses no repository Python imports.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
import pandas as pd


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--repo',type=Path,default=Path.cwd());a=p.parse_args();o=a.output
    spec=importlib.util.spec_from_file_location('density_rebuild',Path(__file__).with_name('rebuild_density_baseline.py'))
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    checks=[]
    def ok(name, condition, **detail):
        assert condition, name
        checks.append(dict(test=name,passed=True,**detail))
    def rows(path):
        with open(path,newline='',encoding='utf-8') as f:return list(csv.DictReader(f))
    def close(x,y):return math.isclose(float(x),float(y),rel_tol=1e-12,abs_tol=1e-12)
    reg=rows(o/'raw_count_unit_registry.csv'); roster=rows(o/'density_station_year_roster.csv'); uses=rows(o/'biomass_roster_count_units.csv')
    manifest=json.loads((o/'input_hash_manifest.json').read_text())
    ok('all_inputs_still_match_hash_manifest',all(hashlib.sha256(Path(k).read_bytes()).hexdigest()==v for k,v in manifest['inputs'].items()),inputs=len(manifest['inputs']))
    raw=next(k for k in manifest['inputs'] if k.endswith('MLPA_kelpforest_swath_2024.csv'))
    allowed={(r['site'],r['year']) for r in roster}; unit_roster=set(); counts=defaultdict(float); seen=set()
    def key(r):return (r['site'],f"{int(r['year']):04}-{int(r['month']):02}-{int(r['day']):02}",r['zone'],r['method'],r['transect'])
    with open(raw,encoding='latin1',newline='') as f:
        for r in csv.DictReader(f):
            if (r['site'],r['year']) not in allowed:continue
            if not (r['campus']=='VRG' or (r['campus']=='UCSB' and r['zone'] not in ('DEEP','SPECIALO'))):continue
            k=key(r);unit_roster.add(k)
            if r['classcode'] in ('MESFRAAD','STRPURAD','PANINT'):
                counts[k+(r['classcode'],)]+=float(r['count']);seen.add(k+(r['classcode'],))
    for r in reg:
        k=(r['site'],r['date'],r['zone'],r['method'],r['transect']); ck=k+(r['classcode'],)
        assert k in unit_roster
        assert close(r['count'],counts[ck]) and close(r['density'],counts[ck]/60)
        assert (r['surveyed_zero']=='True') == (ck not in seen)
    ok('independent_raw_count_and_surveyed_zero_reconstruction',len(reg)==3*len(unit_roster),species_units=len(reg),physical_units=len(unit_roster))
    ids=defaultdict(set)
    for r in reg:ids[(r['site'],r['date'],r['zone'],r['method'],r['transect'])].add(r['physical_unit_id'])
    ok('physical_ids_shared_across_species_not_species_replicates',all(len(s)==1 for s in ids.values()))
    ok('stable_ids_independent_of_input_order',mod.stable_id('unit',['a','b'])==mod.stable_id('unit',list(reversed(['b','a']))))
    ok('delimiter_collision_resistance',mod.stable_id('unit',['a|b','c'])!=mod.stable_id('unit',['a','b|c']))
    v=defaultdict(list)
    for r in uses:
        k=tuple(r[x] for x in ['biomass_series_id','site','status','year','zone','date']);v[k].append(float(r['density']))
    old=[r for r in rows(a.repo/'outputs/pisco_all_months_20261102/visit_values.csv') if r['scenario']==mod.SCENARIO]
    for r in old:
        k=tuple(r[x] for x in ['series_id','site','status','year','zone','date'])
        assert close(statistics.mean(v[k]),r['density'])
    ok('independent_all_biomass_visits_equal_mean_of_count_units',len(v)==len(old),visits=len(old),count_unit_uses=len(uses))
    # Independently apply the hierarchy to every annual variant made from raw counts.
    groups=defaultdict(list)
    for r in rows(o/'density_count_unit_uses.csv'):
        k=tuple(r[x] for x in ['series_id','year','status','site','zone','date','method']);groups[k].append(float(r['density']))
    zones=defaultdict(list)
    for k,vals in groups.items():zones[k[:5]].append(statistics.mean(vals))
    sites=defaultdict(dict)
    for k,vals in zones.items():sites[k[:4]][k[4]]=statistics.mean(vals)
    for variant in ('allmonths_full_available_zones','allmonths_shared_inner_mid'):
        roles=defaultdict(list)
        for k,z in sites.items():
            val=statistics.mean(z.values()) if variant.endswith('available_zones') else statistics.mean([z['INNER'],z['MID']]) if {'INNER','MID'}<=z.keys() else float('nan')
            roles[k[:3]].append(val)
        n=0
        for r in rows(o/'annual_comparison_density_variants.csv'):
            if r['variant']!=variant:continue
            k=(r['series_id'],r['year']); m=statistics.mean(roles[k+('mpa',)]); ref=statistics.mean(roles[k+('reference',)])
            val=math.log((m+float(r['offset']))/(ref+float(r['offset'])))
            if math.isnan(val):assert r['ln_ratio']==''
            else:assert close(val,r['ln_ratio'])
            n+=1
        ok('independent_annual_hierarchy_'+variant,n==84,comparison_years=n)
    lost=rows(o/'recovered_station_year_uses.csv')
    ok('exact_six_prior_years_33_station_uses_restored',len(lost)==33 and len({(r['series_id'],r['year']) for r in lost})==6 and all(float(r['t'])<0 for r in lost))
    identity=rows(o/'original_allmonth_ratio_identity.csv')
    post=[r for r in identity if float(r['t_old'])>0]
    ok('all_63_post_shared_ratios_unchanged',len(post)==63 and all(float(r['absolute_error'])<1e-12 for r in post))
    # Check original all-available-zone scenario independently too, on original roster.
    anns={(r['series_id'],r['year']):r for r in rows(o/'annual_comparison_density_variants.csv') if r['variant']=='allmonths_full_available_zones'}
    oldfull=[r for r in rows(a.repo/'outputs/pisco_all_months_20261102/annual_values.csv') if r['scenario']=='allmonths_oldpool' and r['series_id'] in {s['series_id'] for s in roster}]
    ok('all_78_original_full_zone_allmonth_ratios_unchanged',len(oldfull)==78 and all(close(r['revised'],anns[(r['series_id'],r['year'])]['ln_ratio']) for r in oldfull))
    flags=rows(o/'station_zone_eligibility.csv');bad=[r for r in flags if r['shared_inner_mid_supported']=='False']
    ok('no_silent_drop_of_2008_DEEP_only_reference',len(bad)==3 and all(r['site']=='ROCKY_POINT_S' and r['year']=='2008' and r['available_zones']=='DEEP' for r in bad))
    ok('strict_missing_mean',math.isnan(mod.complete_mean([1,float('nan')])) and math.isnan(mod.complete_mean([])))
    vv=pd.DataFrame({'zone':['INNER','INNER','MID'],'density':[0.,4.,10.]})
    ok('equal_visit_then_equal_zone_not_pooled_visits',mod.station_density(vv,'allmonths_shared_inner_mid')==6.)
    ok('absent_shared_zone_remains_unknown',math.isnan(mod.station_density(vv[vv.zone=='INNER'],'allmonths_shared_inner_mid')))
    ok('no_claim_of_permanent_transect_ids',mod.stable_id('unit',['site','2010-01-01','1'])!=mod.stable_id('unit',['site','2011-01-01','1']))
    result=dict(passed=len(checks),failed=0,checks=checks)
    (o/'tests.json').write_text(json.dumps(result,indent=2)+'\n'); print(json.dumps(result,indent=2))

if __name__=='__main__':main()
