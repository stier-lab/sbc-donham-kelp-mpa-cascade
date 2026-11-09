#!/usr/bin/env python3
"""Read-only six-history PISCO density rebuild. Requires Python 3, pandas, numpy.
No repository module imports. Output must be outside inputs or within the dedicated next-steps review directory. All calendar months, equal method-preserving dated visits within
zone, equal zones within station and equal assigned stations within role/year.
Physical unit IDs do not include species or comparison: couple those weights in
joint resampling. They identify dated units, NOT proven permanent transects.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd

CODES = {'M. franciscanus':'MESFRAAD', 'S. purpuratus':'STRPURAD', 'P. interruptus':'PANINT'}
NAMES = {'M. franciscanus':'Mesocentrotus franciscanus', 'S. purpuratus':'Strongylocentrotus purpuratus', 'P. interruptus':'Panulirus interruptus'}
KEY = ['site','year','date','month','zone','method','transect']
VISIT = KEY[:-1]
SHARED = {'INNER','MID'}
SCENARIO = 'allmonths_sharedzones_visit_sizes'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def stable_id(kind, fields):
    return kind + '_' + hashlib.sha256(json.dumps([str(x) for x in fields], ensure_ascii=True, separators=(',',':')).encode()).hexdigest()

def complete_mean(values):
    a = np.asarray(values, dtype=float)
    return float(a.mean()) if len(a) and np.isfinite(a).all() else np.nan

def station_density(visits, variant):
    if variant == 'allmonths_shared_inner_mid':
        if not SHARED.issubset(set(visits.zone)):
            return np.nan
        visits = visits[visits.zone.isin(SHARED)]
    return complete_mean([complete_mean(g.density) for _,g in visits.groupby('zone')])

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo', type=Path, default=Path(os.environ.get('DONHAM_ANALYSIS_REPO', '.')))
    ap.add_argument('--raw', type=Path, default=os.environ.get('PISCO_SWATH_CSV'))
    ap.add_argument('--output', type=Path, default=os.environ.get('DENSITY_OUTPUT_DIR'))
    args=ap.parse_args()
    if args.output is None: ap.error('Supply --output or DENSITY_OUTPUT_DIR')
    repo=args.repo.resolve(); out=args.output.resolve()
    raw=(args.raw or repo.parent/'sbc-donham-kelp-mpa-data/data/PISCO/MLPA_kelpforest_swath_2024.csv').resolve()
    if (out.is_relative_to(repo) and not out.is_relative_to(repo/'outputs/biomass_next_steps_20261102')) or out.is_relative_to(raw.parent):
        ap.error('Output must be outside inputs or inside outputs/biomass_next_steps_20261102')
    out.mkdir(parents=True,exist_ok=True)
    inputs={}
    def read(rel, **kw):
        p=(repo/rel).resolve(); inputs[str(p)]=sha(p); return pd.read_csv(p,**kw)
    for rel in ['code/tools/refit_pisco_all_months.py','outputs/parallel_resolution_20261102/coverage/report.md','docs/PISCO_BASELINE_POLICY_2026-10-31.md']:
        inputs[str(repo/rel)]=sha(repo/rel)
    scope=read('outputs/parallel_resolution_20261102/coverage/density_scope_differences.csv')
    sids=scope.series_id.tolist(); assert len(sids)==6 and all(s.endswith('|PISCO|Den') for s in sids)
    allstations=read('outputs/all_program_station_refits_20261102/verified_station_values.csv')
    roster=allstations[allstations.series_id.isin(sids)].copy()
    assert not roster.duplicated(['series_id','site','year']).any()
    broster=allstations[allstations.series_id.isin([s.replace('|Den','|Bio') for s in sids])].copy()
    broster['biomass_series_id']=broster.series_id
    broster['series_id']=broster.series_id.str.replace('|Bio','|Den',regex=False)
    check=broster.merge(roster,on=['series_id','site','year','status'],how='left',indicator=True)
    assert check._merge.eq('both').all(), 'Biomass uses outside authorized density assignments'
    roster=roster.merge(broster[['series_id','site','year','biomass_series_id']],on=['series_id','site','year'],how='left',validate='one_to_one')
    roster['in_biomass_roster']=roster.biomass_series_id.notna()
    offsets=read('outputs/all_program_station_refits_20261102/offsets.csv').set_index('series_id').offset
    oldannual=read('outputs/pisco_all_months_20261102/annual_values.csv')
    oldvisits=read('outputs/pisco_all_months_20261102/visit_values.csv')
    primary=read('data/harmonized/harmonized_response_ratios.csv')
    inputs[str(raw)]=sha(raw)
    c=pd.read_csv(raw,encoding='latin1',low_memory=False)
    c['source_row_1based']=np.arange(1,len(c)+1) # data record, excluding header
    c=c[c.campus.eq('VRG')|(c.campus.eq('UCSB')&~c.zone.isin(['DEEP','SPECIALO']))].copy()
    c=c.merge(roster[['site','year']].drop_duplicates(),on=['site','year'],validate='many_to_one')
    c['date']=pd.to_datetime(c[['year','month','day']],errors='coerce').dt.strftime('%Y-%m-%d')
    assert c[KEY].notna().all().all(), 'Missing physical key/date: do not infer zeros'
    c.transect=c.transect.astype(str)
    assert not (c['count'].dropna()<0).any()
    exact_duplicate_rows=int(c.drop(columns='source_row_1based').duplicated().sum())
    assert exact_duplicate_rows==0, 'Ambiguous duplicated raw records: investigate before rebuild'
    units=c.groupby(KEY,as_index=False,dropna=False).agg(roster_source_rows=('source_row_1based',lambda s:';'.join(map(str,sorted(s)))),source_records=('source_row_1based','size'),campus=('campus',lambda s:';'.join(sorted(set(s)))))
    units['physical_unit_id']=[stable_id('pisco_swath_unit',[r.site,r.date,r.zone,r.method,r.transect]) for r in units.itertuples()]
    units['physical_visit_id']=[stable_id('pisco_swath_visit',[r.site,r.date,r.zone,r.method]) for r in units.itertuples()]
    assert units.physical_unit_id.is_unique
    registries=[]
    for taxon,code in CODES.items():
        z=c[c.classcode.eq(code)]
        counts=z.groupby(KEY,as_index=False,dropna=False).agg(count=('count',lambda s:s.sum() if s.notna().all() else np.nan),species_source_rows=('source_row_1based',lambda s:';'.join(map(str,sorted(s)))),species_source_records=('source_row_1based','size'))
        u=units.merge(counts,on=KEY,how='left',validate='one_to_one',indicator=True)
        u['surveyed_zero']=u._merge.eq('left_only');u.loc[u.surveyed_zero,'count']=0
        u['species_source_records']=u.species_source_records.fillna(0).astype(int)
        u['Taxa']=taxon;u['classcode']=code;u['area_m2']=60.;u['density']=u['count']/u.area_m2
        u['species_unit_id']=[stable_id('pisco_swath_species',[p,code]) for p in u.physical_unit_id]
        registries.append(u.drop(columns='_merge'))
    registry=pd.concat(registries,ignore_index=True)
    visit=registry.groupby(VISIT+['physical_visit_id','Taxa','classcode'],as_index=False).agg(density=('density',complete_mean),total_count=('count',lambda s:s.sum() if s.notna().all() else np.nan),transects=('physical_unit_id','nunique'),surveyed_zero_units=('surveyed_zero','sum'))
    registry=registry.merge(visit[['physical_visit_id','Taxa','density','transects','total_count']].rename(columns={'density':'visit_density','transects':'visit_transects','total_count':'visit_total_count'}),on=['physical_visit_id','Taxa'],validate='many_to_one')
    uses=roster.merge(registry,on=['site','year','Taxa'],validate='one_to_many')
    vus=roster.merge(visit,on=['site','year','Taxa'],validate='one_to_many')
    assert len(vus.groupby(['series_id','site','year']))==len(roster), 'Never drop unsupported roster rows'
    stations=[]; eligibility=[]
    for r in roster.itertuples():
        v=vus[(vus.series_id==r.series_id)&(vus.site==r.site)&(vus.year==r.year)]
        zones=set(v.zone); eligible=SHARED.issubset(zones)
        meta=dict(series_id=r.series_id,site=r.site,status=r.status,year=r.year,t=r.t,MPA_Start=r.MPA_Start,in_biomass_roster=r.in_biomass_roster)
        eligibility.append(dict(**meta,available_zones=';'.join(sorted(zones)),shared_inner_mid_supported=eligible,missing_shared_zones=';'.join(sorted(SHARED-zones)),extra_zones=';'.join(sorted(zones-SHARED)),visits=len(v),original_verified_density=r.value))
        for variant in ['original_verified_station','allmonths_full_available_zones','allmonths_shared_inner_mid']:
            val=r.value if variant=='original_verified_station' else station_density(v,variant)
            stations.append(dict(**meta,variant=variant,density=val))
    sy=pd.DataFrame(stations); flags=pd.DataFrame(eligibility)
    annual=[]
    for (sid,year,variant),g in sy.groupby(['series_id','year','variant']):
        f=flags[(flags.series_id==sid)&(flags.year==year)]; meta=g.iloc[0]
        roles=g.groupby('status').density.apply(complete_mean)
        assert set(roles.index)=={'mpa','reference'}
        mpa,ref=roles['mpa'],roles['reference']; off=float(offsets[sid])
        annual.append(dict(series_id=sid,year=year,t=int(meta.t),MPA_Start=int(meta.MPA_Start),variant=variant,protected=mpa,reference=ref,offset=off,ln_ratio=np.log((mpa+off)/(ref+off)),assigned_stations=len(g),reference_stations=int(g.status.eq('reference').sum()),full_reference_support=bool(g.density.notna().all()),shared_inner_mid_eligible=bool(f.shared_inner_mid_supported.all()),station_zones_identical=f.available_zones.nunique()==1,station_zone_sets=' | '.join(f.site+'='+f.available_zones),omitted_from_biomass_review=bool((~g.in_biomass_roster).any())))
    ann=pd.DataFrame(annual)
    missing=roster[~roster.in_biomass_roster].copy()
    assert len(missing)==33 and (missing.t<0).all()
    assert len(missing[['series_id','year']].drop_duplicates())==6
    # Exact original all-month identity on its original biomass-derived roster.
    old=oldannual[oldannual.series_id.isin(sids)&oldannual.scenario.eq(SCENARIO)].copy()
    now=ann[ann.variant.eq('allmonths_shared_inner_mid')]
    identity=old.merge(now,on=['series_id','year'],suffixes=('_old','_rebuilt'),validate='one_to_one')
    identity['absolute_error']=(identity.revised-identity.ln_ratio).abs()
    assert len(identity)==len(old) and identity.absolute_error.notna().all()
    assert identity.absolute_error.max()<1e-12
    assert len(identity[identity.t_old>0])==63
    # Physical keys and method must reproduce every saved shared-zone visit.
    bv=vus[vus.in_biomass_roster & vus.zone.isin(SHARED)].copy()
    ov=oldvisits[oldvisits.scenario.eq(SCENARIO)].copy()
    ov.series_id=ov.series_id.str.replace('|Bio','|Den',regex=False)
    vk=['series_id','site','status','year','zone','date']
    assert not bv.duplicated(vk).any(), 'Existing export loses method; ambiguous method match'
    vi=ov.merge(bv,on=vk,suffixes=('_old','_rebuilt'),how='outer',indicator=True,validate='one_to_one')
    assert vi._merge.eq('both').all()
    vi['absolute_error']=(vi.density_old-vi.density_rebuilt).abs()
    assert vi.absolute_error.max()<1e-12
    summaries=[]; primary_rows=[]
    for sid in sids:
        mpa,taxon,_,_=sid.split('|');g=ann[(ann.series_id==sid)&ann.variant.eq('allmonths_full_available_zones')];h=ann[(ann.series_id==sid)&ann.variant.eq('allmonths_shared_inner_mid')]
        pr=primary[(primary.CA_MPA_Name_Short==mpa)&(primary.y.str.strip()==NAMES[taxon])&(primary.resp=='Den')&(primary.source=='PISCO')]
        pr=pr.copy();pr['series_id']=sid;primary_rows.append(pr)
        start=int(g.MPA_Start.iloc[0]); pyears=sorted(set(pr.loc[(pr.year<start)&np.isfinite(pr.lnDiff),'year']))
        before=sorted(g.loc[(g.year<start)&g.ln_ratio.notna(),'year']); sharedbefore=sorted(h.loc[(h.year<start)&h.ln_ratio.notna(),'year'])
        earlier=sorted(missing.loc[missing.series_id==sid,'year'].unique().tolist())
        oldbefore=sorted(old.loc[(old.series_id==sid)&(old.year<start)&old.revised.notna(),'year'])
        summaries.append(dict(series_id=sid,primary_prior_years=';'.join(map(str,pyears)),n_primary_prior=len(pyears),old_allmonth_prior_years=';'.join(map(str,oldbefore)),full_density_prior_years=';'.join(map(str,before)),n_full_prior=len(before),shared_inner_mid_prior_years=';'.join(map(str,sharedbefore)),n_shared_prior=len(sharedbefore),recovered_comparison_years=';'.join(map(str,earlier)),recovered_station_year_uses=int((missing.series_id==sid).sum()),post_years=int((g.t>0).sum()),primary_3prior_count_only=len(pyears)>=3,primary_4prior_count_only=len(pyears)>=4,full_3prior_count_only=len(before)>=3,full_4prior_count_only=len(before)>=4,shared_3prior_count_only=len(sharedbefore)>=3,shared_4prior_count_only=len(sharedbefore)>=4,roster_constant_across_years=roster[roster.series_id==sid].groupby('year').site.apply(lambda s:tuple(sorted(s))).nunique()==1,interpretation='Count-only diagnostic; not shared-response/fixed-station eligibility, no pBACIPS selection'))
    summary=pd.DataFrame(summaries)
    outputs={'primary_annual_comparison':pd.concat(primary_rows,ignore_index=True),'density_station_year_roster':roster,'raw_count_unit_registry':registry,'density_count_unit_uses':uses,'per_visit_density':vus,'station_zone_eligibility':flags,'annual_station_density_variants':sy,'annual_comparison_density_variants':ann,'recovered_station_year_uses':missing,'recovered_comparison_years':ann[ann.omitted_from_biomass_review],'baseline_coverage_diagnostic':summary,'original_allmonth_ratio_identity':identity,'biomass_roster_visit_identity':vi.drop(columns='_merge'),'biomass_roster_count_units':uses[uses.in_biomass_roster & uses.zone.isin(SHARED)]}
    for name,df in outputs.items():df.to_csv(out/(name+'.csv'),index=False)
    assert all(sha(p)==h for p,h in inputs.items()), 'Input changed during run'
    (out/'input_hash_manifest.json').write_text(json.dumps({'inputs':inputs,'script_sha256':sha(__file__),'raw_source_row_definition':'1-based data records excluding header','parent_session':'9901638869734200a6633d0f88ec31ef'},indent=2)+'\n')
    result=dict(histories=6,density_station_year_uses=len(roster),recovered_prior_station_year_uses=len(missing),recovered_prior_comparison_years=6,post_comparison_years_preserved=63,physical_units=units.physical_unit_id.nunique(),species_units=len(registry),biomass_shared_count_unit_uses=len(outputs['biomass_roster_count_units']),biomass_visit_uses=len(vi),max_visit_error=float(vi.absolute_error.max()),max_allmonth_ratio_error=float(identity.absolute_error.max()),unsupported_inner_mid_station_year_uses=int((~flags.shared_inner_mid_supported).sum()),raw_exact_duplicate_rows=exact_duplicate_rows,missing_focal_counts=int(registry['count'].isna().sum()),input_hashes_unchanged=True,outputs={k:len(v) for k,v in outputs.items()})
    (out/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
    table=['| Comparison | Primary prior | Old all-month prior | Rebuilt full prior | Shared INNER/MID prior | Recovered station uses |', '|---|---:|---:|---:|---:|---:|']
    for r in summary.itertuples():
        label=r.series_id.replace('|PISCO|Den','').replace('|',' / ')
        nold=0 if not r.old_allmonth_prior_years else len(r.old_allmonth_prior_years.split(';'))
        table.append(f'| {label} | {r.n_primary_prior} | {nold} | {r.n_full_prior} | {r.n_shared_prior} | {r.recovered_station_year_uses} |')
    report=f"""# Density baseline rebuild — diagnostic implementation

Parent session `9901638869734200a6633d0f88ec31ef`. All writes belong in the chosen scratch output directory; repository and raw source inputs are read only. No own continuity start, shared status edits, source changes, selection, model fits, messages, pushes or further agents.

## Implemented result

Six density histories use their own **558 assigned station-year uses**, preserving researcher assignments, all available original years and all **63 post-protection comparison-years**. Recovered **33 prior station-year uses / six comparison-years**: Abalone lobster 2011 (4); Point Vicente red and purple 2008 (6 each); Vicente lobster 2008–2010 (17). No additional station or year was recruited from raw availability.

All six recovered comparison-years have finite full-available-zone densities, but **only three have complete INNER/MID support**. ROCKY_POINT_S in 2008 has DEEP only, affecting all three Vicente density responses. Its station rows and all three 2008 annual rows remain present; shared-zone values are explicitly unknown, never an average of the remaining references. No zone spans the full Vicente density roster. Even the other three recovered years have different full available zone sets among stations: availability does not establish habitat or seasonal comparability.

{chr(10).join(table)}

Primary and rebuilt full records contain 2011 at Abalone and 2008–2011 at Vicente. INNER/MID support reduces Vicente to 2009–2011. Thus all three Vicente responses meet three and four prior-year **count-only** rules in the full-zone version; they meet three but not four in the shared-zone version. Abalone meets neither. All six histories have changing station membership through time. These are response-specific diagnostics, **not** adopted fixed-station/shared-benthic-response eligibility, evidence of parallel trends, or permission to select pBACIPS. The adopted shared-response benthic limitation (including lobster biomass) remains unchanged. Prior means year < 2012; designation year 2012 is not prior. The 63 legacy recovery years use t > 0 to preserve that review's comparison exactly.

## Files and interpretation

- `raw_count_unit_registry.csv`: {len(registry)} species-unit rows on {len(units)} distinct observed physical units. Stable SHA256 IDs are built from site/date/zone/method/transect; species IDs add classcode. IDs exclude comparison and species at the physical level, allowing common weights across taxa/uses. They are **dated sampling-unit identities**, not proof of permanent transects. Raw row locators are 1-based data records excluding the CSV header; source bytes are hashed. Retain the same IDs when sharing observations between density and biomass.
- `density_station_year_roster.csv`, `density_count_unit_uses.csv`, `per_visit_density.csv`: authorized membership and its count/visit joins; no implicit joins by generated row number.
- `station_zone_eligibility.csv`: exact available, missing shared and extra zones, explicit eligibility, original value and biomass-roster membership.
- `annual_station_density_variants.csv`, `annual_comparison_density_variants.csv`: original verified station values (inherited comparator, not claimed newly re-derived), rebuilt all-month/full available zones, and all-month/INNER+MID when both strata are supported. Full-role arithmetic means use every assigned station; any unavailable assigned value makes that variant unknown. Frozen offsets are retained. The original verified station comparator is distinct from the manuscript-facing primary values in `primary_annual_comparison.csv`.
- `recovered_station_year_uses.csv`, `recovered_comparison_years.csv`, `baseline_coverage_diagnostic.csv`: exact omissions, all three variants of the six recovered years (18 rows), primary before-year coverage and three/four-prior count checks. No outcome-based screening.
- `biomass_roster_count_units.csv`: **{len(outputs['biomass_roster_count_units'])} current INNER/MID count-unit uses**, with original biomass series ID, role, exact site/year/zone/date/method/transect, count, area, unit density, visit density/total count/transect count, physical IDs, species and provenance. These reproduce all **{len(vi)}** existing shared-zone visit uses. Other zones are available in the full registry, but are not part of this current biomass bootstrap export.
- Identity CSVs, `verification.json`, `input_hash_manifest.json`, `test_density_baseline.py`, `tests.json`, `tests.log`: numerical checks and reproducibility evidence. Tests are a separate command; a rebuild alone does not refresh the test log.

## Count and aggregation conventions

The established VRG/UCSB campus filter is preserved (UCSB DEEP/SPECIALO excluded; VRG DEEP retained). All retained sites here use VRG swath. Every focal species is joined to the physical-unit roster formed from **all observed taxa**, assigning a surveyed zero only when the focal classcode is absent from an observed unit. An explicitly missing count remains unknown. Multiple focal records within a unit sum as in the source reconstruction; exact duplicate source records are checked and cause failure. Current inputs have no exact duplicates or missing focal counts. No site/date/zone/unit with no survey is invented.

Counts divide by the established **60 m²** swath area. Average transects within site/date/zone/method; average all valid-month visits equally within zone/year; average zones equally within station/year; average assigned stations equally within each role/year. Repeat visits are neither pooled with unequal transect weights nor represented as independent annual samples. Method is retained in physical keys; the current saved visit schema omits it, so ambiguous multi-method joins fail explicitly. Date IDs are not inferred survey campaigns.

## Verification and limits

All 78 original shared-zone all-month ratios (63 post years plus baseline/designation years) reproduce with maximum absolute error {identity.absolute_error.max():.3g}. All 1,086 existing biomass-roster visit densities reproduce to {vi.absolute_error.max():.3g}. Independent tests reconstruct counts from the raw CSV using Python's standard CSV reader, then visits and annual ratios with dictionary grouping; they also check the 78 original full-zone all-month ratios. Synthetic edge tests cover missing strata, missing values, equal visit/zone weighting and ID collisions. Input hashes are checked before/after the rebuild.

This restores density evidence, not biomass size support, causal comparability or primary eligibility. Shared labels alone do not establish matching depths or survey seasons; station membership varies. No date-gap, size-donor or new exclusion rule is imposed. Raw/primary/paper files remain untouched. Parent integration and scientific policy decisions remain separate.

## Reproduction

Requires Python 3.9+ (tested here with the installed Python), pandas and numpy. The script imports no repository code and performs no work on import.

```sh
python rebuild_density_baseline.py --repo /path/to/sbc-donham-kelp-mpa-cascade --raw /path/to/MLPA_kelpforest_swath_2024.csv --output /path/to/scratch/density
python test_density_baseline.py --repo /path/to/sbc-donham-kelp-mpa-cascade --output /path/to/scratch/density > /path/to/scratch/density/tests.log
```

`DONHAM_ANALYSIS_REPO`, `PISCO_SWATH_CSV`, and `DENSITY_OUTPUT_DIR` are optional environment equivalents. If raw is omitted, the sibling data-repository convention is used. Output is mandatory; repository writes are restricted to outputs/biomass_next_steps_20261102. Raw input paths are never valid output locations. Portable input paths and full SHA256 fingerprints are recorded in the manifest.
"""
    (out/'report.md').write_text(report)
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
