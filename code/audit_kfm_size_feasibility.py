"""Audit KFM predator sizes against retained density histories; never promote biomass.
Requires numpy/openpyxl. Public inputs are cached in the sibling data repository.
Counts are read from the exact MBON release used by the primary pipeline.
"""
from pathlib import Path
from collections import defaultdict, Counter
from datetime import datetime
import csv, json, hashlib, re, math, os
import numpy as np
import openpyxl
ROOT=Path(__file__).resolve().parents[1]
CACHE=Path(os.environ.get('KFM_SIZE_CACHE', str(ROOT.parent/'sbc-donham-kelp-mpa-data/data/source_archive/kfm_size_feasibility_2026-10-29')))
OUT=ROOT/'tables'; INPUTS=[]
def read(path,encoding='utf-8-sig'):
    path=Path(path);INPUTS.append(path)
    with path.open(encoding=encoding) as f:return list(csv.DictReader(f))
def write(name,rows):
    assert rows,name
    with (OUT/('table_s_kfm_size_'+name+'.csv')).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def number(v):
    try:return float(v)
    except (ValueError,TypeError):return math.nan
primary=read(ROOT/'data/sumstats_final.csv')
retained={(r['MPA'],r['Taxa']) for r in primary if r['Source']=='KFM' and r['Resp']=='Den'}
annual=read(ROOT/'data/harmonized/harmonized_response_ratios.csv')
short={'Bodianus pulcher':'S. pulcher','Semicossyphus pulcher':'S. pulcher','Panulirus interruptus':'P. interruptus'}
years=defaultdict(set)
for r in annual:
    if r['source']=='KFM' and r['resp']=='Den' and r['y'] in short:years[(r['CA_MPA_Name_Short'],short[r['y']])].add(int(r['year']))
mapping=read(ROOT/'data/MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv')
assert len({(r['site_id'],r['CA_MPA_Name_Short']) for r in mapping})==len(mapping)
by_site=defaultdict(list)
for r in mapping:
    if r['data_source']=='kfm':by_site[r['site_id']].append(r)
utils=(ROOT/'code/R/01_utils.R').read_text()
excluded=set(re.findall(r'"(a-k-\d+)"',re.search(r'EXCLUDED_KFM_SITES\s*<-\s*c\((.*?)\)',utils,re.S)[1]))
# Read source rows, not another program's annual averages. RDFC never supplies density here.
counts=defaultdict(list);source_summary=[]
for species,filename,method in [('S. pulcher','SBCMBON_kelp_forest_integrated_fish_20231022.csv','visualfish'),('P. interruptus','SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv','band')]:
    p=ROOT/'data/MBON'/filename;INPUTS.append(p);seen=Counter()
    with p.open() as f:
        for r in csv.DictReader(f):
            if r['data_source']!='kfm' or short.get(r['taxon_name'])!=species:continue
            seen[(r['sample_method'],r['proj_taxon_id'])]+=1
            if r['sample_method']==method:counts[(species,r['site_id'],int(r['date'][:4]))].append(r)
    source_summary.extend(dict(species=species,method=k[0],project_taxon_id=k[1],rows=n,used_for_count_matching=k[0]==method) for k,n in seen.items())
write('count_sources',source_summary)
# Historical source uses stable numeric NPS Site identifiers, not an approximate name join.
sizes=defaultdict(list);size_dates=defaultdict(set);size_source={};duplicates=[]
legacy=read(CACHE/'legacy/Fish size frequency data.csv',encoding='cp1252')
selected=[r for r in legacy if 'pulcher' in r['ScientificName']]
for r in selected:
    key=('S. pulcher',f"a-k-{int(r['Site']):02}",int(r['Year']))
    n=number(r['Count']);length=number(r['TotalLength'])
    assert n>0 and n==int(n) and length>0
    date=datetime.strptime(r['Date'],'%d%b%Y').date().isoformat();assert int(date[:4])==key[2]
    sizes[key].append(dict(n=n,lo=length,mid=length,hi=length,stage=r['ScientificName'].split(',')[-1].strip(),grouped=n>1,range_documented=math.isfinite(number(r['MinLength'])) and math.isfinite(number(r['MaxLength']))))
    size_dates[key].add(date);size_source[key]='Kushner2013_individual_or_group_records'
# Identical rows can be multiple same-sized fish: count and retain, never silently deduplicate.
exact_duplicates=len(selected)-len(set(tuple(r.values()) for r in selected))
# 2015 is a separate, nonoverlapping year, available as binned summaries only.
def workbook(path):
    INPUTS.append(path);s=openpyxl.load_workbook(path,read_only=True,data_only=True).active
    it=iter(s.values);header=next(it);return [dict(zip(header,r)) for r in it]
nh=workbook(CACHE/'2015/qryAppendixI_NHSF_2015.xlsx')
site_names=defaultdict(set)
for r in nh:site_names[r['Site_Name']].add(int(r['Site_Number']))
assert all(len(x)==1 for x in site_names.values())
fs=workbook(CACHE/'2015/qryAppendixH_FSF_2015.xlsx')
groups=defaultdict(list)
for r in fs:
    if 'pulcher' in r['ScientificName']:groups[(r['Site_Name'],r['ScientificName'])].append(r)
for (site,taxon),rows in groups.items():
    assert len(set(r['Cases_N'] for r in rows))==1
    n=rows[0]['Cases_N'];assert n>0
    assert abs(sum(r['Percent'] for r in rows)-100)<1e-5
    key=('S. pulcher',f"a-k-{next(iter(site_names[site])):02}",2015)
    assert key not in size_source or size_source[key]=='NPS2015_binned_summary'
    for r in rows:
        count=r['Percent']*n/100;assert abs(count-round(count))<1e-6
        if count<.5:continue
        nums=[float(x) for x in re.findall(r'\d+',str(r['Bin_Size']))]
        if '<' in str(r['Bin_Size']):lo=r['Minimum_size'];hi=nums[0]-1
        elif len(nums)==2:lo,hi=nums
        else:lo=nums[0];hi=r['Maximum_size']
        lo=max(lo,r['Minimum_size']);hi=min(hi,r['Maximum_size']);assert 0<lo<=hi
        sizes[key].append(dict(n=round(count),lo=lo,mid=(lo+hi)/2,hi=hi,stage=taxon.split(',')[-1].strip(),grouped=False,range_documented=True))
    size_source[key]='NPS2015_binned_summary'
# Inventory lobster size absence in these actual releases, not an assertion about all KFM collections.
inv=read(CACHE/'legacy/Invertebrate Size Frequency.csv',encoding='cp1252')
assert not any(number(r['Species'])==8001 for r in inv)
assert not any('interruptus' in r['Scientific_Name'] for r in nh)
# Conditional precision diagnostic only; no density*size biomass is promoted.
rng=np.random.default_rng(20261029);poolrows=[];poollookup={}
for key,rows in sorted(sizes.items()):
    n=np.array([r['n'] for r in rows]);length=np.array([r['mid'] for r in rows]);mass=.0144*length**3.04
    total=n.sum();mean=float(np.sum(n*mass)/total)
    if size_source[key].startswith('Kushner'):
        # Resample source rows, retaining each grouped count. Not independent-dive inference.
        draws=rng.multinomial(len(rows),np.repeat(1/len(rows),len(rows)),size=2000)
        boot=(draws@(n*mass))/(draws@n)
    else:
        draws=rng.multinomial(int(total),n/total,size=2000);boot=draws@mass/total
    lo,hi=np.quantile(boot,[.025,.975])
    degenerate = len(set(length)) < 2 or total < 2
    if degenerate:lo=hi=math.nan
    po=dict(species=key[0],site_id=key[1],year=key[2],source=size_source[key],n_size_records=len(rows),n_sized=int(total),n_juvenile=int(sum(r['n'] for r in rows if r['stage']=='juvenile')),n_grouped_records=sum(r['grouped'] for r in rows),n_grouped_records_without_range=sum(r['grouped'] and not r['range_documented'] for r in rows),n_size_dates=len(size_dates[key]),size_dates=';'.join(sorted(size_dates[key])),mean_mass_g=mean,conditional_precision_estimable=not degenerate,conditional_lower_g=float(lo),conditional_upper_g=float(hi),conditional_relative_halfwidth=float((hi-lo)/2/mean),max_record_mass_share=float(max(n*mass)/sum(n*mass)),bin_lower_mass_g=float(sum(r['n']*.0144*r['lo']**3.04 for r in rows)/total),bin_upper_mass_g=float(sum(r['n']*.0144*r['hi']**3.04 for r in rows)/total),method_compatible='unverified_roving_sizes_vs_visual_transects',depth_compatible='not_established',precision_interpretation='conditional_on_recorded_sizes_only_not_full_biomass_uncertainty')
    poolrows.append(po);poollookup[key]=po
write('pools',poolrows)
matched=[]
for (species,site,year),rs in sorted(counts.items()):
    for m in by_site[site]:
        mpa=m['CA_MPA_Name_Short']
        if (mpa,species) not in retained or year not in years[(mpa,species)]:continue
        if species=='P. interruptus' and site in excluded:continue
        finite=[r for r in rs if math.isfinite(number(r['count'])) and number(r['area'])>0]
        if not finite:continue
        key=(species,site,year);po=poollookup.get(key);dates={r['date'] for r in finite};sd=size_dates[key]
        transects=defaultdict(list)
        for r in finite:transects[r['transect_id']].append(number(r['count'])/number(r['area']))
        density=float(np.mean([np.mean(v) for v in transects.values()]))
        matched.append(dict(species=species,MPA=mpa,status=m['status'],site_id=site,site_name=m['geolocation'],year=year,n_count_rows=len(finite),n_missing_count_rows=len(rs)-len(finite),sum_observed_counts=sum(number(r['count']) for r in finite),count_density=density,n_count_dates=len(dates),count_dates=';'.join(sorted(dates)),has_size_pool=bool(po),n_sized=po['n_sized'] if po else 0,exact_date_match=bool(sd&dates),size_dates_missing=bool(po) and not bool(sd),all_size_dates_match_count=bool(sd) and sd<=dates,n_sized_juvenile=po['n_juvenile'] if po else 0,conditional_relative_halfwidth=po['conditional_relative_halfwidth'] if po else '',before_designation=year<int(m['MPA_Start'])))
assert len(matched)==len({(r['species'],r['MPA'],r['status'],r['site_id'],r['year']) for r in matched})
write('site_year_matches',matched)
# The fish site membership/weights must reproduce the primary annual density values.
canonical=read(ROOT/'data/harmonized/harmonized_raw_responses.csv')
verification=[]
for r in canonical:
    species=short.get(r['taxon_name']);mpa=r['CA_MPA_Name_Short'];year=int(r['year'])
    if r['source']!='KFM' or species!='S. pulcher' or r['resp']!='Den' or (mpa,species) not in retained:continue
    rows=[x for x in matched if x['species']==species and x['MPA']==mpa and x['year']==year and x['status']==r['status']]
    if not rows:continue
    rebuilt=float(np.mean([x['count_density'] for x in rows]));err=rebuilt-float(r['value'])
    assert abs(err)<1e-10,(mpa,year,r['status'],rebuilt,r['value'])
    verification.append(dict(MPA=mpa,year=year,status=r['status'],primary_density=r['value'],rebuilt_density=rebuilt,difference=err))
write('count_reconstruction',verification)

summary=[]
for species,mpa in sorted({(r['species'],r['MPA']) for r in matched}):
    rs=[r for r in matched if r['species']==species and r['MPA']==mpa]
    complete=[];anypair=[];same_date=[];sizes_years=[];positive_complete=[]
    for year in sorted({r['year'] for r in rs}):
        a=[r for r in rs if r['year']==year];statuses={r['status'] for r in a};assert statuses=={'mpa','reference'},(species,mpa,year,statuses)
        if {r['status'] for r in a if r['has_size_pool']}==statuses:anypair.append(year)
        if all(r['has_size_pool'] for r in a):complete.append(year)
        if any(r['has_size_pool'] for r in a) and all(r['has_size_pool'] or r['sum_observed_counts']==0 for r in a):positive_complete.append(year)
        if all(r['has_size_pool'] and r['all_size_dates_match_count'] for r in a):same_date.append(year)
        if any(r['has_size_pool'] for r in a):sizes_years.append(year)
    withsize=[r for r in rs if r['has_size_pool']]
    summary.append(dict(species=species,MPA=mpa,count_site_years=len(rs),positive_count_site_years=sum(r['sum_observed_counts']>0 for r in rs),size_matched_site_years=len(withsize),positive_count_site_years_without_sizes=sum(r['sum_observed_counts']>0 and not r['has_size_pool'] for r in rs),count_years=len(set(r['year'] for r in rs)),years_any_size_each_side=';'.join(map(str,anypair)),years_all_count_sites_have_sizes=';'.join(map(str,complete)),years_all_positive_count_sites_have_sizes=';'.join(map(str,positive_complete)),years_all_size_dates_match_counts=';'.join(map(str,same_date)),n_before_size_matches=sum(r['before_designation'] and r['has_size_pool'] for r in rs),matched_pools_under10=sum(r['n_sized']<10 for r in withsize),matched_pools_under20=sum(r['n_sized']<20 for r in withsize)))
write('coverage',summary)
INPUTS += [ROOT/'code/R/01_utils.R',Path(__file__)]
write('inputs',[dict(path=str(p.resolve()),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in dict.fromkeys(INPUTS)])
metadata=dict(legacy_sheephead_rows=len(selected),legacy_identical_rows_retained=exact_duplicates,legacy_invert_rows_checked=len(inv),legacy_lobster_size_rows=0,nps2015_lobster_size_rows=0,matched_sites=len({(r['species'],r['site_id']) for r in matched}),no_primary_changes=True,seed=20261029,bootstrap_repetitions=2000)
(OUT/'table_s_kfm_size_audit.json').write_text(json.dumps(metadata,indent=2))
print(json.dumps(metadata,indent=2));print(json.dumps(summary,indent=2))
