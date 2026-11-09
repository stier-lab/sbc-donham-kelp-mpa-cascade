"""Independent, count-preserving comparison of raw and verified KFM urchin sizes."""
from pathlib import Path
from collections import Counter
import pandas as pd
import hashlib,json
O=Path('outputs/kfm_source_comparison_20261102'); O.mkdir(exist_ok=True)
raw_path=Path('data/ALL_sizefreq_2024.csv')
legacy_path=Path('/Users/adrianstier/sbc-donham-kelp-mpa-data/data/source_archive/kfm_size_feasibility_2026-10-29/legacy/Invertebrate Size Frequency.csv')
scott_path=Path('data/external/scott_gabara_2026-11-02/KFM Nat Hab Red Purp 2012-2022.xlsx')
bundle=Path('/Users/adrianstier/sbc-donham-kelp-mpa-data/data/verified_kfm_sizes')
raw=pd.read_csv(raw_path,low_memory=False);legacy=pd.read_csv(legacy_path);scott=pd.read_excel(scott_path)
normalized=pd.read_csv(scott_path.parent/'normalized_source.csv')
# Read original workbook independently; confirm its normalization has preserved all physical fields.
cols=['siteName','scientificName','invertebrateSizeInMillimeters','count']
pd.testing.assert_frame_equal(scott[cols].reset_index(drop=True),normalized[cols].reset_index(drop=True),check_dtype=False)
for c in ['startDate','endDate']:assert pd.to_datetime(scott[c]).equals(pd.to_datetime(normalized[c]))
scott['year']=pd.to_datetime(scott.startDate).dt.year
names={11005:'Strongylocentrotus franciscanus',11006:'Strongylocentrotus purpuratus'}
geo=pd.read_csv('data/MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv')[['site_id','geolocation']].drop_duplicates().set_index('site_id').geolocation.to_dict()
alias={'Gull Island':'Gull Island South','Black Seabass Reef':'Black Sea Bass Reef','Southeast Sea Lion':'SE Sea Lion Rookery','South East Reef':'Southeast Reef'}
lookup={alias.get(v,v):k for k,v in geo.items() if k.startswith('a-k-')}
legacy=legacy[legacy.Species.isin(names)].copy();legacy['site']=legacy.Site.map(lambda n:f'a-k-{int(n):02d}');legacy['taxon']=legacy.Species.map(names);legacy['mm']=(legacy.Size*10).round(6)
scott['site']=scott.siteName.map(lookup)
assert not scott.site.isna().any(),scott.loc[scott.site.isna(),'siteName'].unique()
# Repeated identical measurements are distinct animals. Never deduplicate them.
def groups(df,keys,size,count=None):
 out={}
 for key,g in df.groupby(keys,dropna=False):
  assert g[size].notna().all()
  c=Counter()
  for v,n in zip(g[size],g[count] if count else [1]*len(g)):
   assert pd.notna(n) and n>=0 and int(n)==n
   c[round(float(v),6)]+=int(n)
  out[key]=c
 return out
L=groups(legacy,['site','Year','taxon'],'mm')
S=groups(scott,['site','year','scientificName'],'invertebrateSizeInMillimeters','count')
assert not set(L)&set(S)
sources={**L,**S}
k=raw[raw.campus.eq('KFM')&raw.classcode.isin(names.values())].copy()
R=groups(k,['site','year','classcode'],'size','count')
rep=pd.read_csv(bundle/'replacement.csv');P=groups(rep,['site','year','classcode'],'size','count')
registry=pd.read_csv(bundle/'registry.csv');V={(r.site,r.year,r.classcode):r.verified for r in registry.itertuples()}
rows=[];diff=[]
for key,r in R.items():
 site,year,tax=key;s=sources.get(key);e=P.get(key,r)
 def distance(a,b):return sum(abs(a.get(v,0)-b.get(v,0)) for v in set(a)|set(b))
 origin='legacy' if year<=2011 else 'Scott' if year in set(scott.year) else 'no independent source supplied'
 record=dict(site=site,site_name=geo.get(site,''),year=year,taxon=tax,source=origin,raw_animals=sum(r.values()),source_animals=sum(s.values()) if s is not None else None,raw_exact=s is not None and r==s,effective_exact=s is not None and e==s,effective_replaced=key in P,registry_verified=V[key],raw_absolute_frequency_difference=distance(r,s) if s is not None else None,effective_absolute_frequency_difference=distance(e,s) if s is not None else None)
 rows.append(record)
 if s is not None and r!=s:
  for v in sorted(set(r)|set(s)):
   if r.get(v,0)!=s.get(v,0):diff.append(dict(site=site,year=year,taxon=tax,size_mm=v,raw_animals=r.get(v,0),source_animals=s.get(v,0)))
a=pd.DataFrame(rows);a.to_csv(O/'all_group_comparisons.csv',index=False);pd.DataFrame(diff).to_csv(O/'size_frequency_differences.csv',index=False)
missing=[]
for key,s in sources.items():
 if key not in R:
  site,y,tax=key
  missing.append(dict(site=site,site_name=geo.get(site,''),year=y,taxon=tax,source='legacy' if key in L else 'Scott',animals=sum(s.values()),interpretation='Source-only group; absence is not itself an error or authorization to include'))
m=pd.DataFrame(missing)
def absence(r):
 if r.site=='a-k-21':return 'Excluded Miracle Mile; some measurements occur under Wyckoff in raw file'
 if r.site=='a-k-37':return 'Excluded Southeast Reef; measurements occur under site 14 in raw file'
 if r.site=='a-k-14':return 'Original SE Sea Lion 1985–2004 groups absent from combined file; no automatic addition'
 return 'Separate source-site groups absent from combined file; retain existing eligibility decisions'
m['interpretation']=m.apply(absence,axis=1);m.to_csv(O/'source_only_groups.csv',index=False)
m.groupby(['source','site','site_name','interpretation']).agg(groups=('year','size'),first_year=('year','min'),last_year=('year','max')).reset_index().to_csv(O/'source_only_summary.csv',index=False)
coverage=raw[raw.campus.eq('KFM')].groupby('classcode').agg(rows=('size','size'),first_year=('year','min'),last_year=('year','max')).reset_index();coverage['individual_source_comparison']=coverage.classcode.isin(names.values());coverage.to_csv(O/'taxon_scope.csv',index=False)
# No verified effective group may disagree with its independent original.
assert a.loc[a.registry_verified,'effective_exact'].all()
summary=a.groupby('source').agg(groups=('site','size'),raw_exact=('raw_exact','sum'),effective_exact=('effective_exact','sum'),replacements=('effective_replaced','sum')).reset_index();summary.to_csv(O/'summary.csv',index=False)
protected=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text());assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in protected.items())
verification=dict(workbook_rows=len(scott),workbook_normalization_checked=True,all_groups=len(a),source_only_groups=len(missing),verified_effective_groups=int(a.registry_verified.sum()),unresolved_groups=int((~a.registry_verified).sum()),protected_files_unchanged=len(protected),sources={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [raw_path,legacy_path,scott_path,bundle/'replacement.csv',bundle/'registry.csv']})
(O/'verification.json').write_text(json.dumps(verification,indent=2));print(summary.to_string(index=False));print('Source-only groups',len(missing));print('Effective discrepancies with supplied sources:');print(a[(a.source!='no independent source supplied')&~a.effective_exact][['site','year','taxon','registry_verified']].to_string(index=False))
