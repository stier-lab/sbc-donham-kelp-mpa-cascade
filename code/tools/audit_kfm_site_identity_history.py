from pathlib import Path
import pandas as pd,numpy as np,itertools,json,hashlib,subprocess
O=Path('outputs/kfm_git_origin_20261102');O.mkdir(exist_ok=True)
p=pd.read_csv('data/ALL_sizefreq_2024.csv',low_memory=False);s=pd.read_csv('data/external/scott_gabara_2026-11-02/normalized_source.csv');s['year']=pd.to_datetime(s.startDate).dt.year
sp=['Strongylocentrotus franciscanus','Strongylocentrotus purpuratus'];p=p[p.campus.eq('KFM')&p.classcode.isin(sp)&p.year.isin(s.year.unique())]
geo=pd.read_csv('data/MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv').set_index('site_id');rows=[]
for (site,year,tax),g in p.groupby(['site','year','classcode']):
 target=g.groupby('size')['count'].sum();z=s[s.year.eq(year)&s.scientificName.eq(tax)];hist={n:x.groupby('invertebrateSizeInMillimeters')['count'].sum() for n,x in z.groupby('siteName')}
 def eq(x):return target.subtract(x,fill_value=0).abs().sum()==0
 matches=[n for n,x in hist.items() if eq(x)];pairs=[]
 if not matches:
  for n,m in itertools.combinations(hist,2):
   if eq(hist[n].add(hist[m],fill_value=0)):pairs.append(n+' + '+m)
 rows.append(dict(site=site,year=year,taxon=tax,expected_name=';'.join(sorted(set(geo.loc[[site],'geolocation'].astype(str)))) if site in geo.index else '',matched_source=';'.join(matches or pairs),match_type='single' if matches else 'combined' if pairs else 'unresolved',n_matches=len(matches or pairs),n_animals=g['count'].sum()))
a=pd.DataFrame(rows)
assert len(a)==558 and a.n_matches.eq(1).all()
aliases={'Gull Island':'Gull Island South','Black Seabass Reef':'Black Sea Bass Reef'}
a['identity_status']=a.apply(lambda r:'known source identity error' if r.site in ['a-k-01','a-k-14'] else 'verified matching physical site' if aliases.get(r.expected_name,r.expected_name)==r.matched_source else 'NEW mismatch: investigate',axis=1)
assert not a.identity_status.str.startswith('NEW').any()
a.to_csv(O/'all_kfm_site_year_matches.csv',index=False)
a.groupby(['site','expected_name','matched_source','match_type'],dropna=False).agg(site_species_years=('year','size'),first_year=('year','min'),last_year=('year','max')).reset_index().to_csv(O/'site_identity_summary.csv',index=False)
root=Path('/Users/adrianstier/Library/CloudStorage/GoogleDrive-<account>/My Drive/Stier Lab/People/Emily Donham/Projects/Kelp MPA');copies=[]
for f in [Path('data/ALL_sizefreq_2024.csv'),root/'Github/MPA_Kelp/data/ALL_sizefreq_2024.csv']:
 q=pd.read_csv(f,low_memory=False);k=q[q.campus.eq('KFM')];copies.append(dict(file=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest(),rows=len(q),kfm_rows=len(k),kfm_hist_sha256=hashlib.sha256(k.groupby(['site','year','classcode','size'])['count'].sum().sort_index().to_csv().encode()).hexdigest()))
pd.DataFrame(copies).to_csv(O/'source_copy_comparison.csv',index=False)
print(a.groupby('match_type').size().to_string());print(pd.read_csv(O/'site_identity_summary.csv').to_string(index=False));print(pd.DataFrame(copies).to_string(index=False))
