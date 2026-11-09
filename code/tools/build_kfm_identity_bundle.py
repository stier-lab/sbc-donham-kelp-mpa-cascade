from pathlib import Path
import pandas as pd,numpy as np,hashlib,json
O=Path('outputs/kfm_pipeline_guard_20261102');B=O/'private/bundle';B.mkdir(parents=True,exist_ok=True)
raw=pd.read_csv('data/ALL_sizefreq_2024.csv',low_memory=False);k=raw[raw.campus.eq('KFM')&raw.classcode.isin(['Strongylocentrotus franciscanus','Strongylocentrotus purpuratus'])]
old=pd.read_csv(O/'older_identity_audit.csv');recent=pd.read_csv('outputs/kfm_git_origin_20261102/all_kfm_site_year_matches.csv');keys=set(zip(old.loc[~old.status.eq('UNRESOLVED'),'site'],old.loc[~old.status.eq('UNRESOLVED'),'year'],old.loc[~old.status.eq('UNRESOLVED'),'taxon']))|set(zip(recent.site,recent.year,recent.taxon))
l=pd.read_csv('/Users/adrianstier/sbc-donham-kelp-mpa-data/data/source_archive/kfm_size_feasibility_2026-10-29/legacy/Invertebrate Size Frequency.csv');l['row']=np.arange(len(l))+1;s=pd.read_csv('data/external/scott_gabara_2026-11-02/normalized_source.csv');s['year']=pd.to_datetime(s.startDate).dt.year
resolution_path=Path('outputs/kfm_remaining_sources_20261102/hare_1987_resolution.json')
resolved_hare=False
if resolution_path.exists():
 proof=json.loads(resolution_path.read_text());assert hashlib.sha256(Path(proof['source_file']).read_bytes()).hexdigest()==proof['source_sha256'];assert proof['n']==102 and proof['bin_checks']==20
 keys.add(('a-k-02',1987,'Strongylocentrotus purpuratus'));resolved_hare=True
reg=[];repl=[]
for (site,year,tax),g in k.groupby(['site','year','classcode']):
 verified=(site,year,tax) in keys;reason='Exact source histogram checked' if verified else 'Original source mismatch or year outside supplied originals'
 if verified and (site in ['a-k-01','a-k-14'] or (resolved_hare and site=='a-k-02' and year==1987 and tax=='Strongylocentrotus purpuratus')):
  if year<=2011:
   z=l[l.Site.eq(int(site[-2:]))&l.Year.eq(year)&l.Species.eq(11005 if 'franciscanus' in tax else 11006)]
   records=[(r.Size*10,1,pd.to_datetime(r.Date,format='%d%b%Y'),f'legacy:{r.row}') for r in z.itertuples()]
  else:
   z=s[s.year.eq(year)&s.scientificName.eq(tax)&s.siteName.eq('Wyckoff Ledge' if site=='a-k-01' else 'SE Sea Lion Rookery')]
   records=[(r.invertebrateSizeInMillimeters,r.count,pd.Timestamp(r.startDate),f'scott:{r.excel_row}') for r in z.itertuples()]
  if not records:verified=False;reason='Allowed physical-site replacement unavailable'
  for size,count,date,rid in records:
   d={col:None for col in raw.columns};d.update(campus='KFM',method='original natural habitat size survey',survey_year=year,year=year,month=date.month,day=date.day,site=site,site_new=site,classcode=tax,count=count,size=size,source_record_id=rid);repl.append(d)
 if resolved_hare and site=='a-k-02' and year==1987 and tax=='Strongylocentrotus purpuratus':reason='Legacy102records corroborated by NPS published20-bin table; replaced combined103records'
 reg.append(dict(site=site,year=year,classcode=tax,verified=verified,reason=reason))
r=pd.DataFrame(reg);r.to_csv(B/'registry.csv',index=False);pd.DataFrame(repl).to_csv(B/'replacement.csv',index=False);r[~r.verified].to_csv(O/'unresolved_source_groups.csv',index=False)
meta=dict(source_md5=hashlib.md5(Path('data/ALL_sizefreq_2024.csv').read_bytes()).hexdigest(),files={n:hashlib.md5((B/n).read_bytes()).hexdigest() for n in ['registry.csv','replacement.csv']});(B/'manifest.json').write_text(json.dumps(meta));print('registry',len(r),'verified',r.verified.sum(),'unverified',(~r.verified).sum(),'replacement records',len(repl));print(r[~r.verified].groupby(['site','year']).size().to_string())

# Install only the verified bundle, never overwrite raw files or numerical products.
import subprocess,shutil
subprocess.run(['Rscript','-e',f'm<-jsonlite::fromJSON("{B}/manifest.json");saveRDS(m,"{B}/manifest.rds")'],check=True)
dest=Path('/Users/adrianstier/sbc-donham-kelp-mpa-data/data/verified_kfm_sizes');dest.mkdir(exist_ok=True)
for f in B.iterdir():shutil.copy2(f,dest/f.name)
