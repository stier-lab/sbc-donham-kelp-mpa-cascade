from pathlib import Path
import os,hashlib,json,zipfile
import pandas as pd
O=Path('outputs/kfm_holdings_20261102')
roots=[Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data'),Path('/Users/adrianstier/sbc-donham-kelp-mpa-cascade/data'),Path('/Users/adrianstier/kumagai2024-comparison'),Path('/Users/adrianstier/Library/CloudStorage/GoogleDrive-<account>/My Drive/Stier Lab/People/Emily Donham/Projects/Kelp MPA'),Path('/Users/adrianstier/sbc-donham-kelp-mpa-manuscript/correspondence'),Path('/Users/adrianstier/Documents/ChatGPT/SBC-Data-Catalogue')]
paths=[]
for root in roots:
 for base,dirs,files in os.walk(root):
  dirs[:]=[d for d in dirs if d not in ['.git','node_modules','.venv']]
  for fn in files:
   p=Path(base)/fn;low=str(p).lower()
   if p.suffix.lower() in ['.csv','.xlsx','.zip','.invalid','.rds','.accdb','.mdb'] and any(t in low for t in ['kfm','nhsf','all_sizefreq','sbcmbon','natural_habitat']):paths.append(p)
rows=[];seen={}
for p in paths:
 try:
  target=p.resolve();size=p.stat().st_size
  if target not in seen:seen[target]=hashlib.sha256(p.read_bytes()).hexdigest()
  rows.append(dict(path=str(p),resolved_path=str(target),bytes=size,sha256=seen[target],kind=p.suffix.lower(),valid_archive=zipfile.is_zipfile(p) if p.suffix in ['.zip','.invalid'] else None))
 except Exception as e:rows.append(dict(path=str(p),error=str(e)))
pd.DataFrame(rows).to_csv(O/'file_inventory.csv',index=False)
# Inspect distinct raw tables; cached estimates are not new source measurements.
summaries=[]
for target,h in seen.items():
 if target.suffix.lower() not in ['.csv','.xlsx']:continue
 try:
  if target.suffix=='.xlsx':sheets=pd.read_excel(target,sheet_name=None)
  else:
   try:frame=pd.read_csv(target,low_memory=False)
   except UnicodeDecodeError:frame=pd.read_csv(target,low_memory=False,encoding='latin1')
   sheets={'csv':frame}
  for sheet,x in sheets.items():
   src=next((c for c in ['campus','data_source'] if c in x),None)
   if src:x=x[x[src].astype(str).str.lower().eq('kfm')]
   yearcol=next((c for c in ['year','Year','SurveyYear','survey_year'] if c in x),None)
   datecol=next((c for c in ['date','Date','SurveyDate'] if c in x),None)
   years=pd.to_numeric(x[yearcol],errors='coerce') if yearcol else pd.to_datetime(x[datecol],errors='coerce',format='mixed').dt.year if datecol else pd.Series(dtype=float)
   summaries.append(dict(path=str(target),sheet=sheet,kfm_rows=len(x),first_year=years.min(),last_year=years.max(),has_date_column=bool(datecol),has_site_column=any('site' in c.lower() for c in x),columns=';'.join(x.columns)))
 except Exception as e:summaries.append(dict(path=str(target),error=str(e)))
pd.DataFrame(summaries).to_csv(O/'table_inventory.csv',index=False)
print('Files',len(rows),'resolved targets',len(seen),'content hashes',len(set(seen.values())),'table summaries',len(summaries))
print(pd.DataFrame(summaries)[['path','kfm_rows','first_year','last_year']].to_string(index=False))
(O/'scope.json').write_text(json.dumps({'roots':[str(x) for x in roots],'file_paths':len(rows),'resolved_targets':len(seen),'unique_content_hashes':len(set(seen.values()))},indent=2))
