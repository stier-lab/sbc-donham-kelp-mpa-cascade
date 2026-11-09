from pathlib import Path
from collections import Counter
from decimal import Decimal,InvalidOperation
import pandas as pd,json,hashlib
B=Path('/Users/adrianstier/Library/CloudStorage/GoogleDrive-<account>/My Drive/Stier Lab/People/Emily Donham/Projects/Kelp MPA')
O=Path('outputs/vrg_source_search_20261102');old=pd.read_csv('data/ALL_sizefreq_2024.csv',low_memory=False);v=old[old.campus.eq('VRG')&old.classcode.isin(['MESFRA','STRPUR','PANINT'])].copy();w=pd.read_csv(O/'workspace_SizeFreq.south.csv',low_memory=False)
cols=list(v.columns)
def key(z):
 q=z[cols].copy()
 for c in cols:
  numeric=c in ['survey_year','year','month','day','transect','size','depth','count']
  def norm(x):
   if pd.isna(x) or x=='':return '<NA>'
   if numeric:
    try:return str(Decimal(str(x)).normalize())
    except InvalidOperation:pass
   return str(x)
  q[c]=q[c].map(norm)
 return Counter(map(tuple,q.values))
a=key(v);b=key(w);assert not b-a,'Workspace has unexplained new/different rows'
missing=a-b;rows=[]
for k,n in missing.items():rows.append(dict(zip(cols,k))|{'multiplicity':n})
lost=pd.DataFrame(rows);lost.to_csv(O/'workspace_excluded_records.csv',index=False)
# Record source copy agreement; it verifies preservation, not independent original measurements.
paths=[Path('data/ALL_sizefreq_2024.csv'),B/'data/ALL_sizefreq_2024.csv',B/'Github/MPA_Kelp/data/ALL_sizefreq_2024.csv']
files=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),meaning='combined derived file') for p in paths];assert len({r['sha256'] for r in files})==1
files+= [dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),meaning='archived analysis evidence') for p in [B/'Github/MPA_Kelp/.RData',B/'Github/MPA_Kelp/.Rhistory',B/'Github/MPA_Kelp/code/pBACIPS_PISCO_V10.R',B/'MPA_code_notes.docx']]
pd.DataFrame(files).to_csv(O/'source_inventory.csv',index=False)
v.groupby(['classcode','method']).agg(rows=('count','size'),animals=('count','sum'),first_year=('year','min'),last_year=('year','max'),sites=('site','nunique')).reset_index().to_csv(O/'vrg_coverage.csv',index=False)
lost.groupby(['site','classcode']).multiplicity.sum().reset_index().to_csv(O/'workspace_exclusions.csv',index=False)
# Verify remaining alias assignments against curated physical-site metadata without changing them.
sites=pd.read_csv(B/'Github/MPA_Kelp/data/PISCO/master_site_table_Emilyedit.csv');cross=sites[['site','campus','site_new']].drop_duplicates();assert not cross.duplicated(['site','campus']).any()
z=v.rename(columns={'site_new':'old_site_new'}).merge(cross,on=['site','campus'],how='left',validate='many_to_one');bad=z[z.site_new.isna()|z.old_site_new.ne(z.site_new)];bad.to_csv(O/'vrg_alias_unresolved.csv',index=False)
r=dict(focal_rows=len(v),focal_animals=int(v['count'].sum()),workspace_rows=len(w),workspace_exact_matches=sum(b.values()),removed_rows=sum(missing.values()),removed_sites=lost.site.unique().tolist(),alias_unresolved_rows=len(bad),independent_original_found=False,assembly_script_found=False,adopted=False)
(O/'verification.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2));print(lost.groupby('site').multiplicity.sum().to_string())
