from pathlib import Path
from collections import Counter
import pandas as pd,numpy as np,gzip,hashlib,json,xml.etree.ElementTree as ET
O=Path('outputs/lter_original_sizes_20261102');O.mkdir(exist_ok=True)
f=Path('/Users/adrianstier/sbc-sswd-kelp/data/source_archive/sbc_grit/candidate_data/Annual_Quad_Swath_All_Years_20260923.csv.gz')
meta=Path('/Users/adrianstier/Documents/ChatGPT/SBC-Data-Catalogue/deep-inventory/raw/packages/knb-lter-sbc/19/33/metadata.xml')
b=gzip.decompress(f.read_bytes());expected=ET.parse(meta).find('.//authentication').text;assert hashlib.md5(b).hexdigest()==expected
raw=pd.read_csv(f);raw['source_row']=np.arange(len(raw))+1
old=pd.read_csv('data/ALL_sizefreq_2024.csv',low_memory=False);old=old[old.campus.eq('LTER')&old.classcode.isin(['SFL','SFS','SPL','SPS'])].copy()
# Site names independently decoded in MBON's provider-integrated count table.
site_map={};conflict=False
for chunk in pd.read_csv('data/MBON/SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv',usecols=['data_source','site_id','site_name'],chunksize=200000,low_memory=False):
 for a,c in chunk.loc[chunk.data_source.eq('lter'),['site_id','site_name']].drop_duplicates().itertuples(index=False,name=None):
  assert a not in site_map or site_map[a]==c;site_map[a]=c
pd.DataFrame([{'site':k,'original_site':v} for k,v in site_map.items()]).to_csv(O/'site_crosswalk.csv',index=False)
new=raw[raw.SITE.isin([site_map[k] for k in old.site.unique()])&raw.YEAR.le(old.year.max())&raw.SP_CODE.isin(old.classcode.unique())&raw.COUNT.gt(0)].copy()
new['site']=new.SITE.map({v:k for k,v in site_map.items()});new['year']=new.YEAR;new['month']=new.MONTH;new['day']=pd.to_datetime(new.DATE).dt.day;new['classcode']=new.SP_CODE;new['count']=new.COUNT
new['size']=(new.SIZE.where(new.SIZE.ne(-99999))*10).round(6)
old['size']=old['size'].round(6)
keys=['site','year','month','day','classcode','count']
def counter(d,cols):return Counter(tuple('NA' if pd.isna(v) else v for v in row) for row in d[cols].itertuples(index=False,name=None))
a=counter(old,keys);c=counter(new,keys)
print('Counts/date keys',len(old),len(new),'old-only',sum((a-c).values()),'new-only',sum((c-a).values()))
rows=[]
for code in ['SFL','SFS','SPL','SPS']:
 p=old[old.classcode.eq(code)];q=new[new.classcode.eq(code)];pc=counter(p,keys+['size']);qc=counter(q,keys+['size'])
 rows.append(dict(code=code,combined_rows=len(p),original_rows=len(q),combined_valid_size=int(p['size'].gt(0).sum()),original_valid_size=int(q['size'].gt(0).sum()),exact_rows=sum((pc&qc).values()),combined_unique_sizes=';'.join(map(str,sorted(p['size'].dropna().unique()))) if code in ['SFS','SPS'] else 'see row table'))
small=new.classcode.isin(['SFS','SPS'])
pc=counter(old,keys+['size']);qc=counter(new,keys+['size'])
print(pd.DataFrame(rows).to_string(index=False));print('After cm-to-mm conversion and rounding to6decimals: exact',sum((pc&qc).values()),'old-only',sum((pc-qc).values()),'new-only',sum((qc-pc).values()))
pd.DataFrame(rows).to_csv(O/'comparison.csv',index=False)
# Preserve exact discrepancies; no fabricated one-to-one specimen matching.
for name,diff in [('combined_unmatched',pc-qc),('original_unmatched',qc-pc)]:
 pd.DataFrame([dict(zip(keys+['size','multiplicity'],list(k)+[v])) for k,v in diff.items()],columns=keys+['size','multiplicity']).to_csv(O/(name+'.csv'),index=False)
new.groupby(['site','classcode']).agg(rows=('source_row','size'),measured_rows=('size',lambda x:x.gt(0).sum()),first_year=('year','min'),last_year=('year','max')).reset_index().to_csv(O/'source_support.csv',index=False)
verification=dict(package='knb-lter-sbc.19.33',source=str(f),source_sha256=hashlib.sha256(b).hexdigest(),compressed_sha256=hashlib.sha256(f.read_bytes()).hexdigest(),eml_md5_verified=True,source_rows=len(raw),selected_source_rows=len(new),combined_rows=len(old),counts_and_date_keys_equal=a==c,exact_rows_without_proxy=sum((counter(old,keys+['size'])&counter(new,keys+['size'])).values()),exact_rows_with_dates=sum((pc&qc).values()),small_original_missing_sizes=int((small&new['size'].isna()).sum()),source_unit='cm; average size of individuals counted, not necessarily individual measurements',source_small_sizes='2.49cm values already occur in original source; no added proxy in this comparison')
verification['all_rows_equal_without_day']=counter(old,[c for c in keys if c!='day']+['size'])==counter(new,[c for c in keys if c!='day']+['size'])
assert verification['all_rows_equal_without_day']
(O/'verification.json').write_text(json.dumps(verification,indent=2))
