from pathlib import Path
from collections import Counter
from decimal import Decimal, InvalidOperation
import pandas as pd,hashlib,json,xml.etree.ElementTree as ET
O=Path('outputs/lter_integration_pisco_check_20261102');O.mkdir(exist_ok=True)
base=Path('/Users/adrianstier/Documents/ChatGPT/SBC-Data-Catalogue/pisco-data-catalogue/data');p=base/'raw/PISCO_kelpforest.1.11/PISCO_kelpforest_sizefreq.1.11.csv';sm=base/'metadata/PISCO_kelpforest.1.11/PISCO_kelpforest_sizefreq.1.11.csv.system_metadata.xml'
checks=ET.parse(sm).findall('.//{*}checksum');digest=hashlib.sha256(p.read_bytes()).hexdigest();assert any(x.text==digest for x in checks)
c=pd.read_csv(p,dtype=str,keep_default_na=False);old=pd.read_csv('data/ALL_sizefreq_2024.csv',dtype=str,keep_default_na=False)
focal=['MESFRA','STRPUR','PANINT'];o=old[old.campus.eq('UCSB')&old.method.eq('SBTL_SIZEFREQ_PISCO')&old.classcode.isin(focal)].copy();n=c[c.campus.eq('UCSB')&c.classcode.isin(focal)].copy()
extra=n[pd.to_numeric(n.year).gt(2021)].copy();n=n[pd.to_numeric(n.year).le(2021)].copy()
cols=list(c.columns)
num=['survey_year','year','month','day','transect','count','size','depth']
def normalize(x,original=False):
 x=x.copy()
 if original:
  x['size']=x['size'].map(lambda z: str(Decimal(z)*10) if z not in ['','NA'] else z)
  x['method']=x.method.replace({'SBTL_ROVINGDIVERSIZEFREQ_PISCO':'SBTL_SIZEFREQ_PISCO'})
 x['zone']=x.zone.replace({'':'NA'})
 def numeric(v):
  try:return str(Decimal(v).normalize())
  except InvalidOperation:return v
 for col in num:x[col]=x[col].map(numeric)
 return x
oc=Counter(map(tuple,normalize(o)[cols].values));nc=Counter(map(tuple,normalize(n,True)[cols].values))
assert oc==nc,(sum((oc-nc).values()),sum((nc-oc).values()))
summary=o.groupby('classcode').agg(rows=('count','size'));summary['matched_rows']=summary['rows'];summary.to_csv(O/'pisco_exact_matches.csv')
extra.groupby(['year','classcode']).agg(rows=('count','size')).reset_index().to_csv(O/'pisco_unadopted_extensions.csv',index=False)
alias=o[o.site.eq('SCI_PAINTED_CAVE_E')&o.site_new.eq('SRI_SOUTH_POINT')];assert len(alias)==99
alias.groupby(['site','site_new','year','classcode']).agg(rows=('count','size')).reset_index().to_csv(O/'pisco_alias_conflict.csv',index=False)
report=dict(source=str(p),sha256=digest,system_metadata_checksum_verified=True,raw_rows=len(c),old_roving_rows=len(o),exact_matches=sum(oc.values()),unadopted_new_rows=len(extra),painted_cave_wrong_aggregation_alias_rows=len(alias),vrg_present=bool(c.campus.eq('VRG').any()),normalized_fields=['source cm to mm','method name','blank zone to NA','equivalent numeric serialization'],physical_source_identity_preserved=True,adopted=False)
(O/'pisco_verification.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
