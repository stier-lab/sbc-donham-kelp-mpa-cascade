from pathlib import Path
import pandas as pd,numpy as np,itertools
O=Path('outputs/kfm_pipeline_guard_20261102')
p=pd.read_csv('data/ALL_sizefreq_2024.csv',low_memory=False);legacy=pd.read_csv('/Users/adrianstier/sbc-donham-kelp-mpa-data/data/source_archive/kfm_size_feasibility_2026-10-29/legacy/Invertebrate Size Frequency.csv');legacy['mm']=(legacy.Size*10).round(6)
p=p[p.campus.eq('KFM')&p.classcode.isin(['Strongylocentrotus franciscanus','Strongylocentrotus purpuratus'])&p.year.le(2011)];rows=[]
for (site,year,tax),g in p.groupby(['site','year','classcode']):
 code=11005 if 'franciscanus' in tax else 11006;z=legacy[legacy.Year.eq(year)&legacy.Species.eq(code)];hist={int(n):x.groupby('mm').size() for n,x in z.groupby('Site')};t=g.groupby('size')['count'].sum();eq=lambda x:t.subtract(x,fill_value=0).abs().sum()==0
 matches=[str(n) for n,x in hist.items() if eq(x)]
 if not matches:
  matches=[f'{n}+{m}' for n,m in itertools.combinations(hist,2) if eq(hist[n].add(hist[m],fill_value=0))]
 expected=str(int(site.split('-')[-1]));status='verified intended site' if matches==[expected] else 'known combined 1+21' if site=='a-k-01' and matches==['1+21'] else 'known substituted 37' if site=='a-k-14' and matches==['37'] else 'UNRESOLVED'
 rows.append(dict(site=site,year=year,taxon=tax,source_sites=';'.join(matches),status=status,n_animals=g['count'].sum()))
a=pd.DataFrame(rows);a.to_csv(O/'older_identity_audit.csv',index=False);print(a.groupby('status').size().to_string());print(a[a.status.eq('UNRESOLVED')].to_string(index=False))
