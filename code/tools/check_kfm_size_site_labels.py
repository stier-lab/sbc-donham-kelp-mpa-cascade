from pathlib import Path
import pandas as pd,numpy as np
D=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data');O=Path('outputs/biomass_protocols_20261102')
s=pd.read_csv(D/'ALL_sizefreq_2024.csv',low_memory=False);s=s[s.campus.eq('KFM') & s.classcode.isin(['Strongylocentrotus franciscanus','Strongylocentrotus purpuratus'])].copy();s['code']=s.classcode.map({'Strongylocentrotus franciscanus':11005,'Strongylocentrotus purpuratus':11006})
l=pd.read_csv(D/'source_archive/kfm_size_feasibility_2026-10-29/legacy/Invertebrate Size Frequency.csv');l=l[l.Species.isin([11005,11006])].copy();l['mm']=(l.Size*10).round(6)
rows=[]
for (site,year,code),g in s[s.site.isin(['a-k-01','a-k-14'])&s.year.le(2011)].groupby(['site','year','code']):
 t=g.groupby('size')['count'].sum();q=l[l.Year.eq(year)&l.Species.eq(code)]
 choices=[('labeled_site',[int(site[-2:])])]+([('Wyckoff plus Miracle Mile',[1,21])] if site=='a-k-01' else [('Southeast Reef',[37])])
 for label,sites in choices:
  h=q[q.Site.isin(sites)].groupby('mm').size();err=t.subtract(h,fill_value=0).abs().sum();rows.append(dict(combined_site=site,year=year,species_code=code,candidate=label,legacy_sites=';'.join(map(str,sites)),combined_animals=int(t.sum()),legacy_animals=int(h.sum()),absolute_histogram_difference=float(err),exact_match=err==0))
pd.DataFrame(rows).to_csv(O/'kfm_label_hypotheses.csv',index=False)
print(pd.DataFrame(rows).groupby(['combined_site','candidate']).agg(groups=('year','size'),exact=('exact_match','sum')).to_string())
