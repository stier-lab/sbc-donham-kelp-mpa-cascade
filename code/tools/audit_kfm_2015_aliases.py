from pathlib import Path
import pandas as pd,numpy as np,re,json,hashlib
O=Path('outputs/biomass_exclusion_reconciliation_20261102');D=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data')
x=pd.read_excel(D/'source_archive/kfm_size_feasibility_2026-10-29/2015/qryAppendixI_NHSF_2015.xlsx');s=pd.read_csv(D/'ALL_sizefreq_2024.csv',low_memory=False);s=s[s.campus.eq('KFM')&s.year.eq(2015)]
checks=[];bins=[]
for site,originals in [('a-k-01',[1,21]),('a-k-14',[37])]:
 for sp in ['Strongylocentrotus franciscanus','Strongylocentrotus purpuratus']:
  q=x[x.Site_Number.isin(originals)&x.Scientific_Name.eq(sp)];g=s[s.site.eq(site)&s.classcode.eq(sp)];expected=q.groupby('Bin_size').apply(lambda z:(z.Percent*z.Cases_N/100).sum(),include_groups=False)
  obs=[]
  for label,v in expected.items():
   nums=list(map(int,re.findall(r'\d+',label)));mask=g['size'].lt(nums[0]) if label.startswith('<') else g['size'].gt(nums[0]) if label.startswith('>') else g['size'].between(nums[0],nums[1]);n=g.loc[mask,'count'].sum();obs.append(n);bins.append(dict(combined_site=site,species=sp,bin=label,combined_n=int(n),nps_n=float(v),difference=float(n-v)))
  error=max(abs(np.asarray(obs)-expected.to_numpy()));npsn=q.groupby('Site_Number').Cases_N.first().sum();assert abs(sum(obs)-g['count'].sum())<1e-8
  checks.append(dict(combined_site=site,species=sp,nps_sites=';'.join(map(str,originals)),combined_n=int(g['count'].sum()),nps_n=int(npsn),max_bin_difference=float(error),bins_match=error<1e-6))
c=pd.DataFrame(checks);assert c.bins_match.all();pd.DataFrame(bins).to_csv(O/'2015_bin_checks.csv',index=False);c.to_csv(O/'2015_site_identity_checks.csv',index=False);print(c.to_string(index=False))
# Capture every retained history using the suspect delivered labels. Post2015 identity
# remains unverified; a label alone does not prove that all later records are misassigned.
m=pd.read_csv('outputs/biomass_support_20261102/private/pool_members.csv',low_memory=False);q=m[m.campus.eq('KFM')&m.site.isin(['a-k-01','a-k-14'])]
rows=[]
for (sid,site),g in q.groupby(['series_id','site']):
 rows.append(dict(series_id=sid,delivered_site=site,first_year=int(g.year.min()),last_year=int(g.year.max()),role_years=g[['year','status']].drop_duplicates().shape[0],donor_record_uses=len(g),animals=int(g['count'].sum()),years=';'.join(map(str,sorted(g.year.unique()))),identity_status='Legacy1985–2011 full-frequency match;2015 grouped-frequency confirmation; other years unresolved',required_action='Reconstruct allowed donor samples from original site-specific records; do not deduplicate or remove by slope'))
pd.DataFrame(rows).to_csv(O/'affected_histories.csv',index=False)
