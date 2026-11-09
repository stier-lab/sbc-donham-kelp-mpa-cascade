"""Recover dates only from exact site/year/species size-frequency agreement."""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
D=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data');O=Path('outputs/biomass_protocols_20261102');P=O/'private'
sf=pd.read_csv(D/'ALL_sizefreq_2024.csv',low_memory=False);sf['source_row']=np.arange(1,len(sf)+1)
names={11005:'Strongylocentrotus franciscanus',11006:'Strongylocentrotus purpuratus'}
s=sf[sf.campus.eq('KFM') & sf.classcode.isin(names.values())].copy();s['species_code']=s.classcode.map({v:k for k,v in names.items()});s['site_number']=s.site.str.extract(r'a-k-(\d+)')[0].astype(int);s['mm']=s['size'].round(6)
legacy=D/'source_archive/kfm_size_feasibility_2026-10-29/legacy/Invertebrate Size Frequency.csv'
l=pd.read_csv(legacy);l=l[l.Species.isin(names)].copy();assert l.Measurement.eq('Maximum test diameter').all(),l.Measurement.unique()
l['mm']=(l.Size*10).round(6);l['date']=pd.to_datetime(l.Date,format='%d%b%Y');assert (l.date.dt.year==l.Year).all()
l=l.rename(columns={'Site':'site_number','Year':'year','Species':'species_code'})
k=['site_number','year','species_code'];bins=k+['mm']
a=s.groupby(bins)['count'].sum().rename('combined_n');b=l.groupby(bins).size().rename('legacy_n');j=pd.concat([a,b],axis=1).fillna(0);j['equal']=j.combined_n.eq(j.legacy_n)
check=j.reset_index().groupby(k).agg(histogram_equal=('equal','all'),combined_animals=('combined_n','sum'),legacy_animals=('legacy_n','sum'),size_bins=('mm','size')).reset_index();check['scope']=np.where(check.combined_animals.eq(0),'legacy_only',np.where(check.legacy_animals.eq(0),'combined_only','overlap'))
check.to_csv(O/'kfm_histogram_agreement.csv',index=False)
dates=l.groupby(bins).agg(n_dates=('date','nunique'),first_date=('date','min'),last_date=('date','max')).reset_index()
z=s.merge(check,on=k,validate='many_to_one').merge(dates,on=bins,how='left',validate='many_to_one');z['date_recovered']=z.histogram_equal & z.n_dates.eq(1);z['recovered_date']=z.first_date.where(z.date_recovered)
z[['source_row','site','year','classcode','count','size','histogram_equal','n_dates','recovered_date']].to_csv(P/'source_row_date_overlay.csv',index=False)
# Do not split grouped frequency rows among ambiguous dates or deduplicate animals.
summary=z.groupby(['classcode','year']).agg(source_rows=('source_row','size'),animals=('count','sum'),recovered_rows=('date_recovered','sum')).reset_index();summary.to_csv(O/'kfm_date_recovery_by_year.csv',index=False)
m=pd.read_csv('outputs/biomass_support_20261102/private/pool_members.csv',low_memory=False);m=m[m.campus.eq('KFM')];u=m.merge(z[['source_row','date_recovered','recovered_date','histogram_equal','n_dates']],on='source_row',how='left',validate='many_to_one');u['date_recovered']=u.date_recovered.fillna(False)
rows=[]
for sid,g in u.groupby('series_id'):
 rows.append(dict(series_id=sid,kfm_donor_occurrences=len(g),kfm_animals=int(g['count'].sum()),recovered_occurrences=int(g.date_recovered.sum()),recovered_animals=int(g.loc[g.date_recovered,'count'].sum()),missing_occurrences=int((~g.date_recovered).sum()),years_recovered=';'.join(map(str,sorted(g.loc[g.date_recovered,'year'].unique())))))
pd.DataFrame(rows).to_csv(O/'retained_history_date_support.csv',index=False)
# Compare recovered size dates to all observed count dates, without imposing a cutoff.
mbon=list((D/'MBON').glob('*.csv'))
print('counts candidate files:',[x.name for x in mbon if any(v in x.name.lower() for v in ['density','community','benth','kfm'])])
info=dict(combined_kfm_rows=len(s),combined_kfm_animals=int(s['count'].sum()),overlap_site_year_species=int(check.scope.eq('overlap').sum()),exact_histograms=int((check.histogram_equal & check.scope.eq('overlap')).sum()),recovered_source_rows=int(z.date_recovered.sum()),recovered_animals=int(z.loc[z.date_recovered,'count'].sum()),retained_histories=len(rows),retained_recovered_occurrences=int(u.date_recovered.sum()),retained_unique_recovered_records=int(u.loc[u.date_recovered,'source_row'].nunique()),ambiguous_exact_rows=int((z.histogram_equal & z.n_dates.gt(1)).sum()))
(O/'verification.json').write_text(json.dumps(info,indent=2));print(json.dumps(info,indent=2))
inputs=[legacy,D/'ALL_sizefreq_2024.csv',Path('outputs/biomass_support_20261102/private/pool_members.csv'),Path(__file__)]+list(P.glob('*.html'))+list(P.glob('Table*.csv'))
(O/'input_hashes.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},indent=2))
