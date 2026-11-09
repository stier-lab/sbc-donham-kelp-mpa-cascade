from pathlib import Path
import pandas as pd,numpy as np,json
O=Path('outputs/biomass_protocols_20261102');D=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data')
z=pd.read_csv(O/'private/source_row_date_overlay.csv');z=z[z.recovered_date.notna()].copy();z['date']=pd.to_datetime(z.recovered_date)
c=pd.read_csv(D/'MBON/SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv',low_memory=False);c=c[c.data_source.eq('kfm') & c.taxon_name.isin(['Mesocentrotus franciscanus','Strongylocentrotus purpuratus'])].copy();c['year']=pd.to_datetime(c.date).dt.year;c['classcode']=c.taxon_name.replace({'Mesocentrotus franciscanus':'Strongylocentrotus franciscanus'})
lookup={k:pd.to_datetime(g.date.unique()) for k,g in c.groupby(['site_id','year','classcode'])}
rows=[]
for k,g in z.groupby(['site','year','classcode']):
 dates=lookup.get(k);sd=g.date.unique();gaps=[] if dates is None else [min(abs((pd.Timestamp(d)-x).days) for x in dates) for d in sd]
 rows.append(dict(site=k[0],year=k[1],species=k[2],size_dates=';'.join(str(pd.Timestamp(d).date()) for d in sorted(sd)),count_dates='' if dates is None else ';'.join(str(x.date()) for x in sorted(dates)),n_size_dates=len(sd),max_nearest_gap_days=max(gaps) if gaps else np.nan,all_size_dates_match_counts=bool(gaps) and max(gaps)==0,n_size_rows=len(g)))
r=pd.DataFrame(rows);r.to_csv(O/'kfm_count_size_date_gaps.csv',index=False)
# Only site-years actually supplying retained donor pools are relevant to the review cohort.
m=pd.read_csv('outputs/biomass_support_20261102/private/pool_members.csv',low_memory=False);m=m[m.campus.eq('KFM')];rr=m[['series_id','site','year','classcode']].drop_duplicates().merge(r,left_on=['site','year','classcode'],right_on=['site','year','species'],how='inner',validate='many_to_one');rr.to_csv(O/'retained_count_size_date_gaps.csv',index=False)
# Show discrepancies rather than choosing a source or altering sizes.
check=pd.read_csv(O/'kfm_histogram_agreement.csv');bad=check[check.scope.eq('overlap') & ~check.histogram_equal].copy();bad['site']=bad.site_number.map(lambda n:f'a-k-{n:02d}');bad['classcode']=bad.species_code.map({11005:'Strongylocentrotus franciscanus',11006:'Strongylocentrotus purpuratus'})
bad=bad.merge(m.groupby(['site','year','classcode']).series_id.agg(lambda x:';'.join(sorted(set(x)))).rename('affected_histories').reset_index(),on=['site','year','classcode'],how='left');bad.to_csv(O/'kfm_source_discrepancies.csv',index=False)
info=json.loads((O/'verification.json').read_text());info.update(recovered_site_year_species=len(r),date_groups_with_count_records=int(r.count_dates.ne('').sum()),all_dates_match_groups=int(r.all_size_dates_match_counts.sum()),max_date_gap_days=float(r.max_nearest_gap_days.max()),retained_date_group_uses=len(rr),retained_nonzero_date_gaps=int(rr.max_nearest_gap_days.gt(0).sum()),source_discrepancies_affecting_retained=int(bad.affected_histories.notna().sum()));(O/'verification.json').write_text(json.dumps(info,indent=2));print(json.dumps(info,indent=2));print(rr[rr.max_nearest_gap_days.gt(0)][['site','year','species','max_nearest_gap_days']].drop_duplicates().head(12).to_string(index=False))
