"""Source-record count/size compatibility checks; no new cutoffs or primary edits."""
from pathlib import Path
import pandas as pd,numpy as np
from biomass_source_guard import source_issue
O=Path('outputs/biomass_joint_checks_20261102');D=Path('/Users/adrianstier/sbc-donham-kelp-mpa-data/data')
p=pd.read_csv('outputs/biomass_support_20261102/private/pool_members.csv',low_memory=False)
v=pd.read_csv('outputs/all_program_station_refits_20261102/verified_station_values.csv');v=v[v.Resp.eq('Bio')]
x=pd.read_csv(D/'PISCO/MLPA_kelpforest_swath_2024.csv',encoding='latin1',low_memory=False)
x=x[x.month.between(5,10)&(x.campus.eq('VRG')|(x.campus.eq('UCSB')&~x.zone.isin(['DEEP','SPECIALO'])))].copy()
x['date']=pd.to_datetime(x[['year','month','day']],errors='coerce');x['Source']='PISCO'
# Use the survey roster across taxa, including surveys on which a focal taxon was absent.
e=x[['Source','site','year','date','zone','depth']].drop_duplicates()
y=pd.read_csv(D/'MBON/SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv',low_memory=False)
y=y[(y.data_source.eq('lter')&y.sample_method.eq('quad'))|y.data_source.eq('kfm')].copy();y['date']=pd.to_datetime(y.date,errors='coerce');y['year']=y.date.dt.year;y['Source']=y.data_source.str.upper();y=y.rename(columns={'site_id':'site'});y['zone']=np.nan;y['depth']=np.nan
e=pd.concat([e,y[e.columns]],ignore_index=True).drop_duplicates()
over=pd.read_csv('outputs/biomass_protocols_20261102/private/source_row_date_overlay.csv');over=over[['source_row','recovered_date']].drop_duplicates();assert not over.source_row.duplicated().any()
p=p.merge(over,on='source_row',how='left',validate='many_to_one');p['date']=pd.to_datetime(p[['year','month','day']],errors='coerce');p['date']=p.date.fillna(pd.to_datetime(p.recovered_date,errors='coerce'))
events={k:g for k,g in e.groupby(['Source','site','year'])}
rows=[]
for (sid,yr,st),g in p.groupby(['series_id','year','status']):
 src=g.Source.iloc[0];tx=g.Taxa.iloc[0];sites=set(v.loc[(v.series_id==sid)&(v.year==yr)&(v.status==st),'site'])
 ev=pd.concat([events.get((src,s,yr),e.iloc[:0]) for s in sites],ignore_index=True)
 # Aggregate repeated frequencies with identical matching attributes; preserve multiplicity.
 attrs=['site','date','zone','size','depth'];gg=g.groupby(attrs,dropna=False,as_index=False)['count'].sum()
 total=g['count'].sum();same=known=zero=zone_known=zone_match=bad=depth_known=depth_within=0.;gaps=[];weights=[]
 for r in gg.itertuples():
  issue=source_issue(r.site,tx);bad+=r.count*bool(issue)
  physical=r.site in sites and not issue;same+=r.count*physical
  ee=ev[ev.site.eq(r.site)] if physical else ev.iloc[:0]
  if pd.notna(r.date) and len(ee) and ee.date.notna().any():
   gap=(ee.date-r.date).abs().dt.days.min();known+=r.count;zero+=r.count*(gap==0);gaps.append(gap);weights.append(r.count)
  depths=pd.to_numeric(ee.depth,errors='coerce').dropna();sd=pd.to_numeric(r.depth,errors='coerce')
  if pd.notna(sd) and len(depths):
   depth_known+=r.count;depth_within+=r.count*(depths.min()<=sd<=depths.max())
  zz=set(ee.zone.dropna())
  if pd.notna(r.zone) and zz:
   zone_known+=r.count;zone_match+=r.count*(r.zone in zz)
 rows.append(dict(series_id=sid,Source=src,Taxa=tx,year=yr,status=st,n_animals=total,fraction_verified_label_site=same/total,fraction_same_site_date_evaluable=known/total,fraction_same_site_exact_date=zero/total,max_nearest_date_gap_days=max(gaps) if gaps else np.nan,mean_nearest_date_gap_days=np.average(gaps,weights=weights) if gaps else np.nan,fraction_zone_evaluable=zone_known/total,fraction_matching_recorded_zone=zone_match/total,fraction_blocked_source=bad/total,fraction_numeric_depth_evaluable=depth_known/total,fraction_within_count_depth_range=depth_within/total,min_size_mm=g['size'].min(),max_size_mm=g['size'].max(),fraction_below_25mm=np.average(g['size']<25,weights=g['count']),interpretation='Label-site match excludes known aliases; remaining geography and habitat representativeness not certified. Date gaps conditional on evaluable same-site donors.'))
r=pd.DataFrame(rows);r.to_csv(O/'count_size_matching.csv',index=False)
s=[]
for sid,g in r.groupby('series_id'):
 s.append(dict(series_id=sid,Source=g.Source.iloc[0],Taxa=g.Taxa.iloc[0],role_years=len(g),role_years_with_borrowed_or_unresolved_sizes=int(g.fraction_verified_label_site.lt(1-1e-10).sum()),role_years_without_complete_date_match=int(g.fraction_same_site_exact_date.lt(1-1e-10).sum()),role_years_with_source_exclusion=int(g.fraction_blocked_source.gt(0).sum()),max_nearest_date_gap_days=g.max_nearest_date_gap_days.max(),zone_evaluable_role_years=int(g.fraction_zone_evaluable.gt(0).sum()),action='Withhold joint resampling pending source reconciliation' if g.fraction_blocked_source.gt(0).any() else 'Conditional resampling allowed; representativeness remains separate'))
pd.DataFrame(s).to_csv(O/'count_size_history_summary.csv',index=False)
print(pd.DataFrame(s).groupby('Source').sum(numeric_only=True).to_string())
