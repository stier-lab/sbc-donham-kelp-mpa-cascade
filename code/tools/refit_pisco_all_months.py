"""Review-only all-month PISCO swath refits; six flagged animal comparisons.
Dates, equal repeat-visit means, fixed shared zones and source size matching.
No new sites/years, source edits, paper edits or inferred zero biomass.
"""
from pathlib import Path
import pandas as pd,numpy as np,ast,json,hashlib
from scipy.stats import t as student
O=Path('outputs/pisco_all_months_20261102');D=Path('/Users/adrianstier/sbc-donham-kelp-mpa-data/data')
a=pd.read_csv('outputs/biomass_resolution_20261102/review_annual_input.csv');v=pd.read_csv('outputs/all_program_station_refits_20261102/verified_station_values.csv');sids=list(pd.read_csv('outputs/pisco_size_zone_resolution_20261102/summary.csv').series_id);v=v[v.series_id.isin(sids)].copy()
offset=pd.read_csv('outputs/all_program_station_refits_20261102/offsets.csv').set_index('series_id').offset
c=pd.read_csv(D/'PISCO/MLPA_kelpforest_swath_2024.csv',encoding='latin1',low_memory=False);c=c[c.campus.eq('VRG')|(c.campus.eq('UCSB')&~c.zone.isin(['DEEP','SPECIALO']))].copy();c['date']=pd.to_datetime(c[['year','month','day']],errors='coerce');assert c.date.notna().all();c['transect']=c.transect.astype(str)
# Roster from all observed taxa supplies surveyed zeros for absent focal species.
keys=['site','year','date','month','zone','method','transect'];roster=c[keys].drop_duplicates();groupcols=['site','year','date','month','zone','method']
p=pd.read_csv('outputs/biomass_support_20261102/private/pool_members.csv',low_memory=False);p=p[p.series_id.isin(sids)].copy();p['date']=pd.to_datetime(p[['year','month','day']],errors='coerce');p['mass']=np.select([p.Taxa.eq('M. franciscanus'),p.Taxa.eq('S. purpuratus')],[.00059*p['size']**2.917,.00059*p['size']**2.870],default=.001352821*p['size']**2.913963113)
# Local matching must search raw sizes, including donors absent from the old seasonal pool.
rawsizes=pd.read_csv(D/'ALL_sizefreq_2024.csv',low_memory=False);rawsizes['source_row']=np.arange(1,len(rawsizes)+1);rawsizes=rawsizes[rawsizes.campus.isin(['VRG','UCSB']) & rawsizes.classcode.isin(['MESFRA','STRPUR','PANINT'])].copy();rawsizes['Taxa']=rawsizes.classcode.map({'MESFRA':'M. franciscanus','STRPUR':'S. purpuratus','PANINT':'P. interruptus'});rawsizes=rawsizes[rawsizes['count'].gt(0)&rawsizes['size'].gt(0)&(rawsizes.Taxa.eq('P. interruptus')|rawsizes['size'].ge(25))];rawsizes['date']=pd.to_datetime(rawsizes[['year','month','day']],errors='coerce');rawsizes['mass']=np.select([rawsizes.Taxa.eq('M. franciscanus'),rawsizes.Taxa.eq('S. purpuratus')],[.00059*rawsizes['size']**2.917,.00059*rawsizes['size']**2.870],default=.001352821*rawsizes['size']**2.913963113)
lookup={'M. franciscanus':'MESFRAAD','S. purpuratus':'STRPURAD','P. interruptus':'PANINT'}
T=ast.parse(Path('code/tools/audit_station_coverage_refits.py').read_text());exec(compile(ast.Module(body=[n for n in T.body if isinstance(n,ast.FunctionDef) and n.name in ['fit','slope']],type_ignores=[]),'fit','exec'))
siteout=[];visits=[];support=[];zones=[];excluded=[]
scenarios=['season_repeat_oldpool','allmonths_oldpool','allmonths_sharedzones_oldpool','allmonths_sharedzones_visit_sizes','allmonths_sharedzones_exact_sizes']
for sid in sids:
 meta=a[a.series_id==sid].iloc[0];vv=v[v.series_id==sid];tx=meta.Taxa;pg=p[p.series_id==sid];counts=c[c.classcode==lookup[tx]].groupby(keys,as_index=False,dropna=False)['count'].sum(min_count=1).rename(columns={'count':'n'});rr=roster.merge(counts,on=keys,how='left',validate='one_to_one',indicator=True);rr.loc[rr['_merge']=='left_only','n']=0;rr=rr.drop(columns='_merge');rr['density']=rr.n/60
 # Average transects within a dated zone visit. Each visit has equal weight within zone/year.
 ev=rr.groupby(groupcols,as_index=False).agg(density=('density',lambda x:x.mean() if x.notna().all() else np.nan),transects=('transect','nunique'))
 zsets=[set(ev.loc[(ev.site==r.site)&(ev.year==r.year),'zone']) for r in vv.itertuples()];shared=set.intersection(*zsets);assert shared=={'INNER','MID'}
 zones.append(dict(series_id=sid,constant_shared_zones=';'.join(sorted(shared)),station_years=len(vv),lost_station_years=0))
 for r in vv.itertuples():
  z=ev[(ev.site==r.site)&(ev.year==r.year)].copy();don=pg[(pg.year==r.year)&(pg.status==r.status)];oldmass=np.average(don.mass,weights=don['count']) if len(don) else np.nan
  z['matched_mass']=np.nan;z['exact_mass']=np.nan
  for idx,t in z.iterrows():
   dd=rawsizes[(rawsizes.site==r.site)&(rawsizes.year==r.year)&(rawsizes.Taxa==tx)&(rawsizes.zone==t.zone)&rawsizes.date.notna()];gap=np.nan;ns=0;chosen='';kind='no same-site/zone/year sizes'
   if len(dd):
    gaps=(dd.date-t.date).abs().dt.days;gap=int(gaps.min());dd=dd[gaps.eq(gap)];ns=dd['count'].sum();mass=np.average(dd.mass,weights=dd['count']);z.loc[idx,'matched_mass']=mass
    if gap==0:z.loc[idx,'exact_mass']=mass;kind='same visit'
    else:kind='nearest visit in same site/zone/year'
    chosen=';'.join(sorted(set(dd.date.dt.strftime('%Y-%m-%d'))))
   support.append(dict(series_id=sid,site=r.site,status=r.status,year=r.year,zone=t.zone,count_date=str(t.date.date()),count_month=t.month,density=t.density,size_dates=chosen,gap_days=gap,n_sized=ns,match=kind,retained_shared_zone=t.zone in shared))
  z['biomass_oldpool']=z.density*oldmass;z['biomass_visit']=z.density*z.matched_mass;z['biomass_exact']=z.density*z.exact_mass
  for col in ['biomass_oldpool','biomass_visit','biomass_exact']:z.loc[z.density.eq(0),col]=0
  for sc in scenarios:
   zz=z.copy();col='biomass_oldpool'
   if sc=='season_repeat_oldpool':zz=zz[zz.month.between(5,10)]
   if 'sharedzones' in sc:zz=zz[zz.zone.isin(shared)]
   if sc.endswith('visit_sizes'):col='biomass_visit'
   if sc.endswith('exact_sizes'):col='biomass_exact'
   # Any unsized positive visit makes station biomass unknown. Never average it away.
   means=zz.groupby('zone')[['density',col]].mean();bio=means[col].mean() if len(zz) and zz[col].notna().all() else np.nan;den=means.density.mean() if len(zz) and zz.density.notna().all() else np.nan
   siteout.append(dict(series_id=sid,site=r.site,status=r.status,year=r.year,t=r.t,scenario=sc,biomass=bio,density=den,visits=len(zz),months=';'.join(map(str,sorted(set(zz.month)))),zones=';'.join(sorted(set(zz.zone))),missing_positive_visits=int((zz.density.gt(0)&zz[col].isna()).sum())))
   for tr in zz.itertuples():visits.append(dict(series_id=sid,scenario=sc,site=r.site,status=r.status,year=r.year,zone=tr.zone,date=str(tr.date.date()),doy=tr.date.dayofyear,density=tr.density,biomass=getattr(tr,col)))
sv=pd.DataFrame(siteout);su=pd.DataFrame(support);vi=pd.DataFrame(visits)
sv.to_csv(O/'station_years.csv',index=False);su.to_csv(O/'size_matching.csv',index=False);vi.to_csv(O/'visit_values.csv',index=False);pd.DataFrame(zones).to_csv(O/'shared_zones.csv',index=False)
annual=[];effects=[];coverage=[]
# Retain the original reference station composition for every paired year.
for sid in sids:
 for sc in scenarios:
  for resp,col in [('Bio','biomass'),('Den','density')]:
   originalsid=sid if resp=='Bio' else sid.replace('|Bio','|Den');base=a[(a.series_id==originalsid)&a.year.isin(v[v.series_id==sid].year)].copy();xx=sv[(sv.series_id==sid)&(sv.scenario==sc)];vals=xx.groupby(['year','status'])[col].agg(lambda s:s.mean() if s.notna().all() else np.nan).unstack();zz=base.merge(vals,on='year',validate='one_to_one');zz['new']=np.log((zz.mpa+offset[originalsid])/(zz.reference+offset[originalsid]));zz.loc[~np.isfinite(zz.new),'new']=np.nan
   for row in zz.itertuples():annual.append(dict(series_id=originalsid,scenario=sc,year=row.year,t=row.t,original=row.lnDiff,revised=row.new,protected=row.mpa,reference=row.reference))
   post=zz[zz.t>0];valid=post[post.new.notna()];coverage.append(dict(series_id=originalsid,scenario=sc,original_years=len(post),retained_years=len(valid),lost_years=';'.join(map(str,sorted(set(post.year)-set(valid.year))))))
   for label,df,colfit in [('current_all_years',post,'lnDiff'),('current_matched_years',valid,'lnDiff'),('revised_matched_years',valid,'new')]:
    result=slope(df,colfit)
    if result:effects.append(dict(series_id=originalsid,scenario=sc,comparison=label,**result))
pd.DataFrame(annual).to_csv(O/'annual_values.csv',index=False);pd.DataFrame(effects).to_csv(O/'effects.csv',index=False);pd.DataFrame(coverage).to_csv(O/'coverage.csv',index=False)
# Survey-date summaries use the same equal-visit / equal-zone / equal-site hierarchy.
timing=vi.groupby(['series_id','scenario','year','status','site','zone'],as_index=False).doy.mean().groupby(['series_id','scenario','year','status','site'],as_index=False).doy.mean().groupby(['series_id','scenario','year','status'],as_index=False).doy.mean();timing.to_csv(O/'timing_by_role_year.csv',index=False)
tm=[]
for (sid,sc),g in timing.groupby(['series_id','scenario']):
 w=g.pivot(index='year',columns='status',values='doy').dropna();tm.append(dict(series_id=sid,scenario=sc,mean_protected_minus_reference_days=(w.mpa-w.reference).mean(),max_absolute_role_gap_days=(w.mpa-w.reference).abs().max(),protected_doy_trend_days_per_year=np.polyfit(w.index,w.mpa,1)[0],reference_doy_trend_days_per_year=np.polyfit(w.index,w.reference,1)[0],years=len(w)))
pd.DataFrame(tm).to_csv(O/'timing_summary.csv',index=False)
protected=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text());assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in protected.items())
(O/'verification.json').write_text(json.dumps(dict(histories=6,scenarios=len(scenarios),fits=len(effects),station_scenario_rows=len(sv),size_match_events=len(su),protected_unchanged=len(protected),scope='Existing retained station-years; all calendar months; six affected PISCO animal histories and their density responses'),indent=2));print(pd.DataFrame(coverage).query("scenario=='allmonths_sharedzones_visit_sizes'").to_string(index=False));print(pd.DataFrame(effects).query("scenario=='allmonths_sharedzones_visit_sizes' and comparison=='revised_matched_years'")[['series_id','estimate','low','high','n_years']].to_string(index=False))
