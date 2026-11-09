"""Outcome-blind coverage/depth audit; no outcome model fitted here."""
from pathlib import Path
import os,json,hashlib
import pandas as pd,numpy as np
O=Path('outputs/depth_reference_policy_20261101');D=Path(os.environ.get('DONHAM_DATA_REPO_PISCO',str(Path.home()/'Donham-Stier-CA-MPA-Data-2026/data/PISCO')))
used=[]
def rd(p,**kw):used.append(str(p));return pd.read_csv(p,low_memory=False,**kw)
master=rd(Path('outputs/cross_site_seasons_20261101/candidate_stations.csv'))
a=rd(Path('outputs/pisco_source_reconciliation_20261101/annual_input.csv'));a=a[(a.scenario=='Supported_corrections')&(a.Source=='PISCO')]
starts=a.groupby('MPA').MPA_Start.first().to_dict()
sw=rd(D/'MLPA_kelpforest_swath_2024.csv',encoding='latin1');sw=sw[sw.site.isin(master.site)&(sw.year<=2023)&sw.zone.isin(['INNER','MID'])].copy()
keys=['site','year','survey_year','month','day','zone','transect'];ev=sw.groupby(keys,as_index=False).depth.agg(lambda x:tuple(sorted(set(float(v) for v in x.dropna() if np.isfinite(v) and v>0))))
ev['depth_n']=ev.depth.map(len);ev['depth_min']=ev.depth.map(lambda x:min(x) if x else np.nan);ev['depth_max']=ev.depth.map(lambda x:max(x) if x else np.nan);ev['depth_center']=ev.depth.map(lambda x:float(np.median(x)) if x else np.nan);ev=ev.drop(columns='depth')
# Only completeness flags, not response magnitudes, enter this audit.
classes=['MACPYRAD','MESFRAAD','STRPURAD','PANINT'];f=sw[sw.classcode.isin(classes)].copy();f['bad']=f['count'].isna()|((f.classcode=='MACPYRAD')&f['size'].isna());bad=f.groupby(keys,as_index=False).bad.any();ev=ev.merge(bad,on=keys,how='left').fillna({'bad':False});assert not ev.duplicated(keys).any()
ev['date']=pd.to_datetime(ev[['year','month','day']]);ev=ev.merge(master,on='site');ev=ev[ev.year>ev.MPA.map(starts)]
ev.to_csv(O/'transect_depth_events.csv',index=False)
ev[ev.year!=ev.survey_year].to_csv(O/'campaign_calendar_differences.csv',index=False)
# All stations, not selected longest pairs. Source coordinates retained separately.
map_=rd(D/'master_site_table_Emilyedit.csv',encoding='utf-8-sig');coords=map_[['site','latitude','longitude']].drop_duplicates();conflicts=coords[coords.site.duplicated(False) & coords.site.isin(master.site)];conflicts.to_csv(O/'coordinate_conflicts.csv',index=False);coord=coords[~coords.site.duplicated(False)].set_index('site')
def distance(p,r):
 if p not in coord.index or r not in coord.index:return np.nan
 lat1,lon1=np.deg2rad(coord.loc[p,['latitude','longitude']].astype(float));lat2,lon2=np.deg2rad(coord.loc[r,['latitude','longitude']].astype(float));h=np.sin((lat2-lat1)/2)**2+np.cos(lat1)*np.cos(lat2)*np.sin((lon2-lon1)/2)**2;return float(6371.0088*2*np.arcsin(np.sqrt(h)))
support=[];station=[]
for (m,z),q in ev.groupby(['MPA','zone']):
 pp=q[(q.status=='mpa')&~q.bad].depth_center.dropna();rr=q[(q.status=='reference')&~q.bad].depth_center.dropna();support.append(dict(MPA=m,zone=z,protected_min=pp.min(),protected_max=pp.max(),reference_min=rr.min(),reference_max=rr.max(),common_min=max(pp.min(),rr.min()),common_max=min(pp.max(),rr.max())))
support=pd.DataFrame(support);support.to_csv(O/'common_depth_support.csv',index=False)
ev=ev.merge(support[['MPA','zone','common_min','common_max']],on=['MPA','zone']);ev['in_support']=ev.depth_center.between(ev.common_min,ev.common_max)
for (m,site,status,zone),q in ev.groupby(['MPA','site','status','zone']):
 yearly=q.groupby('year').depth_center.median().dropna();ys=yearly.index.to_numpy();slope=np.polyfit(ys-ys.min(),yearly,1)[0] if len(ys)>=3 else np.nan
 station.append(dict(MPA=m,site=site,status=status,zone=zone,n_events=len(q),known=int(q.depth_center.notna().sum()),multivalued=int((q.depth_n>1).sum()),median_depth=q.depth_center.median(),depth_min=q.depth_min.min(),depth_max=q.depth_max.max(),annual_depth_slope=slope,n_years=len(yearly)))
pd.DataFrame(station).to_csv(O/'station_depth_summary.csv',index=False)
visits=[];matches=[];pairrows=[];summ=[];depthpair=[]
for depth_rule in ['Common_depth','Measured_shallow','Measured_15','Measured_20','Zone_only']:
 q=ev[~ev.bad & (ev.in_support if depth_rule=='Common_depth' else (ev.depth_center.notna() & (ev.depth_center<=({'Measured_15':15,'Measured_20':20}.get(depth_rule,16.25)))) if depth_rule.startswith('Measured') else True)].copy()
 v=q.groupby(['MPA','site','status','year','survey_year','zone','date'],as_index=False).agg(n_transects=('transect','nunique'),depth_center=('depth_center','median'),depth_min=('depth_min','min'),depth_max=('depth_max','max'),missing_depth=('depth_center',lambda x:int(x.isna().sum())))
 v['depth_rule']=depth_rule;visits.append(v);v=v[v.n_transects>=2].copy()
 for m,mm in master.groupby('MPA'):
  vp=v[(v.MPA==m)&(v.status=='mpa')];vr=v[(v.MPA==m)&(v.status=='reference')]
  pairs=vp.merge(vr,on=['MPA','year','survey_year','zone'],suffixes=('_p','_r'))
  pairs['gap_days']=(pairs.date_p-pairs.date_r).abs().dt.days
  # One closest visit pairing per protected/reference pair, year and stratum.
  pairs=pairs.sort_values(['gap_days','date_p','date_r','site_p','site_r']).drop_duplicates(['site_p','site_r','year','survey_year','zone'])
  pairs['depth_rule']=depth_rule;matches.append(pairs)
  for gap in [0,30,60,90]:
   good=pairs if gap==0 else pairs[pairs.gap_days<=gap]
   both=good.groupby(['site_p','site_r','year','survey_year']).filter(lambda x:set(x.zone)=={'INNER','MID'})
   # Keep both strata in one campaign. If two campaigns occur in a calendar year,
   # prefer the tighter joint date match, with campaign ID as deterministic tie breaker.
   ranks=both.groupby(['site_p','site_r','year','survey_year'],as_index=False).agg(max_gap=('gap_days','max'),sum_gap=('gap_days','sum'))
   ranks=ranks.sort_values(['max_gap','sum_gap','survey_year']).drop_duplicates(['site_p','site_r','year'])
   both=both.merge(ranks[['site_p','site_r','year','survey_year']],on=['site_p','site_r','year','survey_year'])
   for minyears in [3,5]:
    minspan=2 if minyears==3 else 5;sc=f'{depth_rule}_d{gap}_n{minyears}'
    selected=[]
    for ps in sorted(mm[mm.status=='mpa'].site):
     for rs in sorted(mm[mm.status=='reference'].site):
      z=both[(both.site_p==ps)&(both.site_r==rs)];ys=sorted(z.year.unique());span=ys[-1]-ys[0] if ys else 0;ok=len(ys)>=minyears and span>=minspan
      # The recommended measured-depth analysis requires overlapping observed
      # ranges in each stratum over the paired record; no per-year exact-depth cut.
      overlap=True
      for zzname in ['INNER','MID']:
       dz=z[z.zone==zzname];pp=dz.depth_center_p.dropna();rr=dz.depth_center_r.dropna()
       if not len(pp) or not len(rr) or max(pp.min(),rr.min())>min(pp.max(),rr.max()):overlap=False
      if depth_rule.startswith('Measured'):ok=ok and overlap
      pairrows.append(dict(MPA=m,scenario=sc,depth_rule=depth_rule,gap=gap,min_years=minyears,protected=ps,reference=rs,n_years=len(ys),span=span,years=';'.join(map(str,ys)),eligible=ok,paired_depth_support_overlap=overlap,distance_km=distance(ps,rs)))
      if ok:
       z=z.copy();z['scenario']=sc;selected.append(z)
    zz=pd.concat(selected) if selected else both.iloc[:0].copy();years=sorted(zz.year.unique());oldyears=sorted(a[(a.MPA==m)&(a.Taxa=='M. pyrifera')&(a.Resp=='Bio')&(a.t>0)].year)
    summ.append(dict(MPA=m,scenario=sc,n_pairs=len(selected),n_protected=zz.site_p.nunique(),n_reference=zz.site_r.nunique(),protected=';'.join(sorted(zz.site_p.unique())),reference=';'.join(sorted(zz.site_r.unique())),n_years=len(years),years=';'.join(map(str,years)),original_years=';'.join(map(str,oldyears)),lost_years=';'.join(map(str,sorted(set(oldyears)-set(years)))),added_years=';'.join(map(str,sorted(set(years)-set(oldyears))))))
    if selected:depthpair.append(zz)
pd.concat(visits).to_csv(O/'visit_coverage.csv',index=False);pd.concat(matches).to_csv(O/'closest_candidate_visits.csv',index=False);pd.DataFrame(pairrows).to_csv(O/'pair_decisions.csv',index=False);pd.DataFrame(summ).to_csv(O/'mpa_decisions.csv',index=False);pd.concat(depthpair).to_csv(O/'selected_pair_zone_years.csv',index=False)
# Depth differences are paired by year/campaign/zone, not a comparison of unrelated global means.
p=pd.concat(matches);p['depth_difference']=p.depth_center_p-p.depth_center_r
p.groupby(['MPA','depth_rule','site_p','site_r','zone']).agg(n_years=('year','nunique'),median_difference=('depth_difference','median'),min_difference=('depth_difference','min'),max_difference=('depth_difference','max'),median_gap=('gap_days','median')).reset_index().to_csv(O/'pair_depth_differences.csv',index=False)
(O/'source_hashes.json').write_text(json.dumps({p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in used},indent=2)+'\n')
print(pd.DataFrame(summ).query("scenario in ['Measured_shallow_d90_n5','Measured_shallow_d60_n5','Zone_only_d90_n5']").to_string(index=False));print(support.to_string(index=False))
