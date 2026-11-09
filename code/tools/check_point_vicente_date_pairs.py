from pathlib import Path
import pandas as pd,numpy as np
O=Path('outputs/point_vicente_seasons_20261101');c=pd.read_csv(O/'survey_dates_depths.csv');grid=pd.read_csv(O/'coverage_grid.csv');rows=[];rules=[]
for season,lo,hi in [('AprOct',4,10),('AprNov',4,11)]:
 for zn,zones in [('INNER_MID',['INNER','MID']),('INNER_MID_OUTER',['INNER','MID','OUTER'])]:
  r=grid[(grid.season==season)&(grid.zones==zn)&(grid.year_window=='full11')].iloc[0];refs=r.reference.split(';');q=c[c.month.between(lo,hi)&c.zone.isin(zones)].copy();q['date']=pd.to_datetime(q[['year','month','day']])
  for site in refs:
   for year in range(2013,2024):
    distances=[]
    for zone in zones:
     p=q[(q.site=='POINT_VICENTE_W')&(q.year==year)&(q.zone==zone)].date.to_numpy();r_=q[(q.site==site)&(q.year==year)&(q.zone==zone)].date.to_numpy();assert len(p) and len(r_);distances.extend(np.abs((p[:,None]-r_[None,:])/np.timedelta64(1,'D')).ravel())
    rows.append(dict(season=season,zones=zn,reference=site,year=year,min_gap_days=min(distances),max_gap_days=max(distances)))
d=pd.DataFrame(rows);d.to_csv(O/'station_pair_date_gaps.csv',index=False)
for (s,z),q in d.groupby(['season','zones']):
 for limit in [30,60,90]:
  refs=q.groupby('reference').max_gap_days.max();ok=refs[refs<=limit];rules.append(dict(season=s,zones=z,max_gap_rule_days=limit,n_references=len(ok),references=';'.join(ok.index)))
pd.DataFrame(rules).to_csv(O/'date_rule_grid.csv',index=False);print(pd.DataFrame(rules).to_string(index=False))
