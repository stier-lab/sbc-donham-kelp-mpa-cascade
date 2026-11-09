"""Apply frozen, result-independent rules and expose every comparison/year."""
from common import *
u=pd.read_parquet(PRIVATE/'measurement_units.parquet');u['row_id']=u.index
c=read_tracked(ROOT/'config/integrated_comparison_scope.csv')
# The existing station lists are starting assignments, not proof that exclusions applied.
ltermap=read_tracked(DATA/'config/lter_physical_sites.csv').set_index('site').original_site.to_dict()
starts=read_source('site_assignments')[['CA_MPA_Name_Short','MPA_Start']].drop_duplicates();assert starts.CA_MPA_Name_Short.is_unique
startmap=starts.set_index('CA_MPA_Name_Short').MPA_Start.to_dict()
# Satellite metadata use abbreviated names; use their recorded start year, not a name guess.
a=read_tracked(ROOT/'outputs/pisco_fish_rules_20261101/approved_annual_input.csv');ss=a[['MPA','MPA_Start']].drop_duplicates();assert ss.MPA.is_unique
for r in ss.itertuples():startmap.setdefault(r.MPA,r.MPA_Start)
assign=[];uses=[];years=[]
for r in c.itertuples():
 if r.exclude:
  assign.append(dict(series_id=r.series_id,source=r.Source,site='',role='',included=False,reason=r.reason));continue
 q=u[(u.Source==r.Source)&(u.Taxa==r.Taxa)]
 stations=[]
 for role,col in [('mpa','mapped_protected_sites'),('reference','mapped_reference_sites')]:
  if r.Source=='Landsat':sites=sorted(s for s in q.site.unique() if s.startswith(r.MPA+'|'+role+'|'))
  else:sites=sorted({ltermap.get(s,s) if r.Source=='LTER' else s for s in str(getattr(r,col)).split(';') if s!='nan'})
  for site in sites:
   omit=blocked(r.Source,site,r.Taxa)
   assign.append(dict(series_id=r.series_id,source=r.Source,site=site,role=role,included=not omit,reason='Response-specific researcher exclusion' if omit else 'Documented assignment; annual availability evaluated separately'))
   if not omit:stations.append((site,role))
 if r.MPA not in startmap:raise ValueError('No verified protection year: '+r.MPA)
 for yr,g in q[q.site.isin([s for s,_ in stations])].groupby('year'):
  z=g.merge(pd.DataFrame(stations,columns=['site','status']),on='site',validate='many_to_one')
  sitezones=z.groupby('site').zone.agg(set);shared=set.intersection(*sitezones) if len(sitezones) else set()
  roles=set(z.status);keep=z[z.zone.isin(shared)].copy();ok=roles=={'mpa','reference'} and len(keep)>0
  years.append(dict(series_id=r.series_id,MPA=r.MPA,Source=r.Source,Taxa=r.Taxa,Resp=r.Resp,year=yr,MPA_Start=int(startmap[r.MPA]),t=yr-int(startmap[r.MPA]),survey_eligible=ok,reason='shared zones and both roles' if ok else 'no common zones or missing role',stations=';'.join(sorted(z.site.unique())),protected_stations=';'.join(sorted(z.loc[z.status=='mpa','site'].unique())),reference_stations=';'.join(sorted(z.loc[z.status=='reference','site'].unique())),shared_zones=';'.join(sorted(shared)),visits=len(keep[['site','date','zone']].drop_duplicates()),depth_equivalence='program-zone matching only' if r.Source=='PISCO' else 'physical depth unavailable; established stations retained',first_survey_date=keep.date.min() if len(keep) else '',last_survey_date=keep.date.max() if len(keep) else ''))
  if ok:
   keep['series_id']=r.series_id;keep['MPA_Start']=int(startmap[r.MPA]);keep['Resp']=r.Resp;uses.append(keep)
pd.DataFrame(assign).to_csv(OUT/'station_assignments.csv',index=False);pd.DataFrame(years).to_csv(OUT/'comparison_year_roster.csv',index=False)
x=pd.concat(uses,ignore_index=True);x.to_parquet(PRIVATE/'measurement_uses.parquet',index=False)
# Coverage sensitivities never change the main roster in response to a result.
sens=[]
for sid,g in x.groupby('series_id'):
 for variant in ['available','constant_stations','constant_zones']:
  z=g.copy();all_years=set(g.year)
  if variant=='constant_stations':
   keep=g.groupby('site').year.agg(set);keep=set(keep[keep.map(lambda v:v==all_years)].index);z=z[z.site.isin(keep)]
  if variant=='constant_zones':
   sets=g.groupby('year').zone.agg(set);shared=set.intersection(*sets);z=z[z.zone.isin(shared)]
  for yr,h in g.groupby('year'):
   k=z[z.year.eq(yr)];sens.append(dict(series_id=sid,variant=variant,year=yr,eligible=set(k.status)=={'mpa','reference'},stations=k.site.nunique(),zones=';'.join(sorted(k.zone.unique())),units=len(k),lost_stations=';'.join(sorted(set(h.site)-set(k.site)))))
pd.DataFrame(sens).to_csv(OUT/'coverage_sensitivity.csv',index=False);write_inputs();print('Comparison/year rows',len(years),'measurement uses',len(x))
