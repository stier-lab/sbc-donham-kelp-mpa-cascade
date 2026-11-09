"""Trace 23 size-zone flags to raw records; isolated matched-year biomass sensitivities."""
from pathlib import Path
import pandas as pd,numpy as np,ast,json,hashlib
from scipy.stats import t as student
O=Path('outputs/pisco_size_zone_resolution_20261102');D=Path('/Users/adrianstier/sbc-donham-kelp-mpa-data/data');J=Path('outputs/biomass_joint_checks_20261102')
f=pd.read_csv(J/'size_zone_mismatches.csv');sids=set(f.series_id)
p=pd.read_csv('outputs/biomass_support_20261102/private/pool_members.csv',low_memory=False);p=p[p.series_id.isin(sids)].copy()
original_sizes=pd.read_csv(D/'ALL_sizefreq_2024.csv',low_memory=False)
unique=p.drop_duplicates('source_row').reset_index(drop=True);original=original_sizes.iloc[unique.source_row.to_numpy()-1].reset_index(drop=True)
for field in ['site','zone']:
 assert unique[field].fillna('__NA__').astype(str).equals(original[field].fillna('__NA__').astype(str)),field
for field in ['year','month','day','count','size']:
 assert np.allclose(pd.to_numeric(unique[field],errors='coerce'),pd.to_numeric(original[field],errors='coerce'),equal_nan=True),field
v=pd.read_csv('outputs/all_program_station_refits_20261102/verified_station_values.csv');v=v[v.series_id.isin(sids)]
c=pd.read_csv(D/'PISCO/MLPA_kelpforest_swath_2024.csv',encoding='latin1',low_memory=False);c['date']=pd.to_datetime(c[['year','month','day']],errors='coerce')
cf=c[c.month.between(5,10)&(c.campus.eq('VRG')|(c.campus.eq('UCSB')&~c.zone.isin(['DEEP','SPECIALO'])))].copy()
fish=pd.read_csv(D/'PISCO/MLPA_kelpforest_fish_2024.csv',encoding='latin1',low_memory=False)
# These are roster records, not just records of focal species: absence is not missing effort.
raw={(s,y):g for (s,y),g in c.groupby(['site','year'])};filtered={(s,y):g for (s,y),g in cf.groupby(['site','year'])}
fishgroups={(s,y):g for (s,y),g in fish.groupby(['site','year'])}
rows=[];pool=[];flags=[]
for (sid,yr,st),g in p.groupby(['series_id','year','status']):
 sites=set(v.loc[(v.series_id==sid)&(v.year==yr)&(v.status==st),'site']);g=g.copy();g['same_site']=g.site.isin(sites);g['zone_match']=False;g['known_zone_mismatch']=False
 for (site,zone),q in g.groupby(['site','zone'],dropna=False):
  if site not in sites:continue
  rr=raw.get((site,yr),c.iloc[:0]);ff=filtered.get((site,yr),c.iloc[:0]);zones=set(ff.zone.dropna());known=pd.notna(zone) and bool(zones)
  if known and zone in zones:g.loc[q.index,'zone_match']=True
  elif known:
   g.loc[q.index,'known_zone_mismatch']=True
   zz=rr[rr.zone.eq(zone)];count_dates=sorted(set(ff.date.dropna().dt.strftime('%Y-%m-%d')));size_dates=pd.to_datetime(q[['year','month','day']],errors='coerce');sd=sorted(set(size_dates.dropna().dt.strftime('%Y-%m-%d')))
   fq=fishgroups.get((site,yr),fish.iloc[:0]);reason='Count zone absent from delivered swath file; unsurveyed versus missing records unresolved'
   if len(zz):reason='Zone present in raw swath but removed by retained month/program filters'
   rows.append(dict(series_id=sid,year=yr,status=st,site=site,size_zone=zone,n_animals=q['count'].sum(),retained_count_zones=';'.join(sorted(zones)),all_raw_count_zones=';'.join(sorted(set(rr.zone.dropna()))),size_dates=';'.join(sd),retained_count_dates=';'.join(count_dates),raw_zone_dates=';'.join(sorted(set(zz.date.dropna().dt.strftime('%Y-%m-%d')))),size_depths=';'.join(map(str,sorted(set(pd.to_numeric(q.depth,errors='coerce').dropna())))),count_depth_min=pd.to_numeric(ff.depth,errors='coerce').min(),count_depth_max=pd.to_numeric(ff.depth,errors='coerce').max(),fish_has_zone=bool(fq.zone.eq(zone).any()),all_size_dates_have_raw_same_zone_count=bool(sd) and set(sd).issubset(set(zz.date.dropna().dt.strftime('%Y-%m-%d'))),raw_zone_months=';'.join(map(str,sorted(set(zz.month)))),reason=reason))
 g['mass']=np.select([g.Taxa.eq('M. franciscanus'),g.Taxa.eq('S. purpuratus')],[.00059*g['size']**2.917,.00059*g['size']**2.870],default=.001352821*g['size']**2.913963113)
 base=np.average(g.mass,weights=g['count'])
 for scenario,mask in [('current',np.ones(len(g),bool)),('remove_known_zone_mismatch',~g.known_zone_mismatch),('same_site_only',g.same_site),('same_site_and_zone',g.zone_match),('same_season',g.month.between(5,10)),('same_site_zone_season',g.zone_match & g.month.between(5,10))]:
  q=g[mask];new=np.average(q.mass,weights=q['count']) if len(q) else np.nan
  pool.append(dict(series_id=sid,year=yr,status=st,scenario=scenario,n_original=g['count'].sum(),n_retained=q['count'].sum(),original_mean_mass=base,new_mean_mass=new,mass_ratio=new/base,known_zone_mismatch_animals=g.loc[g.known_zone_mismatch,'count'].sum()))
details=pd.DataFrame(rows);details.to_csv(O/'source_zone_details.csv',index=False)
check=details.groupby(['series_id','year','status']).n_animals.sum().reset_index();check=check.merge(f,on=['series_id','year','status'],validate='one_to_one');assert len(check)==23;assert np.allclose(check.n_animals_x,(check.fraction_zone_evaluable-check.fraction_matching_recorded_zone)*check.n_animals_y)
po=pd.DataFrame(pool);po.to_csv(O/'pool_changes.csv',index=False)
T=ast.parse(Path('code/tools/audit_station_coverage_refits.py').read_text());exec(compile(ast.Module(body=[n for n in T.body if isinstance(n,ast.FunctionDef) and n.name in ['fit','slope']],type_ignores=[]),'fit','exec'))
a=pd.read_csv('outputs/biomass_resolution_20261102/review_annual_input.csv');inv=pd.read_csv('outputs/biomass_support_20261102/animal_pool_inventory.csv').set_index(['series_id','year','status']);offset=pd.read_csv('outputs/all_program_station_refits_20261102/offsets.csv').set_index('series_id').offset
out=[];annual=[];loss=[]
for sid in sorted(sids):
 z=a[(a.series_id==sid)&(a.t>0)].sort_values('year');base=slope(z,'lnDiff');out.append(dict(series_id=sid,scenario='current_all_years',**base))
 for sc in ['remove_known_zone_mismatch','same_site_only','same_site_and_zone','same_season','same_site_zone_season']:
  lookup=po[(po.series_id==sid)&(po.scenario==sc)].set_index(['year','status']);zz=z.copy();vals={}
  for st in ['mpa','reference']:
   val=np.array([inv.loc[(sid,yr,st),'value'] for yr in zz.year]);ratio=np.array([lookup.loc[(yr,st),'mass_ratio'] if (yr,st) in lookup.index else np.nan for yr in zz.year]);ratio[val==0]=1;vals[st]=val*ratio
  zz['new']=np.log((vals['mpa']+offset[sid])/(vals['reference']+offset[sid]));valid=zz[np.isfinite(zz.new)].copy();loss.append(dict(series_id=sid,scenario=sc,original_years=len(zz),retained_years=len(valid),lost_years=';'.join(map(str,sorted(set(zz.year)-set(valid.year))))))
  for r in valid.itertuples():annual.append(dict(series_id=sid,scenario=sc,year=r.year,t=r.t,original=r.lnDiff,revised=r.new))
  for label,col in [('current_matched_years','lnDiff'),('revised_matched_years','new')]:
   fitout=slope(valid,col)
   if fitout:out.append(dict(series_id=sid,scenario=sc+'__'+label,**fitout))
pd.DataFrame(out).to_csv(O/'effects.csv',index=False);pd.DataFrame(annual).to_csv(O/'annual_sensitivities.csv',index=False);pd.DataFrame(loss).to_csv(O/'coverage_changes.csv',index=False)
protected=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text());assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in protected.items())
inputs=[D/'ALL_sizefreq_2024.csv',D/'PISCO/MLPA_kelpforest_swath_2024.csv',D/'PISCO/MLPA_kelpforest_fish_2024.csv',Path('outputs/biomass_support_20261102/private/pool_members.csv'),Path('outputs/biomass_resolution_20261102/review_annual_input.csv')]
(O/'input_hashes.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},indent=2));(O/'verification.json').write_text(json.dumps(dict(original_size_rows_verified=len(unique),flagged_groups_reproduced=len(check),site_zone_groups=len(details),histories=len(sids),fits=len(out),protected_unchanged=len(protected)),indent=2))
print(details.groupby('reason').agg(zone_groups=('site','size'),animals=('n_animals','sum')).to_string());print(pd.DataFrame(loss).to_string(index=False));print(pd.DataFrame(out)[['series_id','scenario','estimate','low','high']].to_string(index=False))
