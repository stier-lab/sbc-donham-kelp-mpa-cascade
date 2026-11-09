"""Coverage diagnostics only; preserve method identities and outcome-independent choices."""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
O=Path('outputs/coverage_decisions_20261101');I=Path('outputs/dataset_inventory_20261101')
def read(n):return pd.read_csv(I/(n+'.csv'),low_memory=False,dtype={'site_id':str})
def parts(v):return set() if pd.isna(v) or str(v) in ['Not supplied','nan',''] else set(str(v).split(';'))
r=read('raw_station_year_audit');s=read('site_response_coverage');p=read('candidate_pair_coverage');c=read('comparison_eligibility')
ret=set(map(tuple,c.loc[~c.exclude,['Source','MPA','Taxa']].to_numpy()))
lookup={(a.dataset_id,str(a.site_id),a.taxon,a.method,int(a.year)):a for a in r.itertuples()}
site=s.set_index('site_record_id');rows=[]
for a in p.itertuples():
 if (a.program,a.MPA,a.taxon) not in ret:continue
 ps=site.loc[a.protected_record_id];rs=site.loc[a.reference_record_id]
 if any(bool(x) for x in [ps.inherited_pisco_reference_exclusion,rs.inherited_pisco_reference_exclusion,ps.inherited_whole_mpa_exclusion,rs.inherited_whole_mpa_exclusion]):continue
 if a.program=='KFM' and a.taxon!='S. pulcher' and (ps.inherited_kfm_nonfish_exclusion or rs.inherited_kfm_nonfish_exclusion):continue
 for y in sorted(int(z) for z in parts(a.observed_years)):
  x=lookup[(a.dataset_id,str(a.protected_site),a.taxon,a.method,y)];z=lookup[(a.dataset_id,str(a.reference_site),a.taxon,a.method,y)]
  zx,zz=parts(x.zones),parts(z.zones);dx=pd.to_datetime(list(parts(x.dates)),errors='coerce').dropna();dz=pd.to_datetime(list(parts(z.dates)),errors='coerce').dropna()
  gap=min(abs((aa-bb).days) for aa in dx for bb in dz) if len(dx) and len(dz) else np.nan
  rows.append(dict(pair_id=a.pair_id,program=a.program,MPA=a.MPA,dataset_id=a.dataset_id,taxon=a.taxon,method=a.method,protected_site=a.protected_site,reference_site=a.reference_site,year=y,shared_zones=';'.join(sorted(zx&zz)),zones_known=bool(zx and zz),protected_months=x.months,reference_months=z.months,shared_month=bool(parts(x.months)&parts(z.months)),nearest_date_gap_days=gap,measured_depth_both=x.positive_depth_records>0 and z.positive_depth_records>0))
d=pd.DataFrame(rows);d.to_csv(O/'pair_year_timing_zones.csv',index=False)
out=[];zones=[]
for k,g in d.groupby('pair_id'):
 sets=[parts(v) for v in g.shared_zones];common=set.intersection(*sets);union=set.union(*sets)
 out.append(dict(pair_id=k,program=g.program.iloc[0],MPA=g.MPA.iloc[0],taxon=g.taxon.iloc[0],method=g.method.iloc[0],paired_years=len(g),constant_zones_all_years=';'.join(sorted(common)),any_shared_zone_years=sum(bool(v) for v in sets),unknown_zone_years=int((~g.zones_known).sum()),years_with_measured_depth_both=int(g.measured_depth_both.sum()),years_no_shared_month=int((~g.shared_month).sum()),years_gap_over90=int(g.nearest_date_gap_days.gt(90).sum())))
 for zone in sorted(union):
  keep=g[g.shared_zones.map(lambda v:zone in parts(v))];lost=g[~g.index.isin(keep.index)]
  zones.append(dict(pair_id=k,MPA=g.MPA.iloc[0],taxon=g.taxon.iloc[0],method=g.method.iloc[0],zone=zone,keep_years=';'.join(map(str,keep.year)),lost_years=';'.join(map(str,lost.year)),n_keep=len(keep),n_lost=len(lost)))
summary=pd.DataFrame(out);summary.to_csv(O/'constant_zone_options.csv',index=False);pd.DataFrame(zones).to_csv(O/'individual_zone_coverage.csv',index=False)
# Actual approved fish observations: compare reference support, not alternate reference identities.
v=pd.read_csv('outputs/pisco_fish_rules_20261101/station_values.csv');base=pd.read_csv('outputs/pisco_fish_rules_20261101/approved_annual_input.csv');fish=base[base.Source.eq('PISCO')&base.Taxa.eq('S. pulcher')];v=v.merge(fish[['series_id','year']],on=['series_id','year'],validate='many_to_one')
ref=[];annual=[]
for sid,g in v.groupby('series_id'):
 rr=g[g.status.eq('reference')];years=set(g.year);sets={y:set(q.loc[q.value.notna(),'site']) for y,q in rr.groupby('year')};union=set.union(*sets.values());fixed=set.intersection(*[sets.get(y,set()) for y in years]);complete=[y for y in sorted(years) if sets.get(y,set())==union]
 ref.append(dict(series_id=sid,MPA=g.MPA.iloc[0],Resp=g.Resp.iloc[0],n_years=len(years),n_reference_stations=len(union),n_persistent_reference_stations=len(fixed),persistent_reference_stations=';'.join(sorted(fixed)),all_reference_stations=';'.join(sorted(union)),years_all_stations=';'.join(map(str,complete)),n_years_all_stations=len(complete),n_years_lost_requiring_all=len(years)-len(complete)))
 for y,q in rr.groupby('year'):
  av=q.value.mean();fv=q[q.site.isin(fixed)].value.mean() if fixed else np.nan
  annual.append(dict(series_id=sid,year=y,available_reference_mean=av,persistent_reference_mean=fv,n_available=q.value.notna().sum(),n_persistent=len(fixed),difference=fv-av))
pd.DataFrame(ref).to_csv(O/'fish_reference_coverage.csv',index=False);pd.DataFrame(annual).to_csv(O/'fish_reference_means.csv',index=False)
summary.groupby('program').agg(pair_method_taxon_records=('pair_id','size'),no_constant_zone=('constant_zones_all_years',lambda x:sum(not bool(v) for v in x)),all_years_unknown_zones=('unknown_zone_years',lambda x:0)).drop(columns='all_years_unknown_zones').to_csv(O/'program_screen.csv')
# Protect against alias inflation and accidental outcome-dependent loss.
assert not d.duplicated(['pair_id','year']).any()
assert all(z['n_keep']+z['n_lost']==len(d[d.pair_id.eq(z['pair_id'])]) for z in zones)
for f,h in json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text()).items():assert hashlib.sha256(Path(f).read_bytes()).hexdigest()==h
inputs=[I/(n+'.csv') for n in ['raw_station_year_audit','site_response_coverage','candidate_pair_coverage','comparison_eligibility']]+[Path('outputs/pisco_fish_rules_20261101')/n for n in ['station_values.csv','approved_annual_input.csv']]
(O/'input_hashes.json').write_text(json.dumps({str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in inputs},indent=2))
print('Pair/method/taxon records',len(summary),'pair-years',len(d));print(pd.DataFrame(ref).query('n_years_lost_requiring_all>0').to_string(index=False));print(summary.groupby('program').size())
u=pd.read_csv('outputs/pisco_fish_rules_20261101/survey_unit_decisions.csv');u=u[u.complete];links=v[['series_id','MPA','Resp','site','year','status']].drop_duplicates();joined=links.merge(u,on=['site','year'],validate='many_to_many');fy=[]
for (sid,y),g in joined.groupby(['series_id','year']):
 a=g[g.status.eq('mpa')];b=g[g.status.eq('reference')];common=set(a.zone)&set(b.zone)
 gaps=[]
 for zone in common:
  da=pd.to_datetime(a.loc[a.zone.eq(zone),'date']);db=pd.to_datetime(b.loc[b.zone.eq(zone),'date'])
  if len(da) and len(db):gaps.append(min(abs((x-z).days) for x in da for z in db))
 fy.append(dict(series_id=sid,MPA=g.MPA.iloc[0],Resp=g.Resp.iloc[0],year=y,shared_zones=';'.join(sorted(common)),protected_zones=';'.join(sorted(set(a.zone))),reference_zones=';'.join(sorted(set(b.zone))),protected_months=';'.join(map(str,sorted(set(pd.to_datetime(a.date).dt.month)))),reference_months=';'.join(map(str,sorted(set(pd.to_datetime(b.date).dt.month)))),max_zone_nearest_gap_days=max(gaps) if gaps else np.nan))
f=pd.DataFrame(fy);f.to_csv(O/'fish_complete_pass_zone_timing.csv',index=False)
fs=[]
for sid,g in f.groupby('series_id'):
 common=set.intersection(*[parts(x) for x in g.shared_zones]);allz=set.union(*[parts(x) for x in g.shared_zones]);fs.append(dict(series_id=sid,MPA=g.MPA.iloc[0],Resp=g.Resp.iloc[0],n_years=len(g),constant_shared_zones=';'.join(sorted(common)),all_shared_zones=';'.join(sorted(allz)),zone_set_changes=g.shared_zones.nunique()>1,years_gap_over30=int(g.max_zone_nearest_gap_days.gt(30).sum()),years_gap_over60=int(g.max_zone_nearest_gap_days.gt(60).sum()),years_gap_over90=int(g.max_zone_nearest_gap_days.gt(90).sum())))
pd.DataFrame(fs).to_csv(O/'fish_zone_timing_summary.csv',index=False)
print(pd.DataFrame(fs).query("Resp=='Den'").to_string(index=False))
