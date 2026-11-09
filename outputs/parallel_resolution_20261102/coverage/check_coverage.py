"""Independent read-only audit. All artifacts are written beside this script.
Run: python3 /tmp/donham_swarm_20261102/coverage/check_coverage.py
No repository code is imported or executed. No thresholds are adopted.
"""
from pathlib import Path
from itertools import permutations
import pandas as pd, numpy as np, json, hashlib
R=Path('/Users/adrianstier/sbc-donham-kelp-mpa-cascade')
D=Path('/Users/adrianstier/sbc-donham-kelp-mpa-data/data')
O=Path(__file__).resolve().parent
seen={}
def read(p,**kw):
 p=Path(p);p=p if p.is_absolute() else R/p
 seen[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
 return pd.read_csv(p,**kw)
def out(x,n):x.to_csv(O/n,index=False)
v=read('outputs/all_program_station_refits_20261102/verified_station_values.csv')
c=read('outputs/all_program_station_refits_20261102/coverage_summary.csv').set_index('series_id')
rows=[]
for sid,g in v.groupby('series_id'):
 q=g[g.status.eq('reference')&g.value.notna()]; ys=sorted(g.year.unique());post=sorted(g.loc[g.t>0,'year'].unique())
 sets={y:set(q.loc[q.year.eq(y),'site']) for y in ys};fixed=set.intersection(*sets.values());pf=set.intersection(*(sets[y] for y in post)) if post else set();union=set.union(*sets.values())
 a=dict(series_id=sid,Source=g.Source.iloc[0],Taxa=g.Taxa.iloc[0],Resp=g.Resp.iloc[0],histories=1,post_years=len(post),references=len(union),fixed_all=len(fixed),fixed_post=len(pf),changing_full=any(s!=union for s in sets.values()),changing_post=len({tuple(sorted(sets[y])) for y in post})>1,single_reference=len(union)==1,no_full_fixed=not fixed,no_post_fixed=not pf,post_years_full_fixed=len(post) if fixed else 0,post_years_post_fixed=len(post) if pf else 0,full_fixed_stations=';'.join(sorted(fixed)),post_fixed_stations=';'.join(sorted(pf)))
 for new,old in [('references','n_reference_stations'),('fixed_all','n_persistent_references'),('fixed_post','n_post_persistent'),('changing_full','reference_membership_changes'),('changing_post','post_reference_membership_changes')]:assert a[new]==c.loc[sid,old],(sid,new)
 rows.append(a)
x=pd.DataFrame(rows);out(x,'reference_checks.csv')
fields=['histories','post_years','changing_full','changing_post','single_reference','no_full_fixed','no_post_fixed','post_years_full_fixed','post_years_post_fixed']
out(x.groupby(['Source','Taxa','Resp'])[fields].sum().reset_index(),'reference_by_response.csv');out(x.groupby('Source')[fields].sum().reset_index(),'reference_by_program.csv')
# Recompute every station-based arithmetic ratio using supplied fixed offset.
a=read('outputs/all_program_station_refits_20261102/annual_comparisons.csv');offs=read('outputs/all_program_station_refits_20261102/offsets.csv').set_index('series_id').offset
errs=[]
for sid,g in v[v.t>0].groupby('series_id'):
 fixed=set(x.set_index('series_id').loc[sid,'full_fixed_stations'].split(';'))-{''}
 for scenario,ref in [('Available',g[g.status.eq('reference')]),('Persistent',g[g.status.eq('reference')&g.site.isin(fixed)])]:
  p=g[g.status.eq('mpa')].groupby('year').value.mean();r=ref.groupby('year').value.mean();calc=np.log((p+offs[sid])/(r+offs[sid])).dropna();expected=a[(a.series_id==sid)&(a.scenario==scenario)].set_index('year').lnDiff
  assert set(calc.index)==set(expected.index)
  errs.extend(abs(calc-expected).tolist())
# Recompute fish station-zone intersections from complete physical survey units.
u=read('outputs/pisco_fish_rules_20261101/survey_unit_decisions.csv');complete=u[u.complete]
fv=v[v.Source.eq('PISCO')&v.Taxa.eq('S. pulcher')]
zy=[]
for (sid,year),g in fv.groupby(['series_id','year']):
 sets=[set(complete.loc[(complete.site==r.site)&(complete.year==year),'zone']) for r in g.itertuples()]
 zy.append(dict(series_id=sid,year=year,zones=';'.join(sorted(set.intersection(*sets)))))
zy=pd.DataFrame(zy);zold=read('outputs/station_coverage_refits_20261101/station_zone_years.csv').fillna('');j=zy.merge(zold,on=['series_id','year'],validate='one_to_one');assert len(j)==len(zold) and (j.zones==j.zones_shared_every_station).all()
out(zy,'fish_zone_year_checks.csv')
# Raw method-preserving swath and fish event depth audit over retained station-years.
# All-month swath intersections outside six comparisons are diagnostics only.
raw=read(D/'PISCO/MLPA_kelpforest_swath_2024.csv',encoding='latin1',low_memory=False)
raw=raw[raw.campus.eq('VRG')|(raw.campus.eq('UCSB')&~raw.zone.isin(['DEEP','SPECIALO']))].copy()
fish=read(D/'PISCO/MLPA_kelpforest_fish_2024.csv',encoding='latin1',low_memory=False)
fish['date']=pd.to_datetime(fish[['year','month','day']],errors='coerce').dt.strftime('%Y-%m-%d')
fish['transect']=fish.transect.astype(str).str.replace(r'\.0$','',regex=True)
complete=complete.copy();complete['transect']=complete.transect.astype(str).str.replace(r'\.0$','',regex=True)
fish=fish.merge(complete,on=['site','year','date','method','zone','transect'],how='inner',validate='many_to_one')
ks=['site','year','month','day','method','zone','transect']
def events(df):
 df=df.copy();df['positive_depth']=df.depth.where(df.depth>0)
 e=df.groupby(ks,dropna=False).agg(depth=('positive_depth','median'),n_depth=('positive_depth','nunique')).reset_index()
 return e
se=events(raw);fe=events(fish);depth=[]
for sid,g in v[v.Source.eq('PISCO')].groupby('series_id'):
 e=(fe if g.Taxa.iloc[0]=='S. pulcher' else se).merge(g[['site','year','status']].drop_duplicates(),on=['site','year'],validate='many_to_one')
 sets=[set(e.loc[(e.site==r.site)&(e.year==r.year),'zone']) for r in g.itertuples()];shared=set.intersection(*sets)
 depth.append(dict(series_id=sid,Source='PISCO',Taxa=g.Taxa.iloc[0],Resp=g.Resp.iloc[0],station_years=len(g),physical_event_uses=len(e),missing_depth_event_uses=int(e.depth.isna().sum()),multi_depth_event_uses=int(e.n_depth.gt(1).sum()),allmonth_common_zones=';'.join(sorted(shared)),no_allmonth_constant_zone=not shared,depth_min=e.depth.min(),depth_max=e.depth.max(),scope='Complete approved fish units' if g.Taxa.iloc[0]=='S. pulcher' else 'Raw all-month method-filtered swath availability; not approved outside six comparisons'))
depth=pd.DataFrame(depth);out(depth,'pisco_raw_depth_by_history.csv')
out(depth.groupby(['Taxa','Resp']).agg(histories=('series_id','size'),histories_missing_depth=('missing_depth_event_uses',lambda s:int((s>0).sum())),histories_no_constant_zone=('no_allmonth_constant_zone','sum'),event_uses=('physical_event_uses','sum'),missing_depth_event_uses=('missing_depth_event_uses','sum')).reset_index(),'pisco_depth_by_response.csv')
# Verify all-month six-history shared-zone assertion directly from raw event availability.
sz=read('outputs/pisco_all_months_20261102/shared_zones.csv');assert set(depth[depth.series_id.isin(sz.series_id)].allmonth_common_zones)=={'INNER;MID'}
# Independent enumeration of maximum-cardinality/minimum-gap date pairing (no scipy matching).
vis=read('outputs/pisco_all_months_20261102/visit_values.csv');vis=vis[vis.scenario.eq('allmonths_sharedzones_visit_sizes')].copy();vis['date']=pd.to_datetime(vis.date)
pr=read('outputs/pisco_timing_donors_20261102/survey_date_pairs.csv');tc=read('outputs/pisco_timing_donors_20261102/timing_coverage.csv');ta=read('outputs/pisco_timing_donors_20261102/timing_annual.csv')
match_checks=[];cover=[]
for sid,g in vis.groupby('series_id'):
 dsid=sid.replace('|Bio','|Den');start=v[v.series_id==sid].MPA_Start.iloc[0]
 for tol in [0,14,30,60,90,366]:
  good=[]
  for year,z in g.groupby('year'):
   ps=sorted(z.loc[z.status=='mpa','site'].unique());rs=sorted(z.loc[z.status=='reference','site'].unique());covered=0
   for p in ps:
    for r in rs:
     for zone in ['INNER','MID']:
      pp=sorted(z.loc[(z.site==p)&(z.zone==zone),'date']);rr=sorted(z.loc[(z.site==r)&(z.zone==zone),'date'])
      # Exhaustively assign the smaller set to distinct positions in the larger set.
      left,right=(pp,rr) if len(pp)<=len(rr) else (rr,pp)
      possibilities=[]
      for perm in permutations(range(len(right)),len(left)):
       gaps=[abs((dt-right[j]).days) for dt,j in zip(left,perm)];valid=[d for d in gaps if d<=tol];possibilities.append((len(valid),sum(valid)))
      best=min(possibilities,key=lambda w:(-w[0],w[1])) if possibilities else (0,0)
      q=pr[(pr.series_id==dsid)&(pr.tolerance_days==tol)&(pr.year==year)&(pr.zone==zone)&(pr.protected_site==p)&(pr.reference_site==r)]
      assert (len(q),q.gap_days.sum())==best
      assert q.protected_date.is_unique and q.reference_date.is_unique
      assert (abs((pd.to_datetime(q.protected_date)-pd.to_datetime(q.reference_date)).dt.days)==q.gap_days).all()
      covered+=best[0]>0;match_checks.append(best)
   if covered==len(ps)*len(rs)*2 and year>start:good.append(year)
  old=tc[(tc.series_id==dsid)&(tc.tolerance_days==tol)].iloc[0];assert len(good)==old.paired_post_years
  cover.append(dict(series_id=dsid,tolerance_days=tol,post_years_retained=len(good),years=';'.join(map(str,good))))
out(pd.DataFrame(cover),'timing_checks.csv')
# Compare the density roster to the biomass roster used by the six-comparison review.
density_scope=[]
for sid in sz.series_id:
 b=v[v.series_id==sid];d=v[v.series_id==sid.replace('|Bio','|Den')];keys=['site','status','year']
 lost=d.merge(b[keys],on=keys,how='left',indicator=True);lost=lost[lost['_merge']=='left_only']
 density_scope.append(dict(series_id=sid.replace('|Bio','|Den'),full_density_station_years=len(d),six_review_station_years=len(b),omitted_density_station_years=len(lost),omitted_density_years=';'.join(map(str,sorted(set(d.year)-set(b.year)))),omitted_post_station_years=int((lost.t>0).sum()),omitted_prior_station_years=int((lost.t<0).sum()),omitted_designation_station_years=int((lost.t==0).sum())))
out(pd.DataFrame(density_scope),'density_scope_differences.csv')
# Documentation scope audit.
cd=read('outputs/dataset_inventory_20261101/comparison_decisions.csv');st=read('outputs/dataset_inventory_20261101/audit_status_by_comparison.csv');master=read('outputs/dataset_inventory_20261101/audit_master.csv')
six=set(sz.series_id);twelve=six|{s.replace('|Bio','|Den') for s in six};bugs=[]
for r in cd.itertuples():
 rules=r.rule_ids.split(';')
 if r.series_id in twelve and 'R04' not in rules:bugs.append(dict(series_id=r.series_id,issue='R04 missing for implemented density counterpart'))
 if r.Source=='PISCO' and r.Taxa=='S. pulcher' and not r.exclude and 'R05' not in rules:bugs.append(dict(series_id=r.series_id,issue='R05 missing for audited retained fish history'))
 if r.Source=='Landsat' and 'R03' in rules:bugs.append(dict(series_id=r.series_id,issue='Diver shared-zone rule attached to satellite polygon'))
out(pd.DataFrame(bugs),'rule_link_discrepancies.csv')
assert set(st.loc[st.survey_timing.str.contains('All valid months',na=False),'series_id'])==twelve
# Raw schema: confirm no diver-depth measurement field, do not equate MBON height with water depth.
schemas=[]
for f in ['MBON/SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv','MBON/SBCMBON_kelp_forest_integrated_fish_20231022.csv','LTER/Annual_fish_comb_20240307.csv','LTER/Annual_Kelp_All_Years_20240305.csv','LTER/Lobster_Abundance_All_Years_20230831.csv','ALL_sizefreq_2024.csv','MBON/KFM_Macrocystis_RawData_1984-2023.csv','LANDSAT/MPA_Runs_new.csv']:
 h=read(D/f,nrows=0,encoding='latin1');schemas.append(dict(file=f,columns=';'.join(h.columns),depth_named_fields=';'.join(c for c in h.columns if 'depth' in c.lower())))
out(pd.DataFrame(schemas),'source_schema_checks.csv')
# Survey-year labels versus exact date: raw duplicate label groups are evidence, not inferred campaigns.
r=raw[raw.site.isin(v.loc[v.Source.eq('PISCO'),'site'])][['site','year','survey_year','month','day','zone']].drop_duplicates();r['date']=pd.to_datetime(r[['year','month','day']]);campaign=r.groupby(['site','year','survey_year','zone']).agg(n_dates=('date','nunique'),first=('date','min'),last=('date','max')).reset_index();campaign['span_days']=(campaign['last']-campaign['first']).dt.days;out(campaign[campaign.n_dates>1],'survey_year_multi_date_groups.csv')
# Saved older proposal comparison (do not run or adopt it).
old=read('outputs/depth_reference_policy_20261101/mpa_decisions.csv');out(old[old.scenario.isin(['Measured_shallow_d90_n5','Measured_15_d90_n5','Measured_20_d90_n5','Zone_only_d90_n5'])],'historical_depth_proposal_coverage.csv')
for _,z in old[old.scenario.isin(['Measured_shallow_d90_n5','Measured_15_d90_n5','Measured_20_d90_n5'])].groupby('MPA'):
 assert z.years.nunique()==1 and z.reference.nunique()==1 and z.protected.nunique()==1
summary={'station_histories_checked':len(x),'reference_summary_all_fields_match':True,'annual_ratio_cells_checked':len(errs),'max_annual_ratio_error':max(errs),'fish_station_zone_years_checked':len(zy),'raw_pisco_histories_screened':len(depth),'timing_pair_strata_checked':len(match_checks),'timing_pair_rows_checked':len(pr),'timing_coverage_rows_checked':len(tc),'allmonth_status_scope_exactly_twelve_response_histories':True,'rule_link_flags':pd.DataFrame(bugs).issue.value_counts().to_dict(),'raw_survey_year_groups_multiple_dates':int((campaign.n_dates>1).sum()),'raw_survey_year_max_span_days':int(campaign.span_days.max())}
(O/'checks.json').write_text(json.dumps(summary,indent=2)+'\n');(O/'input_hashes.json').write_text(json.dumps(seen,indent=2)+'\n')
print(json.dumps(summary,indent=2));print(x.groupby('Source')[fields].sum().to_string());print(depth.groupby(['Taxa','Resp']).agg(histories=('series_id','size'),missing_depth_histories=('missing_depth_event_uses',lambda s:int((s>0).sum())),no_zone=('no_allmonth_constant_zone','sum')).to_string())
