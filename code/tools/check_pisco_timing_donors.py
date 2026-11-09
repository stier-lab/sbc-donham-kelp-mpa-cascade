"""Review-only date pairing and held-out size-donor validation. No primary changes."""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,ast
from scipy.optimize import linear_sum_assignment
from scipy.stats import t as student
O=Path('outputs/pisco_timing_donors_20261102');B=Path('outputs/pisco_all_months_20261102');D=Path('/Users/adrianstier/sbc-donham-kelp-mpa-data/data')
v=pd.read_csv(B/'visit_values.csv');v=v[v.scenario.eq('allmonths_sharedzones_visit_sizes')].copy();v['date']=pd.to_datetime(v.date);assert not v.duplicated(['series_id','site','year','zone','date']).any()
a=pd.read_csv(B/'annual_values.csv');a=a[a.scenario.eq('allmonths_sharedzones_visit_sizes')&a.series_id.str.endswith('|Den')].copy();a['baseline']=a.revised
of=pd.read_csv('outputs/all_program_station_refits_20261102/offsets.csv').set_index('series_id').offset
T=ast.parse(Path('code/tools/audit_station_coverage_refits.py').read_text());exec(compile(ast.Module(body=[n for n in T.body if isinstance(n,ast.FunctionDef) and n.name in ['fit','slope']],type_ignores=[]),'fit','exec'))
# Match maximum number of pairs within a tolerance, then minimize date gaps.
# Each date is used at most once in a physical station-pair/zone/year; shared sites
# are reused across comparisons, not counted as independent replicates.
pairs=[];coverage=[];annual=[];fits=[]
for sid,g in v.groupby('series_id'):
 dsid=sid.replace('|Bio','|Den');base=a[a.series_id==dsid]
 for tol in [0,14,30,60,90,366]:
  for yr,z in g.groupby('year'):
   ps=sorted(z.loc[z.status=='mpa','site'].unique());rs=sorted(z.loc[z.status=='reference','site'].unique());expected=len(ps)*len(rs)*2;means=[]
   for pp in ps:
    for rr in rs:
     for zone in ['INNER','MID']:
      p=z[(z.site==pp)&(z.zone==zone)].sort_values('date');r=z[(z.site==rr)&(z.zone==zone)].sort_values('date');gap=np.abs((p.date.to_numpy()[:,None]-r.date.to_numpy()[None,:])/np.timedelta64(1,'D'));cost=np.where(gap<=tol,gap,1e7);ii,jj=linear_sum_assignment(cost);keep=gap[ii,jj]<=tol;ii=ii[keep];jj=jj[keep]
      if len(ii):means.append((p.density.to_numpy()[ii].mean(),r.density.to_numpy()[jj].mean()))
      for i,j in zip(ii,jj):pairs.append(dict(series_id=dsid,tolerance_days=tol,year=yr,zone=zone,protected_site=pp,reference_site=rr,protected_date=str(p.iloc[i].date.date()),reference_date=str(r.iloc[j].date.date()),gap_days=int(gap[i,j]),protected_density=p.iloc[i].density,reference_density=r.iloc[j].density))
   b=base[base.year==yr].iloc[0];ok=len(means)==expected;new=np.log((np.mean([m[0] for m in means])+of[dsid])/(np.mean([m[1] for m in means])+of[dsid])) if ok else np.nan
   annual.append(dict(series_id=dsid,tolerance_days=tol,year=yr,t=b.t,baseline=b.baseline,paired=new,expected_station_zone_pairs=expected,matched_station_zone_pairs=len(means),complete=ok))
  q=pd.DataFrame(annual);q=q[(q.series_id==dsid)&(q.tolerance_days==tol)&(q.t>0)];z=q[q.complete];coverage.append(dict(series_id=dsid,tolerance_days=tol,original_post_years=len(q),paired_post_years=len(z),lost_years=';'.join(map(str,q.loc[~q.complete,'year']))))
  for label,qq,col in [('baseline_all',q,'baseline'),('baseline_same_years',z,'baseline'),('date_paired',z,'paired')]:
   ans=slope(qq,col)
   if ans:fits.append(dict(series_id=dsid,tolerance_days=tol,comparison=label,**ans))
pd.DataFrame(pairs).to_csv(O/'survey_date_pairs.csv',index=False);pd.DataFrame(annual).to_csv(O/'timing_annual.csv',index=False);pd.DataFrame(coverage).to_csv(O/'timing_coverage.csv',index=False);pd.DataFrame(fits).to_csv(O/'timing_effects.csv',index=False)
# Raw size samples, kept as frequency records, not deduplicated equal-sized animals.
r=pd.read_csv(D/'ALL_sizefreq_2024.csv',low_memory=False);r=r[r.campus.isin(['VRG','UCSB'])&r.classcode.isin(['MESFRA','STRPUR','PANINT'])&r.zone.isin(['INNER','MID'])].copy();r['Taxa']=r.classcode.map({'MESFRA':'M. franciscanus','STRPUR':'S. purpuratus','PANINT':'P. interruptus'});r=r[r['count'].gt(0)&r['size'].gt(0)&(r.Taxa.eq('P. interruptus')|r['size'].ge(25))];r['date']=pd.to_datetime(r[['year','month','day']],errors='coerce');r=r[r.date.notna()];r['mass']=np.select([r.Taxa.eq('M. franciscanus'),r.Taxa.eq('S. purpuratus')],[.00059*r['size']**2.917,.00059*r['size']**2.870],default=.001352821*r['size']**2.913963113);r['total_mass']=r.mass*r['count']
samples=r.groupby(['site','year','zone','date','Taxa'],as_index=False).agg(n_sized=('count','sum'),total_mass=('total_mass','sum'));samples['mean_mass']=samples.total_mass/samples.n_sized
mapping=v[['series_id','site','status','year']].drop_duplicates();mapping['Taxa']=mapping.series_id.str.split('|').str[1];s=mapping.merge(samples,on=['site','year','Taxa'],validate='many_to_many');s['sample_id']=s[['series_id','site','year','zone','date']].astype(str).agg('|'.join,axis=1);assert s.sample_id.is_unique
# Only validate sample cells coinciding with a counted positive species/zone visit.
counts=v[['series_id','site','year','zone','date','density']];targets=s.merge(counts,on=['series_id','site','year','zone','date'],validate='one_to_one');targets=targets[targets.density>0]
s.to_csv(O/'size_sample_inventory.csv',index=False)
def predict(target,kind,tol):
 d=s[(s.series_id==target.series_id)&(s.year==target.year)]
 if kind=='same_site_other_visit':d=d[(d.site==target.site)&(d.zone==target.zone)&(d.date!=target.date)]
 elif kind=='same_site_other_zone':d=d[(d.site==target.site)&(d.zone!=target.zone)&(d.date==target.date)]
 elif kind=='same_role_other_site':d=d[(d.site!=target.site)&(d.status==target.status)&(d.zone==target.zone)]
 elif kind=='other_role_other_site':d=d[(d.site!=target.site)&(d.status!=target.status)&(d.zone==target.zone)]
 elif kind=='same_role_other_site_any_zone':d=d[(d.site!=target.site)&(d.status==target.status)]
 d=d.copy();d['gap']=(d.date-target.date).abs().dt.days;d=d[d.gap<=tol]
 if d.empty:return None
 # Nearest eligible visit in each donor site, then equal site weights.
 d=d[d.gap.eq(d.groupby('site').gap.transform('min'))]
 site=d.groupby('site',as_index=False).agg(total_mass=('total_mass','sum'),n=('n_sized','sum'));pred=np.mean(site.total_mass/site.n)
 assert kind in ['same_site_other_visit','same_site_other_zone'] or not (d.site==target.site).any()
 return dict(predicted=pred,donor_sites=';'.join(sorted(d.site.unique())),donor_samples=len(d),donor_animals=d.n_sized.sum(),max_gap_days=d.gap.max())
rows=[];attempts=[]
for t in targets.itertuples():
 for kind in ['same_site_other_visit','same_site_other_zone','same_role_other_site','other_role_other_site','same_role_other_site_any_zone']:
  for tol in ([0] if kind=='same_site_other_zone' else [30,60,90,366]):
   ans=predict(t,kind,tol);attempts.append(dict(sample_id=t.sample_id,series_id=t.series_id,strategy=kind,tolerance_days=tol,available=ans is not None))
   if ans:
    err=np.log(ans['predicted']/t.mean_mass);rows.append(dict(sample_id=t.sample_id,series_id=t.series_id,Taxa=t.Taxa,site=t.site,status=t.status,year=t.year,zone=t.zone,date=str(t.date.date()),target_n=t.n_sized,observed=t.mean_mass,strategy=kind,tolerance_days=tol,log_error=err,absolute_log_error=abs(err),**ans))
pr=pd.DataFrame(rows);att=pd.DataFrame(attempts);pr.to_csv(O/'donor_predictions.csv',index=False);att.to_csv(O/'donor_validation_coverage.csv',index=False)
summary=[]
for (tx,kind,tol),z in pr.groupby(['Taxa','strategy','tolerance_days']):
 eligible=targets[targets.Taxa==tx];summary.append(dict(Taxa=tx,strategy=kind,tolerance_days=tol,targets_available=len(eligible),targets_predicted=len(z),target_site_years=z[['site','year']].drop_duplicates().shape[0],median_target_n=z.target_n.median(),median_signed_factor=np.exp(z.log_error.median()),median_absolute_fold_error=np.exp(z.absolute_log_error.median()),fraction_within_25pct=(np.abs(np.exp(z.log_error)-1)<=.25).mean(),fraction_within_50pct=(np.abs(np.exp(z.log_error)-1)<=.5).mean()))
pd.DataFrame(summary).to_csv(O/'donor_validation_summary.csv',index=False)
# Like-for-like performance subset for methods with overlapping held-out samples.
common=[]
for tol in [30,60,90,366]:
 methods=['same_site_other_visit','same_role_other_site','other_role_other_site'];parts=[set(pr.loc[(pr.strategy==k)&(pr.tolerance_days==tol),'sample_id']) for k in methods];ids=set.intersection(*parts)
 for (tx,k),z in pr[(pr.sample_id.isin(ids))&(pr.strategy.isin(methods))&(pr.tolerance_days==tol)].groupby(['Taxa','strategy']):common.append(dict(Taxa=tx,strategy=k,tolerance_days=tol,n_common_targets=len(z),median_absolute_fold_error=np.exp(z.absolute_log_error.median())))
pd.DataFrame(common).to_csv(O/'donor_common_targets.csv',index=False)
# Availability for truly unsupported positive visits, not imputed biomass.
missing=pd.read_csv(B/'size_matching.csv');missing=missing[missing.retained_shared_zone&missing.density.gt(0)&missing.match.str.startswith('no ')].copy();missing['date']=pd.to_datetime(missing.count_date);avail=[]
for t in missing.itertuples():
 for tol in [30,60,90,366]:
  ans=predict(t,'same_role_other_site',tol);avail.append(dict(series_id=t.series_id,site=t.site,status=t.status,year=t.year,zone=t.zone,count_date=t.count_date,tolerance_days=tol,donors_available=ans is not None,max_gap_days=ans['max_gap_days'] if ans else np.nan,donor_sites=ans['donor_sites'] if ans else ''))
pd.DataFrame(avail).to_csv(O/'missing_visit_donor_availability.csv',index=False)
protected=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text());assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in protected.items())
info=dict(timing_fits=len(fits),paired_date_rows=len(pairs),heldout_targets=len(targets),predictions=len(pr),missing_positive_visits=len(missing),protected_unchanged=len(protected));(O/'verification.json').write_text(json.dumps(info,indent=2));print(json.dumps(info,indent=2));print(pd.DataFrame(coverage).pivot(index='series_id',columns='tolerance_days',values='paired_post_years').to_string());print(pd.DataFrame(summary).query("strategy=='same_role_other_site' and tolerance_days==90").to_string(index=False))
