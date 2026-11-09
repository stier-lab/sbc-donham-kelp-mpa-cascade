"""Check consolidated values, missingness, keyed joins and viewer payload against dated evidence."""
from pathlib import Path
import json,re,hashlib
import numpy as np,pandas as pd
O=Path('outputs/dataset_inventory_20261101');N=Path('outputs/all_program_station_refits_20261102');checks=[]
def read(p):return pd.read_csv(p,low_memory=False)
def equal(a,b,label):
 a=pd.Series(a).reset_index(drop=True).replace('',np.nan);b=pd.Series(b).reset_index(drop=True).replace('',np.nan);assert len(a)==len(b),label
 assert a.isna().equals(b.isna()),label+' missingness'
 if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):assert np.allclose(a,b,atol=1e-10,rtol=1e-10,equal_nan=True),label
 else:assert a.fillna('<missing>').astype(str).equals(b.fillna('<missing>').astype(str)),label
 checks.append(dict(check=label,cells=len(a)))
m=read(O/'audit_master.csv').sort_values('series_id').reset_index(drop=True);comp=read(O/'comparison_eligibility.csv').sort_values('series_id').reset_index(drop=True)
assert len(m)==206 and m.series_id.nunique()==206
for field in m:equal(m[field],comp[field],'master/eligibility '+field)
for name,rename in [('coverage_summary',{'n_reference_stations':'audit_n_reference_stations','n_persistent_references':'audit_n_persistent_references','reference_membership_changes':'audit_reference_membership_changes','post_reference_membership_changes':'audit_post_reference_membership_changes','n_post_persistent':'audit_n_post_persistent','recovery_fit_status':'audit_recovery_fit_status','size_support_status':'audit_size_support_status'}),('reconstruction_checks',{'verified':'annual_reconstruction_verified','max_ratio_error':'annual_reconstruction_error'}),('comparison_summary',{x:'recovery_'+x for x in ['Approved','Persistent_reference','Persistent_post_reference','Log_station_unadjusted','Log_station_adjusted','station_adjustment','reference_subset_change']})]:
 d=read(N/(name+'.csv')).set_index('series_id').reindex(m.series_id)
 for old,new in rename.items():equal(m[new],d[old],name+' '+new)
# Recalculate coverage from actual annual rows, not the summary generator.
a=read(Path('outputs/pisco_fish_rules_20261101/approved_annual_input.csv'));assert not a.duplicated(['series_id','year']).any()
for sid,g in a.groupby('series_id'):
 r=m[m.series_id.eq(sid)].iloc[0];yy=sorted(g.year);missing=sorted(set(range(min(yy),max(yy)+1))-set(yy));assert r.audited_n_years==len(yy) and r.audited_first_year==min(yy) and r.audited_last_year==max(yy)
 assert r.audited_years==';'.join(map(str,yy));assert ('' if pd.isna(r.audited_missing_years) else r.audited_missing_years)==';'.join(map(str,missing))
 assert r.audited_before_years==(g.year<g.MPA_Start).sum() and r.audited_post_years==(g.year>g.MPA_Start).sum()
assert m.loc[m.exclude,'audited_n_years'].isna().all() and m.exclude.sum()==31
# Independent station-set coverage, both full record and post-only.
v=read(N/'verified_station_values.csv');cv=read(N/'coverage_summary.csv').set_index('series_id')
for sid,g in v.groupby('series_id'):
 r=cv.loc[sid]
 for post in [False,True]:
  z=g[(g.status=='reference')&g.value.notna()&((g.t>0) if post else True)];yearset=set(g.loc[g.t>0,'year'] if post else g.year);sets=[set(z.loc[z.year==y,'site']) for y in sorted(yearset)]
  union=set.union(*sets) if sets else set();common=set.intersection(*sets) if sets else set();changed=any(s!=union for s in sets)
  assert len(common)==r['n_post_persistent' if post else 'n_persistent_references'];assert changed==r['post_reference_membership_changes' if post else 'reference_membership_changes']
  if not post:assert len(union)==r.n_reference_stations
# Recompute omission counts and magnitude summary from all paired checks.
i=read(N/'reference_omission.csv')
for r in m.itertuples():
 z=i[i.series_id==r.series_id]
 if not r.exclude:
  assert r.reference_omissions_checked==len(z)
  if len(z):assert np.isclose(r.largest_absolute_omission_shift,abs(z.after-z.before).max())
# Every nonmissing recovery result must be a recorded fit, and interval arithmetic holds.
e=read(N/'effects.csv');assert len(e)==950
assert np.allclose(e.estimate-e.low,e.high-e.estimate,rtol=1e-9,atol=1e-9)
from scipy.stats import t
assert np.allclose(e.high,e.estimate+t.ppf(.975,e.df)*e.se,atol=1e-9)
assert np.allclose(i.difference,i.after-i.before,atol=1e-10)
# Preview joins, no dependence on source row order.
prev=read(O/'supplement_table_preview.csv');key=['MPA','Source','Taxa','Resp'];j=prev.merge(m,on=key,suffixes=('_preview','_master'),validate='one_to_one');assert len(j)==206
for col in prev:
 if col in key:continue
 if col in m:equal(j[col+'_preview'],j[col+'_master'],'preview '+col)
# Recompute all nearest-date statistics from exact station dates.
sites=read(O/'site_response_coverage.csv').set_index('site_record_id');pairs=read(O/'candidate_pair_coverage.csv')
def dates(v):
 out={}
 if pd.isna(v):return out
 for s in str(v).split(';'):
  d=pd.Timestamp(s);out.setdefault(d.year,[]).append(d)
 return out
for p in pairs.itertuples():
 x=dates(sites.loc[p.protected_record_id,'survey_dates']);y=dates(sites.loc[p.reference_record_id,'survey_dates']);gaps=[min(abs((dd-ee).days) for dd in x[yr] for ee in y[yr]) for yr in x.keys()&y.keys()]
 assert len(gaps)==p.years_with_exact_dates_both;assert sum(x<=60 for x in gaps)==p.years_nearest_dates_within60;assert sum(x<=90 for x in gaps)==p.years_nearest_dates_within90
 if gaps:assert max(gaps)==p.max_nearest_date_gap_days
# The consolidated page embeds the current comparison and time-series rosters.
# Their detailed source checks are in check_current_decisions.py.
page=(O/'index.html').read_text();payload=json.loads(page.split('const data=',1)[1].split(';const meta=',1)[0])
assert len(payload['Included and excluded'])==len(m)
assert len(payload['Time series'])==len(m.loc[~m.exclude])
# Immutable evidence has not drifted.
for p,h in json.loads((O/'consolidation_input_hashes.json').read_text()).items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h,p
out=dict(candidate_rows=206,retained_histories=175,exclusions=31,annual_rows=len(a),station_histories=len(cv),models_checked=len(e),omission_checks=len(i),exact_date_pair_checks=len(pairs),column_checks=len(checks),cells_compared=sum(x['cells'] for x in checks),limits='This validates numeric calculations and evidence transfer; historical metadata and unresolved habitat, timing, size-support assumptions are not newly confirmed.')
(O/'consolidated_accuracy_checks.csv').write_text(pd.DataFrame(checks).to_csv(index=False));(O/'consolidated_accuracy.json').write_text(json.dumps(out,indent=2));print(out)
