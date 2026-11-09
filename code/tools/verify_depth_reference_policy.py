from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
O=Path('outputs/depth_reference_policy_20261101');e=pd.read_csv(O/'transect_depth_events.csv');v=pd.read_csv(O/'visit_coverage.csv');c=pd.read_csv(O/'closest_candidate_visits.csv');p=pd.read_csv(O/'pair_decisions.csv').fillna({'years':''});m=pd.read_csv(O/'mpa_decisions.csv').fillna({'years':'','protected':'','reference':''});sel=pd.read_csv(O/'selected_pair_zone_years.csv')
assert not e.duplicated(['MPA','site','year','survey_year','month','day','zone','transect']).any()
# Raw depths independently re-read. Dates/transects, not taxon records, define units.
h=json.loads((O/'source_hashes.json').read_text());source=next(k for k in h if k.endswith('MLPA_kelpforest_swath_2024.csv'));raw=pd.read_csv(source,encoding='latin1',low_memory=False);raw=raw[raw.site.isin(e.site)&raw.zone.isin(['INNER','MID'])];assert set(raw.campus)<= {'UCSB','VRG'}
raw['transect']=raw.transect.astype(str);e['transect']=e.transect.astype(str)
keys=['site','year','survey_year','month','day','zone','transect'];r=raw.groupby(keys).depth.apply(lambda x:sorted(set(float(y) for y in x.dropna() if np.isfinite(y) and y>0)))
for z in e.itertuples():
 vals=r.loc[tuple(getattr(z,k) for k in keys)];expected=float(np.median(vals)) if vals else np.nan;assert (np.isnan(expected) and np.isnan(z.depth_center)) or abs(expected-z.depth_center)<1e-10
assert (e.depth_n>=0).all();assert not sel[['MPA','scenario','site_p','site_r','year','zone']].duplicated().any()
for z in c.itertuples():
 assert z.year==pd.Timestamp(z.date_p).year==pd.Timestamp(z.date_r).year
 assert abs((pd.Timestamp(z.date_p)-pd.Timestamp(z.date_r)).days)==z.gap_days
 assert z.n_transects_p>=2 and z.n_transects_r>=2
# Reconstruct pair-year selection and eligibility from the candidate visit ledger.
for rr in p.itertuples():
 q=c[(c.MPA==rr.MPA)&(c.depth_rule==rr.depth_rule)&(c.site_p==rr.protected)&(c.site_r==rr.reference)]
 if rr.gap:q=q[q.gap_days<=rr.gap]
 candidates=[]
 for (year,campaign),g in q.groupby(['year','survey_year']):
  if set(g.zone)=={'INNER','MID'}:candidates.append((int(year),g.gap_days.max(),g.gap_days.sum(),int(campaign),g))
 best={}
 for year,mg,sg,cam,g in sorted(candidates,key=lambda x:(x[1],x[2],x[3])):best.setdefault(year,g)
 years=sorted(best);assert years==[int(x) for x in rr.years.split(';') if x]
 span=years[-1]-years[0] if years else 0;overlap=True
 z=pd.concat(list(best.values())) if best else q.iloc[:0]
 for zone in ['INNER','MID']:
  zz=z[z.zone==zone];pp=zz.depth_center_p.dropna();ref=zz.depth_center_r.dropna()
  if not len(pp) or not len(ref) or max(pp.min(),ref.min())>min(pp.max(),ref.max()):overlap=False
 assert overlap==rr.paired_depth_support_overlap
 ok=len(years)>=rr.min_years and span>=(5 if rr.min_years==5 else 2) and (overlap if rr.depth_rule.startswith('Measured') else True)
 assert ok==rr.eligible
for rr in m.itertuples():
 q=p[(p.MPA==rr.MPA)&(p.scenario==rr.scenario)&p.eligible];years=sorted(set(y for ss in q.years for y in map(int,ss.split(';'))));assert len(q)==rr.n_pairs;assert len(years)==rr.n_years;assert years==[int(x) for x in rr.years.split(';') if x];assert set(q.protected)==set(rr.protected.split(';'))- {''};assert set(q.reference)==set(rr.reference.split(';'))-{''}
# No novel expert-reference assignment is invented by selection.
master=pd.read_csv('outputs/cross_site_seasons_20261101/candidate_stations.csv');assigned=set(zip(master.MPA,master.site,master.status))
for z in p.itertuples():assert (z.MPA,z.protected,'mpa') in assigned and (z.MPA,z.reference,'reference') in assigned
protected=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text())
for file,digest in {**protected,**h}.items():assert hashlib.sha256(Path(file).read_bytes()).hexdigest()==digest,file
for directory in ['cross_site_seasons_20261101','point_vicente_seasons_20261101']:
 for file,digest in json.loads((Path('outputs')/directory/'manifest.json').read_text()).items():assert hashlib.sha256(Path(file).read_bytes()).hexdigest()==digest,file
res=dict(raw_depth_events_checked=len(e),candidate_date_matches_checked=len(c),pair_policy_decisions_checked=len(p),mpa_policy_rows_checked=len(m),novel_reference_assignments=0,protected_files_unchanged=len(protected),previous_audits_unchanged=True,source_hashes_match=True,outcome_models_fitted=0)
(O/'verification.json').write_text(json.dumps(res,indent=2)+'\n');print(res)
