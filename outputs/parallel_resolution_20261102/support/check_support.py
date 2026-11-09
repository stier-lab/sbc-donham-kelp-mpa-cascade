from pathlib import Path
import pandas as pd, numpy as np, json, hashlib
R=Path('/Users/adrianstier/sbc-donham-kelp-mpa-cascade'); D=Path('/Users/adrianstier/sbc-donham-kelp-mpa-data'); O=Path('/tmp/donham_swarm_20261102/support')
examined=[]
def read(p,**kw):
 p=Path(p);examined.append(p);return pd.read_csv(p,**kw)
b=R/'outputs/pisco_biomass_borrowing_20261102'; v=read(b/'visit_biomass_support.csv');v['date']=pd.to_datetime(v.date);v['Taxa']=v.series_id.str.split('|').str[1];u=v[v.method.eq('unsupported')].copy()
saved=read(R/'outputs/pisco_timing_donors_20261102/size_sample_inventory.csv');saved['date']=pd.to_datetime(saved.date)
raw=read(D/'data/ALL_sizefreq_2024.csv',low_memory=False);raw['source_row']=np.arange(1,len(raw)+1);raw['date']=pd.to_datetime(raw[['year','month','day']],errors='coerce')
codes={'M. franciscanus':'MESFRA','S. purpuratus':'STRPUR','P. interruptus':'PANINT'}
r=raw[raw.campus.isin(['VRG','UCSB'])&raw.classcode.isin(codes.values())].copy();r['Taxa']=r.classcode.map({v:k for k,v in codes.items()})
r=r[r['count'].gt(0)&r['size'].gt(0)&(r.Taxa.eq('P. interruptus')|r['size'].ge(25))&r.date.notna()&r.zone.isin(['INNER','MID'])].copy()
r['mass']=np.select([r.Taxa.eq('M. franciscanus'),r.Taxa.eq('S. purpuratus')],[.00059*r['size']**2.917,.00059*r['size']**2.870],default=.001352821*r['size']**2.913963113);r['total_mass']=r.mass*r['count']
physical=['site','year','zone','date','Taxa']; samples=r.groupby(physical,as_index=False).agg(n_sized=('count','sum'),total_mass=('total_mass','sum'));samples['mean_mass']=samples.total_mass/samples.n_sized
mapping=v[['series_id','site','status','year','Taxa']].drop_duplicates();s=mapping.merge(samples,on=['site','year','Taxa'],validate='many_to_many');keys=['series_id',*physical];q=s.merge(saved,on=keys,how='outer',indicator=True,validate='one_to_one');assert q._merge.eq('both').all();assert np.allclose(q.n_sized_x,q.n_sized_y);assert np.allclose(q.mean_mass_x,q.mean_mass_y)
assert not v.duplicated(['series_id','site','year','zone','date']).any()
cal_saved=read(b/'outer_site_calibration.csv');memo={}
def donor(t,kind,outer=None):
 z=s[s.series_id.eq(t.series_id)&s.status.eq(t.status)&s.year.eq(t.year)&s.zone.eq(t.zone)].copy()
 if outer is not None:z=z[z.site.ne(outer)]
 z=z[z.site.ne(t.site)] if kind=='other_site' else z[z.site.eq(t.site)&z.date.ne(t.date)]
 if z.empty:return z
 z['gap']=(z.date-t.date).abs().dt.days
 return z[z.gap.eq(z.groupby('site').gap.transform('min'))]
def calibr(t,kind):
 key=(t.series_id,t.site,t.status,kind)
 if key not in memo:
  rec=[]
  for c in s[s.series_id.eq(t.series_id)&s.status.eq(t.status)&s.site.ne(t.site)].itertuples():
   z=donor(c,kind,t.site)
   if len(z):rec.append((c.site,c.year,int(z.gap.max())))
  memo[key]=rec
 return memo[key]
rows=[];sources=[]
for t in u.itertuples():
 row=t._asdict();row.pop('Index');local=r[r.site.eq(t.site)&r.year.eq(t.year)&r.zone.eq(t.zone)&r.Taxa.eq(t.Taxa)]
 row['raw_local_year_animals']=local['count'].sum();row['raw_same_visit_animals']=local.loc[local.date.eq(t.date),'count'].sum();row['local_dates']=';'.join(sorted(set(local.date.dt.strftime('%Y-%m-%d'))))
 broad=raw[raw.site.eq(t.site)&raw.year.eq(t.year)&raw.zone.eq(t.zone)&raw.classcode.eq(codes[t.Taxa])]
 row['raw_all_local_rows']=len(broad);row['raw_local_invalid_date_rows']=int(broad.date.isna().sum());row['raw_local_under25_rows']=int((broad['size']<25).sum()) if t.Taxa!='P. interruptus' else 0
 aliases=raw[raw.site_new.eq(t.site)&raw.site.ne(t.site)&raw.year.eq(t.year)&raw.zone.eq(t.zone)&raw.classcode.eq(codes[t.Taxa])];row['alternate_site_key_rows']=len(aliases)
 for kind in ['other_visit','other_site']:
  z=donor(t,kind);cal=calibr(t,kind);nc=len(set((x[0],x[1]) for x in cal));maxcal=max((x[2] for x in cal),default=np.nan)
  reason='no donors' if z.empty else 'calibration guard' if nc<2 else 'date-gap guard' if z.gap.max()>maxcal else 'supported'
  row[kind+'_reason']=reason;row[kind+'_cal_siteyears']=nc;row[kind+'_cal_max_gap']=maxcal;row[kind+'_donor_gap']=z.gap.max() if len(z) else np.nan;row[kind+'_donors']=';'.join(sorted(z.site.unique()));row[kind+'_animals']=z.n_sized.sum()
  for a in z.itertuples():sources.append(dict(row_id=t.row_id,kind=kind,donor_site=a.site,donor_date=a.date,zone=a.zone,year=a.year,Taxa=a.Taxa,n_animals=a.n_sized,gap=a.gap))
  cs=cal_saved[cal_saved.series_id.eq(t.series_id)&cal_saved.outer_site.eq(t.site)&cal_saved.status.eq(t.status)&cal_saved.kind.eq(kind)]
  assert len(cs)==len(cal);assert cs.cluster.nunique()==nc
 row['root_reason']='no eligible local or same-role donor data' if row['other_visit_reason']==row['other_site_reason']=='no donors' else 'date-gap guard' if 'date-gap guard' in [row['other_visit_reason'],row['other_site_reason']] else 'calibration guard'
 rows.append(row)
out=pd.DataFrame(rows);annual=read(b/'biomass_annual.csv');out=out.merge(annual[['series_id','year','t']],on=['series_id','year'],validate='many_to_one');out['post']=out.t.gt(0)
# Independent completeness, without pandas skip-NA means.
blocked=set(map(tuple,out[['series_id','year']].to_numpy()));annual['recomputed_supported']=[(sid,yr) not in blocked for sid,yr in zip(annual.series_id,annual.year)];assert (annual.supported==annual.recomputed_supported).all()
post=annual[annual.t.gt(0)];assert len(post)==63 and post.recomputed_supported.sum()==44
out.to_csv(O/'unsupported_use_trace.csv',index=False);pd.DataFrame(sources).to_csv(O/'available_donor_evidence.csv',index=False)
blockedrows=[]
for (sid,year),z in out.groupby(['series_id','year']):
 blockedrows.append(dict(series_id=sid,year=year,post=bool(z.post.iloc[0]),unsupported_uses=len(z),physical_site_dates=len(z[['site','date']].drop_duplicates()),sites=';'.join(sorted(z.site.unique())),roles=';'.join(sorted(z.status.unique())),reasons=';'.join(sorted(z.root_reason.unique())),row_ids=';'.join(map(str,z.row_id))))
pd.DataFrame(blockedrows).to_csv(O/'blocked_comparison_years.csv',index=False)
# Reconstruct every positive unsupported count cell from the raw surveyed transect roster.
c=read(D/'data/PISCO/MLPA_kelpforest_swath_2024.csv',encoding='latin1',low_memory=False);c['date']=pd.to_datetime(c[['year','month','day']],errors='coerce');c=c[c.campus.eq('VRG')|(c.campus.eq('UCSB')&~c.zone.isin(['DEEP','SPECIALO']))]
checks=[]
for t in u.itertuples():
 cc=c[c.site.eq(t.site)&c.year.eq(t.year)&c.zone.eq(t.zone)&c.date.eq(t.date)]; k=['method','transect'];roster=cc[k].drop_duplicates();code=codes[t.Taxa]+('AD' if t.Taxa!='P. interruptus' else ''); counts=cc[cc.classcode.eq(code)].groupby(k,as_index=False)['count'].sum(min_count=1);j=roster.merge(counts,on=k,how='left',indicator=True,validate='one_to_one');j.loc[j._merge.eq('left_only'),'count']=0;den=j['count'].mean()/60
 assert np.isclose(den,t.density);checks.append(dict(row_id=t.row_id,raw_density=den,ledger_density=t.density,transects=len(j),raw_focal_count=j['count'].sum(),raw_focal_rows=len(cc[cc.classcode.eq(code)])))
pd.DataFrame(checks).to_csv(O/'raw_count_checks.csv',index=False)
# Confirm guide visit data reproduces source labels and identity, without edits.
guide=read(R/'outputs/dataset_inventory_20261101/measurement_visit_decisions.csv');print('Guide fields',list(guide))
summary={'unsupported_uses':len(out),'distinct_site_dates':len(out[['site','date']].drop_duplicates()),'distinct_site_zone_dates':len(out[['site','zone','date']].drop_duplicates()),'distinct_species_site_dates':len(out[['Taxa','site','date']].drop_duplicates()),'post_uses':int(out.post.sum()),'post_blocked_comparison_years':len(out[out.post][['series_id','year']].drop_duplicates()),'all_blocked_comparison_years':len(out[['series_id','year']].drop_duplicates()),'same_visit_recoverable':int(out.raw_same_visit_animals.gt(0).sum()),'local_other_date_available':int(out.raw_local_year_animals.gt(0).sum()),'raw_inventory_rows_reproduced':len(s),'root_reasons':out.root_reason.value_counts().to_dict(),'post_root_reasons':out[out.post].root_reason.value_counts().to_dict(),'alternate_site_key_rows':int(out.alternate_site_key_rows.sum())}
(O/'summary.json').write_text(json.dumps(summary,indent=2));(O/'examined_data_files.json').write_text(json.dumps([{'path':str(p),'resolved':str(p.resolve()),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in examined],indent=2))
print(json.dumps(summary,indent=2));print(out.groupby(['series_id','status','root_reason']).size());print(out[out.raw_local_year_animals.gt(0)][['row_id','site','year','Taxa','local_dates','other_visit_reason','other_visit_donor_gap','other_visit_cal_max_gap','other_site_reason']].to_string(index=False))
# Same-visit raw swath sizes and guide equality.
sw=[]
for t in u.itertuples():
 z=c[c.site.eq(t.site)&c.year.eq(t.year)&c.zone.eq(t.zone)&c.date.eq(t.date)&c.classcode.isin([codes[t.Taxa],codes[t.Taxa]+'AD'])]
 sw.append(dict(row_id=t.row_id,positive_direct_size_rows=int(z['size'].gt(0).sum())))
assert sum(x['positive_direct_size_rows'] for x in sw)==0
pd.DataFrame(sw).to_csv(O/'raw_swath_size_checks.csv',index=False)
g=guide.copy();g['date']=pd.to_datetime(g.date)
pd.testing.assert_frame_equal(v.drop(columns='Taxa'),g,check_dtype=False)
# Point Fermin alone: diagnostic local-date exception, preserve all sites/years/zones and exclusions.
t=u[u.row_id.eq(80)].iloc[0];z=donor(t,'other_visit');cal=cal_saved[cal_saved.series_id.eq(t.series_id)&cal_saved.outer_site.eq(t.site)&cal_saved.status.eq(t.status)&cal_saved.kind.eq('other_visit')].copy()
n=cal.groupby('cluster').size();weights=1/cal.cluster.map(n).to_numpy();order=np.argsort(cal.residual.to_numpy());res=cal.residual.to_numpy()[order];w=weights[order];probs=(np.cumsum(w)-.5*w)/w.sum();q=np.interp([.025,.5,.975],probs,res)
mass=(z.groupby('site').total_mass.sum()/z.groupby('site').n_sized.sum()).mean()*np.exp(q[1]);test=v[v.series_id.eq(t.series_id)&v.year.eq(t.year)].copy();test.loc[test.row_id.eq(80),'biomass_point']=t.density*mass
assert test.biomass_point.notna().all()
role=test.groupby(['status','site','zone']).biomass_point.mean().groupby(['status','site']).mean().groupby('status').mean()
newblocked=blocked-{(t.series_id,t.year)};recovered=int(sum((sid,yr) not in newblocked for sid,yr in zip(post.series_id,post.year)))
assert recovered==45
result=dict(scope='DIAGNOSTIC ONLY; no rule adopted or primary mutation',row_id=80,gap_days=294,calibration_max_days=290,excess_days=4,local_sized_animals=int(z.n_sized.sum()),calibration_siteyears=int(cal.cluster.nunique()),calibration_samples=len(cal),recovered_comparison='Abalone Cove SMCA|M. franciscanus|PISCO|Bio',recovered_year=2019,current_supported_post=44,diagnostic_supported_post=recovered,total_post=63,conditional_mass=mass,protected_annual_biomass=role['mpa'],reference_annual_biomass=role['reference'])
(O/'point_fermin_diagnostic.json').write_text(json.dumps(result,indent=2));cal.to_csv(O/'point_fermin_calibration.csv',index=False)
# Source rows supporting only the five same-site other-date candidates.
rr=[]
for t in out[out.raw_local_year_animals.gt(0)].itertuples():
 z=r[r.site.eq(t.site)&r.year.eq(t.year)&r.zone.eq(t.zone)&r.Taxa.eq(t.Taxa)].copy();z['unsupported_row_id']=t.row_id;rr.append(z)
pd.concat(rr).to_csv(O/'local_candidate_source_rows.csv',index=False)
print(json.dumps(result,indent=2))
