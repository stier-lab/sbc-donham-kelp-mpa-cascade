"""Within-visit biomass, source-linked donor pools, explicit unsupported visits.
Reviewed calibrated PISCO borrowing is imported as a bounded branch; its input
hashes must match. It is not silently generalized to new sites or years.
"""
from common import *
import sys
sys.path.insert(0,str(ROOT/'code/tools'))
from rebuild_borrowing_support import Engine,key
x=pd.read_parquet(PRIVATE/'measurement_uses.parquet')
s=read_tracked(DATA/'data/processed/source_build/analysis_sizes.csv');s['Taxa']=s.classcode.map(TAXA);s['Source']=s.campus.map({'UCSB':'PISCO','VRG':'PISCO','KFM':'KFM','LTER':'LTER'});s=s[s.Source.isin(['KFM','PISCO'])&s.Taxa.notna()&s['size'].gt(0)&s['count'].gt(0)].copy();s=s[~((s.Source=='PISCO')&s.Taxa.isin(['M. franciscanus','S. purpuratus'])&s['size'].lt(25))];s['zone']=s.zone.fillna('not_recorded');s['date']=pd.to_datetime(s.date,errors='coerce').fillna(pd.to_datetime(s[['year','month','day']],errors='coerce')).dt.strftime('%Y-%m-%d');s=s[s.date.notna()];s['mass']=[mass(t,z) for t,z in zip(s.Taxa,s['size'])]
s=s.loc[~pd.Series([blocked(a,b,t) for a,b,t in zip(s.Source,s.site,s.Taxa)],index=s.index)]
# Original kelp plant counts retain individual identity and actual month conversion.
k=read_source('kfm_kelp');k['_row']=np.arange(1,len(k)+1);k=k.drop_duplicates(subset=[c for c in k if c!='_row']);assert not k.duplicated(['SiteNumber','SurveyDate','Marker']).any();k['site']='a-k-'+k.SiteNumber.astype(str).str.zfill(2);k['date']=pd.to_datetime(k.SurveyDate,format='mixed').dt.strftime('%Y-%m-%d');k['year']=pd.to_datetime(k.date).dt.year;k['month']=pd.to_datetime(k.date).dt.month
cal=read_source('kelp_calibration');cal=cal[(cal.SP_CODE=='MAPY')&(cal.model==1)].set_index('MONTH').slope
k['mass']=k.Stipe_Count*k.month.map(cal)*1000;k['count']=1.;k['Taxa']='M. pyrifera';k['Source']='KFM';k['zone']='not_recorded';k['source_record_id']='kfm_kelp:'+SOURCES['kfm_kelp']['sha256']+':'+k['_row'].astype(str);k=k[~k.site.isin(KFM_EXCLUDED)&k.mass.ge(0)]
s=pd.concat([s,k],ignore_index=True);pcols=['Source','site','year','date','zone','Taxa'];pools=[];members=[]
for vals,g in s.groupby(pcols):
 pool=sid('size',*vals);mu=np.average(g.mass,weights=g['count']);n=g['count'].sum();pools.append(dict(zip(pcols,vals),pool_id=pool,mean_mass=mu,total_mass=float((g.mass*g['count']).sum()),n_sized=n,source_rows=len(g),singleton=len(g)==1,zero_empirical_variance=g.mass.nunique()==1))
 for r in g.itertuples():members.append(dict(pool_id=pool,source_record_id=r.source_record_id,mass=r.mass,frequency=r.count))
imp=read_tracked(PRIVATE/'lobster_imputation_members.csv')
for pool,g in imp.groupby('pool_id'):
 vals=dict(Source='LTER',site=g.site.iloc[0],year=int(g.date.iloc[0][:4]),date=g.date.iloc[0],zone='not_recorded',Taxa='P. interruptus')
 mu=np.average(g.mass,weights=g.frequency);pools.append(dict(**vals,pool_id=pool,mean_mass=mu,total_mass=float((g.mass*g.frequency).sum()),n_sized=g.frequency.sum(),source_rows=len(g),singleton=len(g)==1,zero_empirical_variance=g.mass.nunique()==1))
 members.extend(g[['pool_id','source_record_id','mass','frequency']].to_dict('records'))
p=pd.DataFrame(pools);p.to_csv(OUT/'size_pool_inventory.csv',index=False);pd.DataFrame(members).to_parquet(PRIVATE/'size_pool_members.parquet',index=False)
# Calibrated borrowing only for the six previously reviewed PISCO histories.
known=read_tracked(ROOT/'outputs/pisco_biomass_borrowing_20261102/visit_biomass_support.csv');approved=set(known.series_id);reviewed_targets=set(map(tuple,known[['series_id','site','year','zone','date']].drop_duplicates().itertuples(index=False,name=None)))
asg=read_tracked(OUT/'station_assignments.csv');asgb=asg[asg.included&asg.series_id.isin(approved)]
sp=p.merge(asgb[['series_id','site','role']].rename(columns={'role':'status'}),on='site',how='inner');sp['sample_key']=[key(r) for r in sp.itertuples()]
sp=sp[sp.Source.eq('PISCO')&sp.Taxa.eq(sp.series_id.str.split('|').str[1])]
cal_scope=read_tracked(ROOT/'outputs/pisco_timing_donors_20261102/size_sample_inventory.csv');sp=sp.merge(cal_scope[['series_id','site','status','year','zone','date']].drop_duplicates(),on=['series_id','site','status','year','zone','date'],how='inner',validate='one_to_one')
vis=x[x.Resp.eq('Bio')].groupby(['series_id','Source','site','status','year','date','zone','Taxa'],dropna=False).agg(density=('density',strict_mean)).reset_index();eng=Engine(sp,vis[vis.series_id.isin(approved)])
eng.calibration.to_csv(OUT/'borrowing_calibration.csv',index=False)
# Separate observed calibration scatter from size-sampling noise before adding
# donor draws. This avoids counting donor measurement error twice.
member_groups={k:g for k,g in pd.DataFrame(members).groupby('pool_id')};sample_lookup=sp.set_index('sample_key');caldraw={};variance_cache={};variance_rows=[]
def sample_draw(pool):
 if pool not in caldraw:
  g=member_groups[pool];w=np.array([np.random.default_rng(int(hashlib.sha256(str(r.source_record_id).encode()).hexdigest()[:15],16)).gamma(r.frequency,1,400) for r in g.itertuples()]);caldraw[pool]=g.mass.to_numpy()@w/w.sum(axis=0)
 return caldraw[pool]
def calibration_variance(series,site,role,kind):
 ck=(series,site,role,kind)
 if ck not in variance_cache:
  c=eng.cal[ck];e=[]
  for rr in c.itertuples():
   observed=sample_lookup.loc[rr.calibration_target_key];don=sample_lookup.loc[json.loads(rr.donor_keys)];by_site=[]
   for _,dd in don.groupby('site'):
    by_site.append(np.average([sample_draw(v.pool_id) for v in dd.itertuples()],axis=0,weights=dd.n_sized))
   predicted=np.mean(by_site,axis=0);e.append(np.log(sample_draw(observed.pool_id)/rr.observed)-np.log(predicted/rr.predicted))
  w=1/c.cluster.map(c.groupby('cluster').size()).to_numpy();w=w/w.sum();res=c.residual.to_numpy();center=np.sum(w*res);observed=float(np.sum(w*(res-center)**2));ee=np.array(e);ee=ee-(w@ee)[None,:];sampling=float(np.mean(np.sum(w[:,None]*ee**2,axis=0)));latent=max(0.,observed-sampling);shared=observed/c.cluster.nunique();variance_cache[ck]=(latent,shared,center)
  variance_rows.append(dict(series_id=series,site=site,status=role,kind=kind,calibration_groups=c.cluster.nunique(),observed_variance=observed,size_sampling_variance=sampling,latent_prediction_variance=latent,shared_calibration_variance=shared,mean_log_residual=center,latent_boundary=latent==0,scope='conditional empirical size draws; singleton variance unestimable'))
 return variance_cache[ck]
lookup={v:g for v,g in p.groupby(['Source','site','year','zone','Taxa'])};ledger=[];attempts=[]
for r in vis.itertuples():
 use=x.series_id.eq(r.series_id)&x.site.eq(r.site)&x.year.eq(r.year)&x.date.eq(r.date)&x.zone.eq(r.zone);z=x.loc[use];unknown=z.biomass.isna().any();chosen=[];basis='same_count_record';gap=0;error_key='';error_sd=0.;shared_sd=0.
 if unknown:
  local=lookup.get((r.Source,r.site,r.year,r.zone,r.Taxa),p.iloc[:0]);q=local[local.date.eq(r.date)]
  if q.empty and r.Source=='KFM' and r.Taxa in ['M. franciscanus','S. purpuratus','M. pyrifera'] and len(local):
   gaps=(pd.to_datetime(local.date)-pd.Timestamp(r.date)).abs().dt.days
   if gaps.min()<=POLICY['kfm_same_site_size_gap_days']:q=local[gaps.eq(gaps.min())];gap=int(gaps.min())
  if len(q):chosen=q.pool_id.tolist();basis='same_site_near_date_sizes' if gap else 'same_site_zone_date_sizes'
  elif r.series_id in approved and r.status=='reference' and (r.series_id,r.site,r.year,r.zone,r.date) in reviewed_targets:
   for kind in ['other_visit','other_site']:
    a=eng.attempt(r,kind,0);attempts.append(dict(series_id=r.series_id,**a))
    if a['eligible']:
     d=eng.donors(r,kind);chosen=d.pool_id.tolist();basis='calibrated_'+kind;gap=a['max_gap'];error_key='|'.join([r.series_id,r.site,r.status,kind]);latent,shared,center=calibration_variance(r.series_id,r.site,r.status,kind);error_sd=float(np.sqrt(latent));shared_sd=float(np.sqrt(shared));break
  if chosen:
   donor=p.set_index('pool_id').loc[chosen];mu=float(donor.mean_mass.mean())
   if basis.startswith('calibrated_'):mu=float(a['predicted']*np.exp(center+.5*latent))
   # Each target visit receives one mass estimate before any annual averaging.
   need=use&x.biomass.isna();x.loc[need,'biomass']=x.loc[need,'density']*mu;x.loc[need,'size_basis']=basis;x.loc[need,'donor_pool']=json.dumps(chosen)
  else:basis='unsupported_positive_count_or_missing_count'
 ledger.append(dict(series_id=r.series_id,Source=r.Source,site=r.site,status=r.status,year=r.year,date=r.date,zone=r.zone,Taxa=r.Taxa,density=r.density,biomass_supported=bool(x.loc[use,'biomass'].notna().all()),basis=basis,donor_pools=json.dumps(chosen),gap_days=gap,prediction_error_key=error_key,prediction_error_sd=error_sd,prediction_shared_sd=shared_sd,prediction_target_key='|'.join(map(str,[r.Source,r.site,r.year,r.zone,r.Taxa]))))
pd.DataFrame(variance_rows).to_csv(OUT/'borrowing_variance.csv',index=False);pd.DataFrame(ledger).to_csv(OUT/'visit_biomass_support.csv',index=False);pd.DataFrame(attempts).to_csv(OUT/'donor_attempts.csv',index=False);x.to_parquet(PRIVATE/'supported_measurement_uses.parquet',index=False);write_inputs();print('Biomass visits',len(ledger),'supported',sum(z['biomass_supported'] for z in ledger))
