"""Read-only bounded audit; all outputs beside this script. No primary fits."""
from pathlib import Path
import pandas as pd, numpy as np, json, hashlib, subprocess
ROOT=Path('/Users/adrianstier/sbc-donham-kelp-mpa-cascade')
DATA=Path('/Users/adrianstier/sbc-donham-kelp-mpa-data')
OUT=Path(__file__).resolve().parent
seen={}
def read(rel,root=ROOT,**kwargs):
 p=root/rel;seen[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();return pd.read_csv(p,**kwargs)
def save(x,name):x.to_csv(OUT/name,index=False)
p='outputs/pisco_biomass_borrowing_20261102/'
c=read(p+'outer_site_calibration.csv');v=read(p+'visit_biomass_support.csv');cv=read(p+'heldout_interval_validation.csv');s=read('outputs/pisco_timing_donors_20261102/size_sample_inventory.csv')
keys=['series_id','outer_site','status','kind']
cal=c.groupby(keys).agg(residual_rows=('residual','size'),cal_sites=('cal_site','nunique'),cal_years=('cal_year','nunique'),cal_siteyears=('cluster','nunique'),min_gap=('gap_days','min'),max_gap=('gap_days','max')).reset_index()
save(cal,'calibration_groups.csv')
b=v[v.method.str.startswith('borrowed')].copy();b['kind']=b.method.str.replace('borrowed_','');b=b.merge(cal,left_on=['series_id','site','status','kind'],right_on=keys,validate='many_to_one');b['below_min_calibration_gap']=b.max_gap_x<b.min_gap;b['above_max_calibration_gap']=b.max_gap_x>b.max_gap_y
save(b,'borrowed_support_diagnostic.csv')
assert not b.above_max_calibration_gap.any()
assert not (c.cal_site==c.outer_site).any()
assert all(r.outer_site not in r.donor_sites.split(';') for r in c.itertuples())
cv=cv[cv.available.astype(str).str.lower().eq('true')].copy();cv['hit']=cv.covered.astype(str).str.lower().eq('true');cv['siteyear']=cv.site+'|'+cv.year.astype(str)
validation=cv.groupby(['series_id','kind']).agg(target_visits=('hit','size'),heldout_sites=('site','nunique'),heldout_years=('year','nunique'),heldout_siteyears=('siteyear','nunique'),coverage=('hit','mean')).reset_index();save(validation,'validation_groups.csv')
# Reciprocal transfer rows are not separate evidence; reverse predictions induce +/- log ratios.
t=c[c.kind.eq('other_visit')];recip=[]
for key,g in t.groupby(keys+['cal_site','cal_year','cal_zone']):
 rr=np.sort(g.residual.to_numpy());recip.append(dict(zip(keys+['cal_site','cal_year','cal_zone'],key),rows=len(g),max_reciprocal_error=float(np.max(np.abs(rr+rr[::-1])))))
save(pd.DataFrame(recip),'temporal_reciprocity.csv')
# Reuse of the same donor among target sites in a species/year.
overlap=[]
for key,g in b[b.kind.eq('other_site')].groupby(['series_id','year']):
 sites={site:set(';'.join(z.donors).split(';')) for site,z in g.groupby('site')}
 for i,x in enumerate(sites):
  for y in list(sites)[i+1:]:
   common=sites[x]&sites[y]
   if common:overlap.append(dict(series_id=key[0],year=key[1],target_a=x,target_b=y,common_donors=';'.join(sorted(common))))
save(pd.DataFrame(overlap),'shared_donor_target_pairs.csv')
# Descriptive lag-one residual association only, never an independent significance test.
lag=[]
for key,g in c[c.kind.eq('other_site')].groupby(keys):
 yy=g.groupby(['cal_site','cal_year']).residual.mean().reset_index();pairs=yy.merge(yy.assign(cal_year=yy.cal_year+1),on=['cal_site','cal_year'],suffixes=('_now','_previous'))
 lag.append(dict(zip(keys,key),pairs=len(pairs),lag1_correlation=pairs.residual_now.corr(pairs.residual_previous) if len(pairs)>2 else np.nan))
save(pd.DataFrame(lag),'descriptive_temporal_correlation.csv')
# Independent raw mean/multiplicity rebuild for the 624 measured visit uses.
r=read('data/ALL_sizefreq_2024.csv',root=DATA,low_memory=False);r=r[r.campus.isin(['VRG','UCSB'])&r.classcode.isin(['MESFRA','STRPUR','PANINT'])].copy();r['Taxa']=r.classcode.map({'MESFRA':'M. franciscanus','STRPUR':'S. purpuratus','PANINT':'P. interruptus'});r=r[r['count'].gt(0)&r['size'].gt(0)&(r.Taxa.eq('P. interruptus')|r['size'].ge(25))];r['date']=pd.to_datetime(r[['year','month','day']],errors='coerce').dt.strftime('%Y-%m-%d');r['mass']=np.select([r.Taxa.eq('M. franciscanus'),r.Taxa.eq('S. purpuratus')],[.00059*r['size']**2.917,.00059*r['size']**2.870],default=.001352821*r['size']**2.913963113);r['total_mass']=r.mass*r['count'];r['mass2']=r.mass**2*r['count'];physical=['site','year','zone','date','Taxa'];local=r.groupby(physical).agg(n=('count','sum'),total_mass=('total_mass','sum'),mass2=('mass2','sum'),size_classes=('size','nunique')).reset_index();local['mean']=local.total_mass/local.n;local['bayes_mean_variance']=(local.mass2/local.n-local['mean']**2).clip(lower=0)/(local.n+1)
m=v[v.method.eq('measured_same_visit')].copy();m['Taxa']=m.series_id.str.split('|').str[1];m=m.merge(local,on=physical,validate='many_to_one');m['raw_mean_error']=m.mass-m['mean'];assert m.raw_mean_error.abs().max()<1e-8;save(m[['series_id']+physical+['n','size_classes','bayes_mean_variance','raw_mean_error']],'measured_sample_check.csv')
# Independent coefficient and saved-draw regression recomputation.
e=read(p+'biomass_effects.csv');a=read(p+'biomass_annual.csv');reg=[]
for row in e.itertuples():
 z=a[(a.series_id==row.series_id)&(a.t>0)&a.biomass_lnratio.notna()];X=np.c_[np.ones(len(z)),z.t];coef=np.linalg.lstsq(X,z.biomass_lnratio,rcond=None)[0];reg.append({'series_id':row.series_id,'slope_error':coef[1]-row.estimate,'years':len(z)})
assert max(abs(z['slope_error']) for z in reg)<1e-10
save(pd.DataFrame(reg),'slope_checks.csv')
for prefix in [p,'outputs/biomass_joint_checks_20261102/']:
 dd=read(prefix+'private/verification_draws.csv');kk=['series_id','draw']+(['scheme'] if 'scheme' in dd else []);errs=[]
 for _,z in dd.groupby(kk):errs.append(abs(np.linalg.lstsq(np.c_[np.ones(len(z)),z.t],z.lnDiff,rcond=None)[0][1]-z.expected_slope.iloc[0]))
 assert max(errs)<1e-10
# Kelp raw plant identity, visit replication and expected mass, all relevant count units.
j='outputs/biomass_joint_checks_20261102/'
k=read(j+'kelp_expected_biomass_check.csv');kp=read('data/MBON/KFM_Macrocystis_RawData_1984-2023.csv',root=DATA);raw_n=len(kp);dup=int(kp.duplicated().sum());kp=kp.drop_duplicates();assert not kp.duplicated(['SiteNumber','SurveyDate','Marker']).any();kp['site']='a-k-'+kp.SiteNumber.astype(str).str.zfill(2);kp['date']=pd.to_datetime(kp.SurveyDate,format='mixed');kp['year']=kp.date.dt.year;kp['month']=kp.date.dt.month
co=read('data/external/sbc1275/Macrocystis_Biomass_Relationships_20261016.csv');co=co[(co.SP_CODE=='MAPY')&(co.model==1)].set_index('MONTH');kp['mass']=kp.Stipe_Count*kp.month.map(co.slope)*1000
ks=read(j+'private/kfm.sum.csv');ks=ks[(ks.taxon_name=='Macrocystis pyrifera')&(ks.area==10)];ks['den']=ks['count']/ks.area;km=kp.groupby(['site','year']).agg(plants=('mass','size'),dates=('date','nunique'),mean_mass=('mass','mean')).reset_index();kd=ks.groupby(['site_id','year']).den.mean().reset_index().rename(columns={'site_id':'site'});z=k.merge(km,on=['site','year'],how='left').merge(kd,on=['site','year']);z['rebuilt_expected']=np.where(z.den==0,0,z.den*z.mean_mass);z['error']=z.rebuilt_expected-z.expected_biomass;assert z.error.abs().max()<1e-8;save(z,'kelp_independent_check.csv')
used=kp.merge(k[['site','year']].drop_duplicates(),on=['site','year']);ud=used.groupby(['site','year','date']).agg(plants=('mass','size'),mean_mass=('mass','mean')).reset_index();save(ud,'kelp_visit_support.csv')
# Bounded two-level weighting illustration for the four multi-visit site-years.
rng=np.random.default_rng(20261102);sens=[]
for key,g in used.groupby(['site','year']):
 if g.date.nunique()<2:continue
 codes,_=pd.factorize(g.date);w=rng.exponential(size=(len(g),2000));wc=rng.exponential(size=(g.date.nunique(),2000));iid=(g.mass.to_numpy()@w)/w.sum(0);two=(g.mass.to_numpy()@(w*wc[codes]))/(w*wc[codes]).sum(0)
 sens.append(dict(site=key[0],year=key[1],plants=len(g),dates=g.date.nunique(),plant_only_sd=iid.std(ddof=1),two_level_sd=two.std(ddof=1)))
save(pd.DataFrame(sens),'kelp_cluster_sensitivity.csv')
summary=dict(calibration_rows=len(c),calibration_models=len(cal),borrowed_visits=len(b),borrowed_temporal_visits=int(b.kind.eq('other_visit').sum()),temporal_below_min_gap=int(b.loc[b.kind.eq('other_visit'),'below_min_calibration_gap'].sum()),cross_site_below_min_gap=int(b.loc[b.kind.eq('other_site'),'below_min_calibration_gap'].sum()),overlapping_donor_target_pairs=len(overlap),measured_visits=len(m),singleton_measured_visits=int(m.n.eq(1).sum()),one_size_class_measured_visits=int(m.size_classes.eq(1).sum()),kelp_raw_rows=raw_n,kelp_exact_duplicates=dup,kelp_unique_plant_keys=len(kp),kelp_comparison_station_years=len(k),kelp_unique_count_station_years=len(k[['site','year']].drop_duplicates()),kelp_supported_plant_station_years=len(z.dropna(subset=['plants'])[['site','year']].drop_duplicates()),kelp_multi_visit_station_years=len(sens),raw_mean_max_error=m.raw_mean_error.abs().max(),kelp_expected_max_error=z.error.abs().max())
(OUT/'diagnostic_summary.json').write_text(json.dumps(summary,indent=2));(OUT/'examined_data_hashes.json').write_text(json.dumps(seen,indent=2));print(json.dumps(summary,indent=2));print(validation.to_string(index=False));print(pd.DataFrame(sens).to_string(index=False))
# Independent aggregation of source count exports to 4,236 saved reconstruction rows.
counts=read(j+'count_reconstruction.csv');tax={'Macrocystis pyrifera':'M. pyrifera','Mesocentrotus franciscanus':'M. franciscanus','Strongylocentrotus purpuratus':'S. purpuratus','Panulirus interruptus':'P. interruptus','MESFRAAD':'M. franciscanus','STRPURAD':'S. purpuratus','PANINT':'P. interruptus'}
kraw=read(j+'private/kfm.sum.csv');kraw=kraw[kraw.taxon_name.ne('Macrocystis pyrifera')|kraw.area.eq(10)].copy();kraw['Taxa']=kraw.taxon_name.map(tax);kraw['density']=kraw['count']/kraw.area;kd=kraw.groupby(['site_id','year','Taxa']).density.mean().reset_index().rename(columns={'site_id':'site'});kd['Source']='KFM'
lraw=read(j+'private/lter.sub.csv');lraw=lraw[lraw.taxon_name.ne('Macrocystis pyrifera')].copy();lraw['Taxa']=lraw.taxon_name.map(tax);ld=lraw.groupby(['site_id','year','Taxa','transect_id']).den.mean().groupby(['site_id','year','Taxa']).mean().reset_index().rename(columns={'site_id':'site','den':'density'});ld['Source']='LTER'
praw=read(j+'private/Swath.site.PISCO.sum.int.csv',dtype={'transect':str});praw=praw[praw.y.isin(['MESFRAAD','STRPURAD','PANINT'])].copy();praw['Taxa']=praw.y.map(tax);pdens=praw.groupby(['site','year','Taxa','zone','transect'])['count'].mean().groupby(['site','year','Taxa','zone']).mean().groupby(['site','year','Taxa']).mean().div(60).rename('density').reset_index();pdens['Source']='PISCO'
allc=pd.concat([kd,ld,pdens]);counts['Taxa']=counts.series_id.str.split('|').str[1];counts['Source']=counts.series_id.str.split('|').str[2];chk=counts.merge(allc,on=['Source','site','year','Taxa'],validate='many_to_one');assert len(chk)==len(counts);chk['independent_error']=chk.density-chk.observed_density;assert chk.independent_error.abs().max()<1e-7
save(chk,'independent_count_checks.csv');summary['independent_count_checks']=len(chk);summary['count_max_error']=chk.independent_error.abs().max();(OUT/'diagnostic_summary.json').write_text(json.dumps(summary,indent=2));(OUT/'examined_data_hashes.json').write_text(json.dumps(seen,indent=2));print('Independent count checks:',len(chk),'max error',summary['count_max_error'])
