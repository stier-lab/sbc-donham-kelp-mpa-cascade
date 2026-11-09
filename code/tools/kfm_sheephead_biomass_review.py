"""Restricted KFM review: nonlinear joint count/size propagation; no primary edits."""
from pathlib import Path
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import pandas as pd, numpy as np, json, hashlib
from scipy.stats import norm,rankdata
from scipy.stats import t
from scipy.stats import shapiro
from scipy.linalg import cholesky,solve_triangular
from joint_borrowed_biomass import fit
A=Path(__file__).resolve().parents[2]; D=A.parent/'sbc-donham-kelp-mpa-data'; O=A/'outputs/kfm_sheephead_biomass_20261103'; O.mkdir(exist_ok=True)
B=4000; rng=np.random.default_rng(20261103);copula_rng=np.random.default_rng(20261105)
paths=[A/'outputs/fish_lter_reconciliation_20261102/fish_count_size_matches.csv',D/'data/source_archive/fish_lter_reconciliation_20261102/chis_kelp_forest_fish_size.csv',A/'data/MBON/SBCMBON_kelp_forest_integrated_fish_20231022.csv',D/'config/kfm_physical_sites.csv']
roster=pd.read_csv(paths[0]); s=pd.read_csv(paths[1]);s=s[s.scientificName.str.contains('pulcher',na=False)].copy();s['source_row']=s.index+2
s['site_id']=s.siteName.map(pd.read_csv(paths[3]).set_index('source_name').site);s['year']=s.startDate.str[:4].astype(int);s['mass']=.0144*s.totalLengthInCentimeters**3.04
raw=pd.read_csv(paths[2],low_memory=False);raw=raw[raw.data_source.eq('kfm')&raw.sample_method.eq('visualfish')&raw.taxon_name.str.contains('pulcher',na=False)].copy();raw['count']=pd.to_numeric(raw['count'],errors='coerce');raw['year']=raw.date.str[:4].astype(int);raw['density']=raw['count']/raw.area
assert raw['count'].dropna().ge(0).all()
# Restrict every calibration record to the existing response-specific cohort.
import argparse
ap=argparse.ArgumentParser();ap.add_argument('--count-model',choices=['gamma','empirical'],default='gamma');ap.add_argument('--pool-strength',type=float,default=2);ap.add_argument('--count-prior',type=float,default=.5);ap.add_argument('--dependence',type=float,default=0);ap.add_argument('--variant',default='');args=ap.parse_args()
if args.variant:O=O/args.variant;O.mkdir(exist_ok=True)
s=s.merge(roster[['site_id','year']].drop_duplicates(),on=['site_id','year'],validate='many_to_one')
ss=0.;df=0
for _,g in s.groupby(['site_id','startDate']):
 n=g.speciesCount.to_numpy();m=g.mass.to_numpy();mu=np.average(m,weights=n);ss+=np.sum(n*(m/mu-1)**2);df+=n.sum()-1
pooled_cv2=ss/df
# Quasi-Poisson dispersion within count visits, calibrated on retained site-years only.
rc=raw.merge(roster[['site_id','year']].drop_duplicates(),on=['site_id','year'],validate='many_to_one')
pearson=0.;cdf=0
for _,g in rc.groupby(['site_id','date','transect_id']):
 if g['count'].isna().any():continue
 rate=g['count'].sum()/g.area.sum();expected=rate*g.area
 if rate>0 and len(g)>1:pearson+=np.sum((g['count']-expected)**2/expected);cdf+=len(g)-1
phi=max(1.,pearson/cdf)
size_cache={};count_draws={};mass_draws={};biomass_draws={};biomass_count_draws={};points={};ledger=[];size_use=[];visit_ledger=[]
for (site,year),rr in roster.groupby(['site_id','year']):
 r=rr.iloc[0];c=raw[raw.site_id.eq(site)&raw.year.eq(year)];assert c['count'].notna().all() and c.area.gt(0).all()
 g=s[s.site_id.eq(site)&s.year.eq(year)];dates=set(c.date);cp=c.groupby('transect_id').density.mean().mean();assert np.isclose(cp,r.count_density,atol=1e-12)
 assert c.transect_id.nunique()==1 # actual source structure, no fictitious independent transects
 cds=[];mds=[];bds=[];bcs=[];bps=[];mps=[];lps=[];hps=[];weights=[];modes=[];nfish=0;nrows=0;used=set()
 for date,cv in c.groupby('date'):
  v=g[g.startDate.eq(date)];mode='exact'
  if v.empty and not g.empty:
   gap=(pd.to_datetime(g.startDate)-pd.Timestamp(date)).dt.days.abs();v=g[gap.eq(gap.min())] if gap.min()<=7 else g.iloc[:0];mode='near_date' if len(v) else 'unsupported'
  if v.empty:mode='unsupported'
  modes.append(mode);density=cv.density.mean();wt=len(cv)/len(c);weights.append(wt)
  if args.count_model=='gamma':
   # Conditional quasi-Poisson rate posterior; Gamma shape increment; raw-count equivalent is phi times the increment.
   cd=rng.gamma(cv['count'].sum()/phi+args.count_prior,phi/cv.area.sum(),size=B)
  else:
   w=rng.exponential(size=(len(cv),B));cd=cv.density.to_numpy()@w/w.sum(axis=0)
  cds.append(cd)
  if v.empty:
   visit_ledger.append(dict(site_id=site,year=year,count_date=date,size_date='',mode=mode,count_density=density,n_fish=0));continue
  assert v.startDate.nunique()==1
  n=v.speciesCount.to_numpy();m=v.mass.to_numpy();mu=np.average(m,weights=n);nfish+=n.sum();nrows+=len(v)
  # Source row weights retain grouped fish frequencies; pseudo-record dispersion prevents singleton certainty.
  size_key=(site,v.startDate.iloc[0])
  if size_key in size_cache:md=size_cache[size_key].copy()
  else:
   w=rng.exponential(size=(len(v),B));md=(n*m)@w/(n@w)
   strength=args.pool_strength;vextra=strength*pooled_cv2/((max(len(v)-1,0)+strength)*max(len(v),1)) if strength else 0
   md*=np.exp(rng.normal(-.5*np.log1p(vextra),np.sqrt(np.log1p(vextra)),B))
   size_cache[size_key]=md.copy()
  if args.dependence:
   zsize=norm.ppf((rankdata(md)-.5)/B);latent=args.dependence*zsize+np.sqrt(1-args.dependence**2)*copula_rng.normal(size=B)
   aligned=np.empty(B);aligned[np.argsort(latent)]=np.sort(cd);cd=aligned;cds[-1]=cd
  mds.append(md);bds.append(cd*md);bcs.append(cd*mu);bps.append(density*mu);mps.append(mu)
  lps.append(density*np.average(.0144*v.minimumLengthInCentimeters.fillna(v.totalLengthInCentimeters)**3.04,weights=n));hps.append(density*np.average(.0144*v.maximumLengthInCentimeters.fillna(v.totalLengthInCentimeters)**3.04,weights=n))
  visit_ledger.append(dict(site_id=site,year=year,count_date=date,size_date=v.startDate.iloc[0],mode=mode,count_density=density,n_fish=int(n.sum()),mean_mass_g=mu,biomass_g_m2=density*mu,annual_weight=wt))
  for row in v.itertuples():
   if row.source_row not in used:size_use.append(dict(site_id=site,year=year,date=row.startDate,source_row=row.source_row,n_fish=row.speciesCount,mode=mode));used.add(row.source_row)
 count_draws[(site,year)]=np.average(cds,axis=0,weights=weights)
 supported=len(bds)==len(cds)
 if supported:
  biomass_draws[(site,year)]=np.average(bds,axis=0,weights=weights);biomass_count_draws[(site,year)]=np.average(bcs,axis=0,weights=weights);mass_draws[(site,year)]=np.average(mds,axis=0,weights=weights)
  bp=np.average(bps,weights=weights);ml=np.average(lps,weights=weights);mh=np.average(hps,weights=weights);meanmass=bp/cp if cp else np.average(mps,weights=weights)
  points[(site,year)]=(cp,meanmass,ml/cp if cp else meanmass,mh/cp if cp else meanmass)
 else:meanmass=np.nan;points[(site,year)]=(cp,np.nan,np.nan,np.nan)
 mode='unsupported' if not supported else 'near_date' if 'near_date' in modes else 'exact'
 ledger.append(dict(site_id=site,year=year,mode=mode,n_fish=int(s.loc[s.source_row.isin(used),'speciesCount'].sum()),n_size_rows=len(used),n_size_visits=s.loc[s.source_row.isin(used),'startDate'].nunique(),n_size_visit_uses=len(mds),count_density=cp,mean_mass_g=meanmass,n_transects=1,n_count_rows=len(c),supported=supported))
ledger=pd.DataFrame(ledger);ledger.to_csv(O/'site_year_support.csv',index=False);pd.DataFrame(size_use).to_csv(O/'size_record_use.csv',index=False);pd.DataFrame(visit_ledger).to_csv(O/'visit_support.csv',index=False)
s[['site_id','year','startDate','source_row','speciesCount']].to_csv(O/'variance_calibration_membership.csv',index=False)
rr=roster.merge(ledger[['site_id','year','mode','supported']],on=['site_id','year'],validate='many_to_one')
cov=rr.groupby(['MPA','year']).agg(complete=('supported','all'),n_sites=('site_id','size'),roles=('status','nunique')).reset_index();cov['eligible']=cov.complete&cov.roles.eq(2)&cov.year.ge(2007)
cov['expected_sites']=cov.MPA.map(rr.groupby('MPA').site_id.nunique());cov['all_expected_sites_present']=cov.n_sites.eq(cov.expected_sites)
cov.to_csv(O/'year_eligibility.csv',index=False)
annual=[];results=[];influence=[];diagnostics=[];fixed=[];bootchecks=[]
def estimate(y,S,years,rho):
 x=years-years.mean();X=np.column_stack([np.ones(len(x)),x]);scale=max(np.std(y),1e-8);K=rho**np.abs(years[:,None]-years[None,:]);f=fit(y/scale,X,S/scale**2,[K]);beta=f['beta']*scale;se=np.sqrt(np.diag(f['bc']))*scale;crit=t.ppf(.975,len(y)-2)
 residual=y-X@beta
 white=solve_triangular(cholesky(f['C']*scale**2,lower=True),residual,lower=True)
 return dict(normality_p=float(shapiro(white).pvalue),C=f['C']*scale**2,beta=beta,residual=residual,slope=beta[1],lower=beta[1]-crit*se[1],upper=beta[1]+crit*se[1],se=se[1],process_variance=f['q'][0]*scale**2,boundary=f['q'][0]<1e-8,upper_boundary=f['upper_boundary'],lag1=float(np.corrcoef(residual[:-1],residual[1:])[0,1]),fitted=X@beta)
for mpa,base in rr.groupby('MPA'):
 yy=cov[cov.MPA.eq(mpa)&cov.eligible].year.to_numpy();g=base[base.year.isin(yy)].copy();assert len(yy)>=5
 # Entire original role membership retained within each eligible year.
 for scenario in ['matched','exact_only','constant_stations']:
  z=g.copy()
  if scenario=='exact_only':
   bad=z[z['mode'].eq('near_date')].year.unique();z=z[~z.year.isin(bad)]
  if scenario=='constant_stations':
   for role in ['mpa','reference']:
    gg=z[z.status.eq(role)];keep=gg.groupby('site_id').year.nunique();keep=keep[keep.eq(gg.year.nunique())].index;z=z[~z.status.eq(role)|z.site_id.isin(keep)]
  years=np.sort(z.year.unique()); n=len(years);assert z.groupby('year').status.nunique().eq(2).all()
  vals={k:[] for k in ['density','biomass','low_mass','high_mass']};draws={k:[] for k in ['density','biomass','biomass_count_only']};means={k:[] for k in ['density','biomass']}
  for yr in years:
   rolevals={};roledraw={}
   for role in ['mpa','reference']:
    gg=z[z.year.eq(yr)&z.status.eq(role)];pt=[];bd=[];dd=[];bc=[]
    for row in gg.itertuples():
     key=(row.site_id,yr);cp,mu,lo,hi=points[key];pt.append([cp,0 if cp==0 else cp*mu,0 if cp==0 else cp*lo,0 if cp==0 else cp*hi]);dd.append(count_draws[key]);bd.append(biomass_draws[key]);bc.append(biomass_count_draws[key])
    rolevals[role]=np.mean(pt,axis=0);roledraw[role]=[np.mean(dd,axis=0),np.mean(bd,axis=0),np.mean(bc,axis=0)]
    if scenario=='matched':annual.append(dict(MPA=mpa,year=yr,status=role,n_sites=len(gg),sites=';'.join(gg.site_id),density=rolevals[role][0],biomass_g_m2=rolevals[role][1],size_low_biomass=rolevals[role][2],size_high_biomass=rolevals[role][3],density_lower=np.quantile(roledraw[role][0],.025),density_upper=np.quantile(roledraw[role][0],.975),biomass_lower=np.quantile(roledraw[role][1],.025),biomass_upper=np.quantile(roledraw[role][1],.975)))
   for j,k in enumerate(vals):vals[k].append(rolevals['mpa'][j]-rolevals['reference'][j])
   for j,k in enumerate(draws):draws[k].append(roledraw['mpa'][j]-roledraw['reference'][j])
   for j,k in enumerate(means):means[k].append((rolevals['mpa'][j]+rolevals['reference'][j])/2)
  for k in ['density','biomass_count_only','biomass']:
   response='biomass' if k.startswith('biomass') else 'density';y=np.array(vals[response]);S=np.cov(np.array(draws[k]));normalizer=np.mean(means[response]);assert normalizer>0
   for rho in [0.,.3,.6]:
    f=estimate(y,S,years,rho);factor=100 if response=='density' else .1 # fish/100m2 or kg/100m2
    row=dict(MPA=mpa,scenario=scenario,response=k,rho=rho,n_years=n,years=';'.join(map(str,years)),first_year=years.min(),last_year=years.max(),slope=f['slope']*factor,lower=f['lower']*factor,upper=f['upper']*factor,units='fish/100m2/year' if response=='density' else 'kg/100m2/year',scaled_slope_percent=100*f['slope']/normalizer,scaled_lower_percent=100*f['lower']/normalizer,scaled_upper_percent=100*f['upper']/normalizer,process_variance=f['process_variance'],boundary=f['boundary'],upper_boundary=f['upper_boundary'],residual_lag1=f['lag1'])
    results.append(row)
    if scenario=='matched' and rho==0:
     # Parametric studentized interval: refit process variance on each annual draw.
     # Conditional on estimated sampling covariance and Gaussian working process.
     brng=np.random.default_rng(20261104); scores=[]; boundary_n=0
     for bdraw in brng.multivariate_normal(f['fitted'],f['C'],size=199 if args.variant else 499):
      bf=estimate(bdraw,S,years,rho);scores.append((bf['slope']-f['slope'])/bf['se']);boundary_n+=bf['boundary']
     qlo,qhi=np.quantile(scores,[.025,.975]);row['bootstrap_lower']=(f['slope']-qhi*f['se'])*factor;row['bootstrap_upper']=(f['slope']-qlo*f['se'])*factor
     bootchecks.append(dict(MPA=mpa,response=k,n_bootstrap=199 if args.variant else 499,boundary_fits=boundary_n,normality_p=f['normality_p']))
     # Equal-year estimator separates size-driven changes from precision weighting.
     x=years-years.mean();w=x/(x@x);b=float(w@y)
     fixed.append(dict(MPA=mpa,response=k,slope=b*factor,scaled_slope_percent=100*b/normalizer,sampling_se=float(np.sqrt(w@S@w))*factor,interpretation='Equal-year point estimate; sampling SE excludes annual process uncertainty'))

     for i,yr in enumerate(years):
      take=np.arange(n)!=i; ff=estimate(y[take],S[np.ix_(take,take)],years[take],rho);influence.append(dict(MPA=mpa,response=k,omitted_year=yr,scaled_slope_percent=100*ff['slope']/normalizer))
     diagnostics.append(dict(MPA=mpa,response=k,min_cov_eigenvalue=float(np.linalg.eigvalsh(S).min()),bootstrap_draws=B,zero_role_years=sum(v[response]==0 for v in []),residual_lag1=f['lag1']))
  # Conservative source-range sensitivity, same observations; not confidence limits.
  for k in ['low_mass','high_mass']:
   y=np.array(vals[k]);f=estimate(y,np.cov(draws['biomass']),years,0.)
   results.append(dict(MPA=mpa,scenario=scenario,response=k,rho=0,n_years=n,years=';'.join(map(str,years)),slope=f['slope']*.1,lower=f['lower']*.1,upper=f['upper']*.1,units='kg/100m2/year',scaled_slope_percent=100*f['slope']/np.mean(means['biomass'])))
pd.DataFrame(fixed).to_csv(O/'equal_year_slopes.csv',index=False);pd.DataFrame(bootchecks).to_csv(O/'bootstrap_diagnostics.csv',index=False)
pd.DataFrame(annual).to_csv(O/'annual_responses.csv',index=False);pd.DataFrame(results).to_csv(O/'trends.csv',index=False);pd.DataFrame(influence).to_csv(O/'leave_one_year.csv',index=False);pd.DataFrame(diagnostics).drop(columns='zero_role_years').to_csv(O/'diagnostics.csv',index=False)
(O/'inputs.json').write_text(json.dumps([dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths],indent=2))
(O/'run.json').write_text(json.dumps(dict(seed=20261103,draws=B,pooled_relative_mass_variance=pooled_cv2,count_dispersion=phi,count_model=args.count_model,pool_strength=args.pool_strength,count_prior=args.count_prior,count_size_dependence=args.dependence,primary_changed=False),indent=2))
print(pd.DataFrame(results).query("scenario=='matched' and rho==0 and response in ['density','biomass_count_only','biomass']")[['MPA','response','n_years','slope','lower','upper','scaled_slope_percent']].to_string(index=False))
