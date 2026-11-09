"""Conditional joint measurement resampling; immutable primary inputs.
Permanent KFM plots/LTER transects share weights through time. PISCO units
share weights within site/year/zone/transect. Source size records retain identity
across all uses. Draw ranges are NOT total confidence intervals.
"""
from pathlib import Path
import ast,json,hashlib
import numpy as np,pandas as pd
from scipy import sparse
from scipy.stats import t as student
from biomass_source_guard import source_issue,require_allowed
O=Path('outputs/biomass_joint_checks_20261102');B=Path('outputs/biomass_support_20261102');S=Path('outputs/all_program_station_refits_20261102');D=Path('/Users/adrianstier/sbc-donham-kelp-mpa-data/data')
N=1000;rng=np.random.default_rng(20261102)
a=pd.read_csv('outputs/biomass_resolution_20261102/review_annual_input.csv');v=pd.read_csv(S/'verified_station_values.csv');off=pd.read_csv(S/'offsets.csv').set_index('series_id').offset
T=ast.parse(Path('code/tools/audit_station_coverage_refits.py').read_text());exec(compile(ast.Module(body=[n for n in T.body if isinstance(n,ast.FunctionDef) and n.name in ['fit','slope']],type_ignores=[]),'fit','exec'))
taxa={'Macrocystis pyrifera':'M. pyrifera','Mesocentrotus franciscanus':'M. franciscanus','Strongylocentrotus purpuratus':'S. purpuratus','Panulirus interruptus':'P. interruptus','MESFRAAD':'M. franciscanus','STRPURAD':'S. purpuratus','PANINT':'P. interruptus'}
k=pd.read_csv(O/'private/kfm.sum.csv');k=k[k.taxon_name.ne('Macrocystis pyrifera')|k.area.eq(10)].copy();k['site']=k.site_id;k['Taxa']=k.taxon_name.map(taxa);k['stratum']='all';k['unit']='KFM|'+k.site+'|'+k.replicate_id.astype(str)+'|'+k.area.astype(str);k['Source']='KFM';k['density']=k['count']/k.area
l=pd.read_csv(O/'private/lter.sub.csv');l=l[l.taxon_name.ne('Macrocystis pyrifera')].copy();l['Taxa']=l.taxon_name.map(taxa);l['site']=l.site_id;l['Source']='LTER';l=l.groupby(['Source','site','year','Taxa','transect_id'],as_index=False).den.mean();l['density']=l.den;l['stratum']='all';l['unit']='LTER|'+l.site+'|'+l.transect_id.astype(str)
p=pd.read_csv(O/'private/Swath.site.PISCO.sum.int.csv',dtype={'transect':str});p=p[p.y.isin(['MESFRAAD','STRPURAD','PANINT'])].copy();p['Taxa']=p.y.map(taxa);p=p.groupby(['site','year','Taxa','zone','transect'],as_index=False)['count'].mean();p['density']=p['count']/60;p['Source']='PISCO';p['stratum']=p.zone;p['unit']='PISCO|'+p.site+'|'+p.year.astype(str)+'|'+p.zone+'|'+p.transect.astype(str)
cols=['Source','site','year','Taxa','stratum','unit','density'];units=pd.concat([k[cols],l[cols],p[cols]],ignore_index=True)
assert units.Taxa.notna().all();assert units.groupby(cols[:-1]).density.nunique().le(1).all();units=units.drop_duplicates(cols[:-1]).reset_index(drop=True);assert units.density.notna().all()
# A physical spatial unit receives one weight across every use and taxon.
uid,labels=pd.factorize(units.unit);uw=rng.exponential(size=(len(labels),N));units['uid']=uid
count_draw={};count_mean={};count_info=[]
for key,g in units.groupby(['Source','site','year','Taxa']):
 means=[];draw=[]
 for _,q in g.groupby('stratum'):
  w=uw[q.uid];means.append(q.density.mean());draw.append((q.density.to_numpy()@w)/w.sum(axis=0))
 m=np.mean(means);count_mean[key]=m;count_draw[key]=np.mean(draw,axis=0)/m if m>0 else np.ones(N)
 count_info.append(dict(zip(['Source','site','year','Taxa'],key),count_mean=m,n_count_units=g.unit.nunique(),n_strata=g.stratum.nunique()))
pd.DataFrame(count_info).to_csv(O/'count_unit_support.csv',index=False)
# Admit only count histories that reproduce the retained station density.
density=v[v.Resp.eq('Den')];checks=[];valid={}
for r in density.itertuples():
 key=(r.Source,r.site,r.year,r.Taxa)
 if key not in count_mean:continue
 error=count_mean[key]-r.value;ok=abs(error)<1e-7*max(1,abs(r.value));valid[key]=valid.get(key,True) and ok
 checks.append(dict(series_id=r.series_id,site=r.site,year=r.year,observed_density=r.value,reconstructed_density=count_mean[key],error=error,verified=ok))
pd.DataFrame(checks).to_csv(O/'count_reconstruction.csv',index=False)
# Animal donor draws; a single frequency row represents count animals, not one animal.
pm=pd.read_csv(B/'private/pool_members.csv',low_memory=False);pm['issue']=[source_issue(s,t) for s,t in zip(pm.site,pm.Taxa)];blocked=set(pm.loc[pm.issue.ne(''),'series_id']);pm=pm[~pm.series_id.isin(blocked)].copy()
pm['mass']=np.select([pm.Taxa.eq('M. franciscanus'),pm.Taxa.eq('S. purpuratus')],[.00059*pm['size']**2.917,.00059*pm['size']**2.870],default=.001352821*pm['size']**2.913963113)
keys=['series_id','year','status'];pools=pm[keys].drop_duplicates().reset_index(drop=True);pools['pool_id']=np.arange(len(pools));pm=pm.merge(pools,on=keys,validate='many_to_one')
raw=pm[['source_row','count','mass']].drop_duplicates().sort_values('source_row').reset_index(drop=True);assert not raw.source_row.duplicated().any();j=pd.Series(np.arange(len(raw)),index=raw.source_row);pm['j']=pm.source_row.map(j)
W=sparse.csr_matrix((pm['count'],(pm.pool_id,pm.j)),shape=(len(pools),len(raw)));means=(W@raw.mass.to_numpy())/np.asarray(W.sum(axis=1)).ravel();weights=rng.gamma(raw['count'].to_numpy()[:,None],1,size=(len(raw),N))/raw['count'].to_numpy()[:,None];sr=((W@(weights*raw.mass.to_numpy()[:,None]))/(W@weights))/means[:,None];del weights
size_ratio={tuple(r):sr[i] for i,r in enumerate(pools[keys].to_numpy())}
# Kelp plant identity check: exact duplicate removal is justified by unique marker keys.
plants=pd.read_csv(D/'MBON/KFM_Macrocystis_RawData_1984-2023.csv');nraw=len(plants);ndup=int(plants.duplicated().sum());plants=plants.drop_duplicates();assert not plants.duplicated(['SiteNumber','SurveyDate','Marker']).any()
plants['site']='a-k-'+plants.SiteNumber.astype(str).str.zfill(2);plants['date']=pd.to_datetime(plants.SurveyDate,format='mixed');plants['year']=plants.date.dt.year;plants['month']=plants.date.dt.month
cal=pd.read_csv('data/external/sbc1275/Macrocystis_Biomass_Relationships_20261016.csv');cal=cal[(cal.SP_CODE=='MAPY')&(cal.model==1)].set_index('MONTH').slope
plants['mass']=plants.Stipe_Count*plants.month.map(cal)*1000
plant_ratio={};plant_mean={};plant_info=[]
ks=pd.read_csv(O/'private/kfm.site.csv');ks=ks[ks.taxon_name.eq('Macrocystis pyrifera')&ks.area.eq(10)];ks['date']=pd.to_datetime(ks.date)
for key,g in plants.groupby(['site','year']):
 if source_issue(key[0],'M. pyrifera','KFM_Macrocystis_RawData_1984-2023.csv'):continue
 assert g.mass.notna().all() and g.mass.ge(0).all();w=rng.exponential(size=(len(g),N));m=g.mass.mean();plant_mean[key]=m;plant_ratio[key]=(g.mass.to_numpy()@w)/w.sum(axis=0)/m if m>0 else np.ones(N)
 dates=ks.loc[(ks.site_id==key[0])&(ks.year==key[1]),'date'].drop_duplicates();gaps=[abs((dates-d).dt.days).min() for d in g.date.unique()] if len(dates) else []
 plant_info.append(dict(site=key[0],year=key[1],plants=len(g),survey_dates=g.date.nunique(),mean_mass_g=m,max_nearest_count_date_gap=max(gaps) if gaps else np.nan))
pd.DataFrame(plant_info).to_csv(O/'kelp_plant_support.csv',index=False)
# Three paired resampling schemes on identical years. No temporal residual resampling.
effects=[];admission=[];baseline=[];deterministic=[];drawout=[]
pooled=set(pools.series_id)|blocked
candidates=v[(v.Resp=='Bio')&(v.series_id.isin(pooled)|((v.Source=='KFM')&(v.Taxa=='M. pyrifera')))]
for sid,allz in candidates.groupby('series_id'):
 z=a[(a.series_id==sid)&(a.t>0)].sort_values('year');base=slope(z,'lnDiff')
 if base is None:admission.append(dict(series_id=sid,admitted=False,reason='Fewer than three observed post years'));continue
 if sid in blocked:admission.append(dict(series_id=sid,admitted=False,reason='Excluded-site alias in size pool; awaiting original records'));continue
 vv=allz[allz.year.isin(z.year)];kelp=vv.Taxa.iloc[0]=='M. pyrifera';bad=[]
 for r in vv.itertuples():
  key=(r.Source,r.site,r.year,r.Taxa)
  if key not in count_draw:bad.append('Missing count units');continue
  if not kelp and not valid.get(key,False):bad.append('Count aggregation does not reproduce retained density')
  if kelp and r.value>0 and (r.site,r.year) not in plant_ratio:bad.append('Missing kelp plants')
 if bad:admission.append(dict(series_id=sid,admitted=False,reason='; '.join(sorted(set(bad)))));continue
 admission.append(dict(series_id=sid,admitted=True,reason='Conditional measurement audit; site/date representativeness not certified'))
 baseline.append(dict(series_id=sid,**base))
 for scheme in ['count_only','size_only','joint']:
  ys=[]
  for yr in z.year:
   vals={}
   for st in ['mpa','reference']:
    zz=vv[(vv.year==yr)&(vv.status==st)];arr=[]
    for r in zz.itertuples():
     key=(r.Source,r.site,r.year,r.Taxa);cr=count_draw[key] if scheme!='size_only' else np.ones(N)
     ss=plant_ratio.get((r.site,r.year),np.ones(N)) if kelp else size_ratio.get((sid,yr,st),np.ones(N))
     if r.value>0 and not kelp:assert (sid,yr,st) in size_ratio
     arr.append(r.value*cr*(ss if scheme!='count_only' else 1))
     if scheme=='joint' and kelp:
      deterministic.append(dict(series_id=sid,site=r.site,year=yr,current_biomass=r.value,expected_biomass=count_mean[key]*plant_mean.get((r.site,yr),0)))
    assert len(arr);vals[st]=np.mean(arr,axis=0)
   ys.append(np.log((vals['mpa']+off[sid])/(vals['reference']+off[sid])))
  ys=np.array(ys);X=np.column_stack([np.ones(len(z)),z.t]);b=(np.linalg.pinv(X)@ys)[1];effects.append(dict(series_id=sid,scheme=scheme,estimate=base['estimate'],temporal_low=base['low'],temporal_high=base['high'],draw_mean=b.mean(),measurement_draw_low=np.quantile(b,.025),measurement_draw_high=np.quantile(b,.975),measurement_draw_sd=b.std(ddof=1),n_years=len(z),draws=N))
  # Save bounded individual draws for an independent regression check.
  for d in range(3):
   for yy,tt,val in zip(z.year,z.t,ys[:,d]):drawout.append(dict(series_id=sid,scheme=scheme,draw=d,year=yy,t=tt,lnDiff=val,expected_slope=b[d]))
pd.DataFrame(effects).to_csv(O/'conditional_joint_results.csv',index=False);pd.DataFrame(admission).to_csv(O/'resampling_admission.csv',index=False);pd.DataFrame(baseline).to_csv(O/'baseline_models.csv',index=False);pd.DataFrame(drawout).to_csv(O/'private/verification_draws.csv',index=False);pd.DataFrame(deterministic).drop_duplicates().to_csv(O/'kelp_expected_biomass_check.csv',index=False)
protected=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text());assert all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in protected.items())
info=dict(draws=N,seed=20261102,admitted_histories=sum(r['admitted'] for r in admission),held_histories=sum(not r['admitted'] for r in admission),blocked_alias_histories=len(blocked),raw_kelp_rows=nraw,exact_duplicate_kelp_rows=ndup,unique_plant_keys=len(plants),count_checks=len(checks),count_check_failures=sum(not r['verified'] for r in checks),shared_size_rows=len(raw),protected_files_unchanged=len(protected))
(O/'verification.json').write_text(json.dumps(info,indent=2));print(json.dumps(info,indent=2));print(pd.DataFrame(admission).to_string(index=False))
