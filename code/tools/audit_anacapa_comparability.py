"""Source-record audit; isolated diagnostics, never alters primary inputs."""
from pathlib import Path
import os,json,hashlib
import numpy as np
import pandas as pd
root=Path.cwd();dr=Path(os.environ.get('MPA_DATA_REPO',str(root.parent/'Donham-Stier-CA-MPA-Data-2026')))/'data';prior=root/'outputs/upstream_sites_20261101';out=root/'outputs/anacapa_comparability_20261101';out.mkdir(exist_ok=True)
files=[dr/'PISCO/MLPA_kelpforest_swath_2024.csv',dr/'PISCO/master_site_table_Emilyedit.csv',dr/'MBON/SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv',dr/'MBON/KFM_Macrocystis_RawData_1984-2023.csv',dr/'MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv']
(out/'source_hashes.json').write_text(json.dumps({str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in files},indent=2))
ps=pd.read_csv(files[0],encoding='latin1',low_memory=False);mapping=pd.read_csv(files[1]);mapping=mapping[mapping.CA_MPA_Name_Short=='Anacapa Island SMCA'].drop_duplicates('site');ks=pd.read_csv(files[4]);ks=ks[ks.CA_MPA_Name_Short=='Anacapa Island SMCA'].drop_duplicates('site_id');mapping=mapping[mapping.site.isin(ps.site)]
ks.to_csv(out/'kfm_stations.csv',index=False);mapping.to_csv(out/'pisco_stations.csv',index=False)
ps=ps[ps.site.isin(mapping.site)&(ps.campus=='UCSB')&~ps.zone.isin(['DEEP','SPECIALO'])&ps.month.between(5,10)].copy()
# Preserve actual survey events even if no adult kelp recorded.
event=['site','year','month','zone','transect'];vis=ps.groupby(event,as_index=False).agg(depth=('depth','median'),depth_min=('depth','min'),depth_max=('depth','max'),day=('day','min'));kelp=ps[ps.classcode=='MACPYRAD'].copy();assert kelp['count'].notna().all() and kelp['size'].notna().all();kelp['plants']=kelp['count'];kelp['stipes']=kelp['count']*kelp['size']
calpath=next((root/'data/external/sbc1275').glob('*.csv'));cal=pd.read_csv(calpath);cal=cal[cal.model==1].set_index('MONTH').slope;kelp['biomass']=kelp.stipes*kelp.month.map(cal)*1000
sums=kelp.groupby(event,as_index=False)[['plants','stipes','biomass']].sum();pv=vis.merge(sums,on=event,how='left',validate='one_to_one');pv[['plants','stipes','biomass']]=pv[['plants','stipes','biomass']].fillna(0)/60;pv=pv.merge(mapping[['site','site_status']],on='site',validate='many_to_one').rename(columns={'site_status':'status'});pv['Source']='PISCO';pv.to_csv(out/'pisco_transects.csv',index=False)
k=pd.read_csv(files[2],low_memory=False);k=k[k.site_id.isin(ks.site_id)&(k.taxon_name=='Macrocystis pyrifera')&(k.sample_method=='5mquad')].copy();k['count']=pd.to_numeric(k['count'],errors='raise');k['date']=pd.to_datetime(k.date);k['year']=k.date.dt.year;k['month']=k.date.dt.month;assert set(k.proj_taxon_id)=={'t-k-035','t-k-036'};k=k.merge(ks[['site_id','status']],on='site_id',validate='many_to_one')
keys=['site_id','year','date','month','status','replicate_id','area'];kv=k.groupby(keys,as_index=False)['count'].sum();assert not kv.duplicated(['site_id','year','replicate_id']).any();kv['plants']=kv['count']/kv.area;kv.to_csv(out/'kfm_quadrats.csv',index=False)
st=pd.read_csv(files[3]).drop_duplicates();st=st[st.SiteNumber.isin([11,31,32,34])].copy();st['site_id']='a-k-'+st.SiteNumber.astype(int).astype(str).str.zfill(2);st['date']=pd.to_datetime(st.SurveyDate,format='%m/%d/%Y %H:%M:%S');st['year']=st.date.dt.year;st['month']=st.date.dt.month;st['mass']=st.Stipe_Count*st.month.map(cal)*1000
st.groupby(['site_id','year'],as_index=False).agg(n=('Stipe_Count','size'),mean_stipes=('Stipe_Count','mean'),mean_mass=('mass','mean'),months=('month',lambda x:';'.join(map(str,sorted(x.unique()))))).to_csv(out/'kfm_size_pools.csv',index=False)
# Exact expectation of the existing empirical bootstrap (avoids Monte Carlo noise).
kannual=kv.groupby(['site_id','status','year'],as_index=False).agg(plants=('plants','mean'),quadrats=('plants','size'),count=('count','sum'));kannual=kannual.merge(st.groupby(['site_id','year'],as_index=False).agg(mean_stipes=('Stipe_Count','mean'),mean_mass=('mass','mean')),on=['site_id','year'],how='left',validate='one_to_one');kannual['stipes']=kannual.plants*kannual.mean_stipes;kannual['biomass']=kannual.plants*kannual.mean_mass;kannual.loc[kannual.plants==0,['stipes','biomass']]=0;kannual['Source']='KFM';kannual=kannual.rename(columns={'site_id':'site'});kannual.to_csv(out/'kfm_exact_annual.csv',index=False)
def pann(x):
 x=x.groupby(['site','status','year','zone','transect'],as_index=False)[['plants','stipes','biomass']].mean();x=x.groupby(['site','status','year','zone'],as_index=False)[['plants','stipes','biomass']].mean();return x.groupby(['site','status','year'],as_index=False)[['plants','stipes','biomass']].mean()
pa=pann(pv);pa['Source']='PISCO';pa.to_csv(out/'pisco_annual.csv',index=False)
check=pa.merge(pd.read_csv(prior/'pisco_selected_station_values.csv').query("MPA=='Anacapa Island SMCA' and Taxa=='M. pyrifera' and Resp=='Bio'"),on=['site','year','status']);assert len(check)>20 and np.max(abs(check.biomass-check.value))<1e-8
pd.concat([pa,kannual],ignore_index=True).to_csv(out/'station_annual.csv',index=False)
# Geometry: nominal coordinates measure proximity, not identical transects.
coords=[]
for _,a in ks.iterrows():
 for _,b in mapping.iterrows():
  lat1,lat2=np.radians([a.Lat,b.latitude]);dl=np.radians(b.longitude-a.Lon);dp=lat2-lat1;dist=6371*2*np.arcsin(np.sqrt(np.sin(dp/2)**2+np.cos(lat1)*np.cos(lat2)*np.sin(dl/2)**2));coords.append(dict(kfm_site=a.site_id,kfm_name=a.geolocation,pisco_site=b.site,kfm_status=a.status,pisco_status=b.site_status,distance_km=dist))
pd.DataFrame(coords).to_csv(out/'station_distances.csv',index=False)
# Prespecified descriptive restrictions; each retains both statuses in actual paired years.
scenarios={'All mapped stations':(kannual,pa),'Fixed Lighthouse reference':(kannual[kannual.site.isin(['a-k-31','a-k-34'])],pa[pa.site.str.contains('WEST_ISLE|LIGHTHOUSE')]),'Nearest coordinate stations':(kannual[kannual.site.isin(['a-k-31','a-k-11','a-k-34'])],pa[pa.site.isin(['ANACAPA_WEST_ISLE_CEN','ANACAPA_ADMIRALS_W','ANACAPA_LIGHTHOUSE_REEF_E'])]),'PISCO inner zone':(kannual,pann(pv[pv.zone=='INNER'])),'PISCO outer zone':(kannual,pann(pv[pv.zone=='OUTER'])),'PISCO 7–15 m':(kannual,pann(pv[pv.depth.between(7,15)]))}
floors={}
for source,d in [('KFM',kannual),('PISCO',pa)]:
 for variable in ['plants','stipes','biomass']:
  a=d.groupby(['year','status'])[variable].mean();floors[source,variable]=a[a>0].min()/2
# Bracket the depth restriction rather than choose a cutoff by its result.
for lo,hi in [(5,15),(6,15),(8,15),(7,12)]:
 scenarios[f'PISCO {lo}–{hi} m']=(kannual,pann(pv[pv.depth.between(lo,hi)]))
results=[];annual=[]
for scenario,(kd,pd_) in scenarios.items():
 for source,d in [('KFM',kd),('PISCO',pd_)]:
  for variable in ['plants','stipes','biomass']:
   a=d.groupby(['year','status'])[variable].mean().unstack();floor=floors[source,variable];a=a.loc[a.index.isin(range(2007,2013))].dropna(subset=['mpa','reference']);assert len(a)>0
   # Zero correction from full annual domain; these diagnostics do not redefine primary weights.
   a['lnRR']=np.log((a.mpa+floor)/(a.reference+floor));results.append(dict(scenario=scenario,Source=source,variable=variable,years=';'.join(map(str,a.index)),n=len(a),floor=floor,mean_lnRR=a.lnRR.mean(),mean_mpa=a.mpa.mean(),mean_reference=a.reference.mean(),zero_mpa=int((a.mpa==0).sum()),zero_reference=int((a.reference==0).sum())))
   a=a.reset_index();a['scenario']=scenario;a['Source']=source;a['variable']=variable;annual.append(a)
pd.DataFrame(results).to_csv(out/'diagnostic_contrasts.csv',index=False);pd.concat(annual).to_csv(out/'diagnostic_annual.csv',index=False)
# Dates and footprint breadth, including depths on all observed species (not only kelp presence).
pv[pv.year.between(2007,2012)].groupby(['site','year','status','zone'],as_index=False).agg(transects=('transect','size'),min_depth=('depth_min','min'),max_depth=('depth_max','max'),months=('month',lambda x:';'.join(map(str,sorted(x.unique()))))).to_csv(out/'pisco_coverage.csv',index=False)
kv[kv.year.between(2007,2012)].groupby(['site_id','year','status'],as_index=False).agg(quadrats=('replicate_id','size'),area=('area','sum'),plants=('count','sum'),date=('date','first')).to_csv(out/'kfm_coverage.csv',index=False)
# Match station counterpart month; count eligible same-month role/year comparisons.
month_rows=[]
for kid,pname in [('a-k-31','WEST_ISLE'),('a-k-11','ADMIRALS'),('a-k-34','LIGHTHOUSE')]:
 for year in range(2007,2013):
  km=set(kv[(kv.site_id==kid)&(kv.year==year)].month);pm=set(pv[pv.site.str.contains(pname)&(pv.year==year)].month);month_rows.append(dict(kfm=kid,pisco=pname,year=year,kfm_months=str(sorted(km)),pisco_months=str(sorted(pm)),overlap=';'.join(map(str,sorted(km&pm)))))
pd.DataFrame(month_rows).to_csv(out/'month_overlap.csv',index=False)
print(pd.DataFrame(results).query("variable=='biomass'").to_string(index=False));print('PISCO reconstruction maximum error:',np.max(abs(check.biomass-check.value)))

# Every scenario uses its source/response's original zero correction.
checks={'pisco_reconstruction_max_error':float(np.max(abs(check.biomass-check.value))), 'pisco_reconstructed_cells':len(check),'no_same_month_both_roles':not any(any(r['year']==year and r['kfm']=='a-k-31' and r['overlap'] for r in month_rows) and any(r['year']==year and r['kfm'] in ['a-k-11','a-k-34'] and r['overlap'] for r in month_rows) for year in range(2007,2013))}
protected=json.loads((prior/'protected_hashes.json').read_text())
checks['source_files_unchanged']=all(hashlib.sha256(Path(f).read_bytes()).hexdigest()==h for f,h in json.loads((out/'source_hashes.json').read_text()).items())
(out/'verification.json').write_text(json.dumps(checks,indent=2))
