from pathlib import Path
import pandas as pd,numpy as np,re,json,hashlib
O=Path('outputs/point_vicente_seasons_20261101');D=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data/PISCO');used=[]
def rd(p,**kw):used.append(str(p));return pd.read_csv(p,low_memory=False,**kw)
a=rd(Path('outputs/pisco_source_reconciliation_20261101/annual_input.csv'));a=a[a.scenario=='Supported_corrections'].copy();m='Point Vicente SMCA';old=a[(a.MPA==m)&(a.Source=='PISCO')];sw=rd(D/'MLPA_kelpforest_swath_2024.csv',encoding='latin1');fi=rd(D/'MLPA_kelpforest_fish_2024.csv',encoding='latin1');sw=sw[(sw.campus=='VRG')|((sw.campus=='UCSB')&~sw.zone.isin(['DEEP','SPECIALO']))];fi=fi[fi.campus.isin(['VRG','UCSB'])&(fi.level!='CAN')];master=rd(D/'master_site_table_Emilyedit.csv',encoding='utf-8-sig');master=master[(master.CA_MPA_Name_Short==m)&(master.BaselineRegion=='SOUTH')].rename(columns={'CA_MPA_Name_Short':'MPA','site_status':'status'})[['site','status']].drop_duplicates();assert not master.site.duplicated().any()
valid=pd.concat([sw[sw.month.between(5,10)][['site','year']],fi[['site','year']]]).drop_duplicates().groupby('site').size();utils=D.parent.parent/'code/R/01_utils.R';used.append(str(utils));exc=re.findall(r'"([^"]+)"',re.search(r'EXCLUDED_REFERENCE_SITES <- c\((.*?)\)',utils.read_text(),re.S).group(1));master=master[master.site.isin(valid[valid>=5].index)&~master.site.isin(exc)&master.status.isin(['mpa','reference'])];sw=sw.merge(master,on='site');sw=sw[sw.year.between(2013,2023)].copy()
ev=['site','status','year','month','day','zone','transect'];events=sw[ev].drop_duplicates();assert not events.duplicated(['site','year','month','zone','transect']).any()
taxa={'MACPYRAD':('M. pyrifera','Bio'),'MESFRAAD':('M. franciscanus','Den'),'STRPURAD':('S. purpuratus','Den'),'PANINT':('P. interruptus','Den')};f=sw[sw.classcode.isin(taxa)].copy();f['value']=f['count'];k=f.classcode=='MACPYRAD';f.loc[k,'value']=f.loc[k,'count']*f.loc[k,'size'];bad=f[f.value.isna()][ev+['classcode','count','size']];bad.to_csv(O/'missing_measurements.csv',index=False)
# A present but unmeasured record invalidates its response/event; absent taxa are zeros.
ct=f.groupby(ev+['classcode'],as_index=False).value.agg(lambda x:x.sum(skipna=False));ct['recorded']=True;g=events.merge(pd.DataFrame({'classcode':list(taxa)}),how='cross').merge(ct,on=ev+['classcode'],how='left');g.loc[g.recorded.isna(),'value']=0;g['Taxa']=g.classcode.map(lambda c:taxa[c][0]);g['Resp']=g.classcode.map(lambda c:taxa[c][1]);co=rd(Path('data/external/sbc1275/Macrocystis_Biomass_Relationships_20261016.csv'));co=co[co.model==1].set_index('MONTH').slope;g.value/=60;g.loc[g.Resp=='Bio','value']*=g.loc[g.Resp=='Bio','month'].map(co)*1000
# A zone is usable only when every observed transect/visit has all four focal measurements.
coverage=sw.groupby(['site','status','year','month','day','zone']).agg(transects=('transect','nunique'),depth_min=('depth','min'),depth_max=('depth','max')).reset_index();coverage.to_csv(O/'survey_dates_depths.csv',index=False)
g.to_csv(O/'event_responses.csv',index=False)
windows={'MayOct':(5,10),'AprOct':(4,10),'AprNov':(4,11),'MarNov':(3,11),'AllMonths':(1,12)};yearsets={'full11':list(range(2013,2024)),'old10':[2013]+list(range(2015,2024)),'late9':list(range(2015,2024))};zone_sets={'INNER_MID':['INNER','MID'],'INNER_MID_OUTER':['INNER','MID','OUTER']};zoneparts=[];grid=[];foot=[];roleparts=[];pairs=[]
def aggregate(q):
 # Equal visit means within transects, equal transects within zones, equal zones within stations.
 return q.groupby(['site','status','year','zone','transect','Taxa','Resp'],as_index=False).value.agg(lambda x:x.mean(skipna=False)).groupby(['site','status','year','zone','Taxa','Resp'],as_index=False).value.agg(lambda x:x.mean(skipna=False))
for season,(lo,hi) in windows.items():
 z=aggregate(g[g.month.between(lo,hi)]);z['season']=season;zoneparts.append(z)
 for yn,years in yearsets.items():
  for zn,zones in zone_sets.items():
   q=z[z.year.isin(years)&z.zone.isin(zones)].copy();complete=q.groupby(['site','year']).filter(lambda x:len(x)==len(zones)*4 and x.value.notna().all());sites=complete.groupby('site').filter(lambda x:set(x.year)==set(years)).site.unique();q=q[q.site.isin(sites)];roles=set(q.status);ok=roles=={'mpa','reference'};name=f'{season}_{zn}_{yn}';grid.append(dict(scenario=name,season=season,zones=zn,year_window=yn,eligible=ok,years=';'.join(map(str,years)),protected=';'.join(sorted(q[q.status=='mpa'].site.unique())),reference=';'.join(sorted(q[q.status=='reference'].site.unique())),n_stations=len(sites)))
   if not ok:continue
   q['scenario']=name;foot.append(q)
   rr=q.groupby(['site','status','year','Taxa','Resp'],as_index=False).value.mean().groupby(['status','year','Taxa','Resp'],as_index=False).value.mean();rr['scenario']=name;roleparts.append(rr)
   timing=coverage[coverage.site.isin(sites)&coverage.year.isin(years)&coverage.zone.isin(zones)&coverage.month.between(lo,hi)]
   for y,yy in timing.groupby('year'):
    pp=yy[yy.status=='mpa'];re_=yy[yy.status=='reference'];pmon=sorted(pp.month.unique());rmon=sorted(re_.month.unique());pdates=pd.to_datetime(pp[['year','month','day']]);rdates=pd.to_datetime(re_[['year','month','day']]);gap=np.abs((pdates.to_numpy()[:,None]-rdates.to_numpy()[None,:])/np.timedelta64(1,'D'))
    pairs.append(dict(scenario=name,year=y,protected_months=';'.join(map(str,pmon)),reference_months=';'.join(map(str,rmon)),same_month_possible=bool(set(pmon)&set(rmon)),minimum_month_gap=min(abs(x-y) for x in pmon for y in rmon),minimum_day_gap=float(gap.min()),maximum_day_gap=float(gap.max()),protected_depth_min=pp.depth_min.min(),protected_depth_max=pp.depth_max.max(),reference_depth_min=re_.depth_min.min(),reference_depth_max=re_.depth_max.max()))
pd.concat(zoneparts).to_csv(O/'zone_values.csv',index=False);pd.DataFrame(grid).to_csv(O/'coverage_grid.csv',index=False);pd.concat(foot).to_csv(O/'grid_footprints.csv',index=False);pd.concat(roleparts).to_csv(O/'grid_role_values.csv',index=False);pd.DataFrame(pairs).to_csv(O/'date_alignment.csv',index=False)
# Validate original May-Oct available-zone annual values on the original years before changing support.
raw=rd(Path('data/harmonized/harmonized_raw_responses.csv')).rename(columns={'CA_MPA_Name_Short':'MPA','source':'Source','resp':'Resp'});nm={'Macrocystis pyrifera':'M. pyrifera','Mesocentrotus franciscanus':'M. franciscanus','Strongylocentrotus purpuratus':'S. purpuratus','Panulirus interruptus':'P. interruptus'};raw['Taxa']=raw.taxon_name.str.strip().map(nm);legacy=g[g.month.between(5,10)].copy();legacy.value=legacy.value.fillna(0);lz=aggregate(legacy);lr=lz.groupby(['site','status','year','Taxa','Resp'],as_index=False).value.mean().groupby(['status','year','Taxa','Resp'],as_index=False).value.mean();checks=[];offsets={}
legacy_annual=rd(Path('outputs/researcher_scope_20261101/annual_input.csv'));legacy_annual=legacy_annual[(legacy_annual.MPA==m)&(legacy_annual.Source=='PISCO')]
for tx,resp in taxa.values():
 h=legacy_annual[(legacy_annual.Taxa==tx)&(legacy_annual.Resp==resp)&(legacy_annual.t>0)].sort_values('year');r=raw[(raw.MPA==m)&(raw.Source=='PISCO')&(raw.Taxa==tx)&(raw.Resp==resp)];j=lr[(lr.Taxa==tx)&(lr.Resp==resp)].merge(r,on=['year','status','Taxa','Resp'],suffixes=('_new','_old'));assert (j.value_new-j.value_old).abs().max()<1e-7
 v=lr[(lr.Taxa==tx)&(lr.Resp==resp)].pivot(index='year',columns='status',values='value').reindex(h.year);rat=np.exp(h.lnDiff.to_numpy());imp=np.divide(v.mpa.to_numpy()-rat*v.reference.to_numpy(),rat-1,out=np.full(len(rat),np.nan),where=abs(rat-1)>1e-6);offset=float(np.nanmedian(imp));er=max(abs(np.log((v.mpa+offset)/(v.reference+offset)).to_numpy()-h.lnDiff.to_numpy()));assert er<1e-7;offsets[(tx,resp)]=offset;checks.append(dict(Taxa=tx,Resp=resp,offset=offset,max_raw_error=(j.value_new-j.value_old).abs().max(),max_ratio_error=er,n_raw=len(j)))
pd.DataFrame(checks).to_csv(O/'reconstruction.csv',index=False)
# Fixed Apr-Oct full11 footprints anchor all main contrasts; do not reselect controls for a time-window test.
rrparts=[];modelparts=[];fparts=[]
allz=pd.concat(zoneparts)
for zn in zone_sets:
 anchor=f'AprOct_{zn}_full11';eligible=pd.DataFrame(grid).set_index('scenario').loc[anchor,'eligible']
 if not eligible:continue
 ff=pd.concat(foot);ff=ff[ff.scenario==anchor][['site','status','zone']].drop_duplicates()
 for season,yn in [('AprOct','full11'),('AprOct','old10'),('AprOct','late9'),('MayOct','late9'),('AprNov','full11'),('MarNov','full11'),('AllMonths','full11')]:
  sc=f'{season}_{zn}_{yn}_anchored';z=allz[(allz.season==season)&allz.year.isin(yearsets[yn])].merge(ff,on=['site','status','zone']);assert z.value.notna().all();assert len(z)==len(ff)*len(yearsets[yn])*4;z['scenario']=sc;fparts.append(z)
  roles=z.groupby(['site','status','year','Taxa','Resp'],as_index=False).value.mean().groupby(['status','year','Taxa','Resp'],as_index=False).value.mean();roles['scenario']=sc;rrparts.append(roles)
  for tx,resp in taxa.values():
   h=old[(old.Taxa==tx)&(old.Resp==resp)];v=roles[(roles.Taxa==tx)&(roles.Resp==resp)].pivot(index='year',columns='status',values='value');q=pd.DataFrame([h.iloc[0].to_dict() for _ in v.index]);q['year']=v.index.to_numpy();q['t']=q.year-q.MPA_Start;q['time']=q.t;q['BA']='After';q['lnDiff']=np.log((v.mpa.to_numpy()+offsets[(tx,resp)])/(v.reference.to_numpy()+offsets[(tx,resp)]));q['scenario']=sc;q['baseline_design_verified']=False;modelparts.append(q)
pd.concat(fparts).to_csv(O/'anchored_footprints.csv',index=False);pd.concat(rrparts).to_csv(O/'anchored_role_values.csv',index=False);mpatch=pd.concat(modelparts);mpatch.drop(columns=['mpa','reference','Diff'],errors='ignore').to_csv(O/'local_annual.csv',index=False)
# Rerun pooling for the full11 and old10 anchored variants only; late9 fits are explicitly local time-window diagnostics.
parts=[a.assign(scenario='Previous_review')]
for sc in mpatch.scenario.unique():
 if 'late9' in sc:continue
 q=a.copy();patch=mpatch[mpatch.scenario==sc];ids=set(patch.series_id);q=q[~(q.series_id.isin(ids)&(q.t>0))];q=pd.concat([q,patch],ignore_index=True);q.scenario=sc;parts.append(q)
pd.concat(parts).drop(columns=['mpa','reference','Diff'],errors='ignore').to_csv(O/'annual_input.csv',index=False)
(O/'source_hashes.json').write_text(json.dumps({p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in used},indent=2));print(pd.DataFrame(grid).query("season=='AprOct'").to_string(index=False));print(pd.DataFrame(pairs).query("scenario=='AprOct_INNER_MID_full11'").to_string(index=False))
