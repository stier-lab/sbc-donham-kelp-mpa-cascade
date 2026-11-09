"""Isolated, source-based seasonal and timing coverage audit. See output PLAN.md."""
from pathlib import Path
import os,re,json,hashlib
import pandas as pd
import numpy as np
O=Path('outputs/cross_site_seasons_20261101')
D=Path(os.environ.get('DONHAM_DATA_REPO_PISCO',str(Path.home()/'Donham-Stier-CA-MPA-Data-2026/data/PISCO')))
used=[]
def rd(p,**kw):
 used.append(str(p));return pd.read_csv(p,low_memory=False,**kw)
a=rd(Path('outputs/pisco_source_reconciliation_20261101/annual_input.csv'));a=a[a.scenario=='Supported_corrections'].copy()
mpas=['Abalone Cove SMCA','Gull Island SMR','Point Dume SMCA','Point Vicente SMCA','South Point SMR']
taxa={'MACPYRAD':('M. pyrifera','Bio'),'MESFRAAD':('M. franciscanus','Den'),'STRPURAD':('S. purpuratus','Den'),'PANINT':('P. interruptus','Den')}
sw=rd(D/'MLPA_kelpforest_swath_2024.csv',encoding='latin1');fi=rd(D/'MLPA_kelpforest_fish_2024.csv',encoding='latin1')
sw=sw[(sw.campus=='VRG')|((sw.campus=='UCSB')&~sw.zone.isin(['DEEP','SPECIALO']))].copy();fi=fi[fi.campus.isin(['VRG','UCSB'])&(fi.level!='CAN')]
master=rd(D/'master_site_table_Emilyedit.csv',encoding='utf-8-sig');master=master[(master.BaselineRegion=='SOUTH')&master.CA_MPA_Name_Short.isin(mpas)].rename(columns={'CA_MPA_Name_Short':'MPA','site_status':'status'})[['site','MPA','status']].drop_duplicates()
assert not master.duplicated(['site','MPA']).any()
u=D.parent.parent/'code/R/01_utils.R';used.append(str(u));exc=re.findall(r'"([^"]+)"',re.search(r'EXCLUDED_REFERENCE_SITES <- c\((.*?)\)',u.read_text(),re.S).group(1))
valid=pd.concat([sw[sw.month.between(5,10)][['site','year']],fi[['site','year']]]).drop_duplicates().groupby('site').size()
master=master[master.site.isin(valid[valid>=5].index)&~master.site.isin(exc)&master.status.isin(['mpa','reference'])];master.to_csv(O/'candidate_stations.csv',index=False)
sw=sw[sw.site.isin(master.site)&(sw.year<=2023)].copy()
ev=['site','year','month','day','zone','transect'];events=sw[ev].drop_duplicates()
# Do not collapse distinct dated visits into one transect observation.
f=sw[sw.classcode.isin(taxa)].copy();f['value']=f['count'];k=f.classcode=='MACPYRAD';f.loc[k,'value']=f.loc[k,'count']*f.loc[k,'size']
f[f.value.isna()][ev+['classcode','count','size']].to_csv(O/'missing_measurements.csv',index=False)
ct=f.groupby(ev+['classcode'],as_index=False).value.agg(lambda x:x.sum(skipna=False));ct['present']=True
g=events.merge(pd.DataFrame({'classcode':list(taxa)}),how='cross').merge(ct,on=ev+['classcode'],how='left');g.loc[g.present.isna(),'value']=0;g['Taxa']=g.classcode.map(lambda x:taxa[x][0]);g['Resp']=g.classcode.map(lambda x:taxa[x][1]);g.value/=60
co=rd(Path('data/external/sbc1275/Macrocystis_Biomass_Relationships_20261016.csv'));co=co[co.model==1].set_index('MONTH').slope;g.loc[g.Resp=='Bio','value']*=g.loc[g.Resp=='Bio','month'].map(co)*1000
g.to_csv(O/'event_responses.csv',index=False)
dates=sw.groupby(['site','year','month','day','zone'],as_index=False).agg(depth_min=('depth','min'),depth_max=('depth','max'),n_transects=('transect','nunique'));dates.to_csv(O/'survey_dates_depths.csv',index=False)
def aggregate(q):
 return q.groupby(['site','year','zone','transect','Taxa','Resp'],as_index=False).value.agg(lambda x:x.mean(skipna=False)).groupby(['site','year','zone','Taxa','Resp'],as_index=False).value.agg(lambda x:x.mean(skipna=False))
windows={'MayOct':(5,10),'AprOct':(4,10),'AprNov':(4,11),'MarNov':(3,11),'AllMonths':(1,12)};zonesets={'IM':['INNER','MID'],'IMO':['INNER','MID','OUTER']}
zone={s:aggregate(g[g.month.between(lo,hi)]) for s,(lo,hi) in windows.items()}
pd.concat([z.assign(season=s) for s,z in zone.items()]).to_csv(O/'zone_values.csv',index=False)
# Reproduce legacy role values and infer frozen offsets from original paired histories.
orig=rd(Path('outputs/researcher_scope_20261101/annual_input.csv'));recon=[];offsets={}
legacy=aggregate(g[g.month.between(5,10)].fillna({'value':0})).merge(master,on='site')
raw=rd(Path('data/harmonized/harmonized_raw_responses.csv')).rename(columns={'CA_MPA_Name_Short':'MPA','source':'Source','resp':'Resp'})
nm={'Macrocystis pyrifera':'M. pyrifera','Mesocentrotus franciscanus':'M. franciscanus','Strongylocentrotus purpuratus':'S. purpuratus','Panulirus interruptus':'P. interruptus'};raw['Taxa']=raw.taxon_name.str.strip().map(nm)
new=rd(Path('outputs/upstream_sites_20261101/unpaired_annual_values.csv'));raw=pd.concat([raw[['MPA','Source','Taxa','Resp','year','status','value']],new[['MPA','Source','Taxa','Resp','year','status','value']]]).drop_duplicates(['MPA','Source','Taxa','Resp','year','status'],keep='last')
for m in mpas:
 lr=legacy[legacy.MPA==m].groupby(['site','status','year','Taxa','Resp'],as_index=False).value.mean().groupby(['status','year','Taxa','Resp'],as_index=False).value.mean()
 for tx,resp in taxa.values():
  h=orig[(orig.MPA==m)&(orig.Source=='PISCO')&(orig.Taxa==tx)&(orig.Resp==resp)&(orig.t>0)].sort_values('year');v=lr[(lr.Taxa==tx)&(lr.Resp==resp)].pivot(index='year',columns='status',values='value').reindex(h.year)
  rat=np.exp(h.lnDiff.to_numpy());im=np.divide(v.mpa.to_numpy()-rat*v.reference.to_numpy(),rat-1,out=np.full(len(rat),np.nan),where=abs(rat-1)>1e-6);off=float(np.nanmedian(im));err=max(abs(np.log((v.mpa+off)/(v.reference+off)).to_numpy()-h.lnDiff.to_numpy()));assert err<1e-7,(m,tx,err)
  r=raw[(raw.MPA==m)&(raw.Source=='PISCO')&(raw.Taxa==tx)&(raw.Resp==resp)&raw.year.isin(h.year)];j=lr[(lr.Taxa==tx)&(lr.Resp==resp)].merge(r,on=['status','year','Taxa','Resp'],suffixes=('_new','_old'));er=(j.value_new-j.value_old).abs().max();assert er<1e-7,(m,tx,er)
  offsets[(m,tx,resp)]=off;recon.append(dict(MPA=m,Taxa=tx,Resp=resp,offset=off,n=len(h),raw_max_error=er,ratio_max_error=err))
pd.DataFrame(recon).to_csv(O/'reconstruction.csv',index=False)
# Coverage only: select best pair before any fitted model or outcome comparison.
ledger=[];rules=[];selected=[];pairyears=[];foot=[];roles=[];model=[]
for m in mpas:
 mm=master[master.MPA==m];ha=a[(a.MPA==m)&(a.Source=='PISCO')];start=int(ha.MPA_Start.iloc[0]);md=dates.merge(mm,on='site');post=md[md.year>start]
 target=sorted(set(post[post.status=='mpa'].year)&set(post[post.status=='reference'].year));oldyears=sorted(ha[(ha.Taxa=='M. pyrifera')&(ha.Resp=='Bio')&(ha.t>0)].year)
 for season,(lo,hi) in windows.items():
  zz=zone[season];cc=dates[dates.month.between(lo,hi)].copy();cc['date']=pd.to_datetime(cc[['year','month','day']]);lookup={key:q.date.to_numpy() for key,q in cc.groupby(['site','year','zone'])}
  for zn,zones in zonesets.items():
   z=zz[zz.year.isin(target)&zz.zone.isin(zones)&zz.site.isin(mm.site)]
   validrows=z.groupby(['site','year']).filter(lambda x:len(x)==4*len(zones) and x.value.notna().all());keys=set(map(tuple,validrows[['site','year']].drop_duplicates().to_numpy()))
   candidates=[]
   for psite in sorted(mm[mm.status=='mpa'].site):
    for rsite in sorted(mm[mm.status=='reference'].site):
     gap={}
     for year in target:
      validp=(psite,year) in keys;validr=(rsite,year) in keys
      if validp and validr:
       gaps=[]
       for zname in zones:
        pd_=lookup[(psite,year,zname)];rd_=lookup[(rsite,year,zname)];gaps.extend(np.abs((pd_[:,None]-rd_[None,:])/np.timedelta64(1,'D')).ravel())
       gap[year]=float(max(gaps))
      pairyears.append(dict(MPA=m,season=season,zones=zn,protected=psite,reference=rsite,year=year,protected_complete=validp,reference_complete=validr,max_gap_days=gap.get(year,np.nan)))
     for limit in [0,30,60,90]:
      ys=sorted(y for y,d in gap.items() if limit==0 or d<=limit)
      row=dict(MPA=m,season=season,zones=zn,limit=limit,protected=psite,reference=rsite,n=len(ys),years=';'.join(map(str,ys)),span=ys[-1]-ys[0] if ys else 0,first=ys[0] if ys else 9999,full_history=set(ys)==set(target),target_years=';'.join(map(str,target)))
      ledger.append(row);candidates.append(row)
   for limit in [0,30,60,90]:
    cc_=sorted([x for x in candidates if x['limit']==limit],key=lambda x:(-x['n'],-x['span'],x['first'],x['protected'],x['reference']));best=cc_[0];ys=[int(y) for y in best['years'].split(';') if y];sc=f'{season}_{zn}_d{limit}';full=sum(x['full_history'] for x in cc_)
    rules.append(dict(best,scenario=sc,n_full_pairs=full,n_candidate_pairs=len(cc_),original_years=';'.join(map(str,oldyears)),lost_original=';'.join(map(str,sorted(set(oldyears)-set(ys)))),added_years=';'.join(map(str,sorted(set(ys)-set(oldyears)))),eligible_slope=len(ys)>=3))
    if len(ys)<3:continue
    fz=z[z.site.isin([best['protected'],best['reference']])&z.year.isin(ys)].merge(mm,on='site');assert len(fz)==2*len(zones)*len(ys)*4;fz['scenario']=sc;foot.append(fz)
    rr=fz.groupby(['MPA','status','year','Taxa','Resp'],as_index=False).value.mean();rr['scenario']=sc;roles.append(rr)
    for tx,resp in taxa.values():
     h=ha[(ha.Taxa==tx)&(ha.Resp==resp)];v=rr[(rr.Taxa==tx)&(rr.Resp==resp)].pivot(index='year',columns='status',values='value');q=pd.DataFrame([h.iloc[0].to_dict() for y in v.index]);q['year']=v.index.to_numpy();q['t']=q.year-q.MPA_Start;q['time']=q.t;q['BA']='After';q['baseline_design_verified']=False;q['lnDiff']=np.log((v.mpa.to_numpy()+offsets[(m,tx,resp)])/(v.reference.to_numpy()+offsets[(m,tx,resp)]));q['scenario']=sc;model.append(q)
# Anchor an additional six timing contrasts to one coverage-chosen pair per MPA.
# This prevents reference reselection from masquerading as an effect of the month window.
anchors=[dict(r) for r in rules if r['scenario']=='AprNov_IM_d90']
for anchor in anchors:
 m=anchor['MPA'];mm=master[master.MPA==m];ha=a[(a.MPA==m)&(a.Source=='PISCO')]
 oldyears=set(int(y) for y in anchor['original_years'].split(';'))
 for season in ['AprOct','AprNov']:
  for limit in [0,60,90]:
   r=next(dict(v) for v in ledger if v['MPA']==m and v['season']==season and v['zones']=='IM' and v['limit']==limit and v['protected']==anchor['protected'] and v['reference']==anchor['reference'])
   ys=[int(y) for y in r['years'].split(';') if y];sc=f'Anchor_{season}_IM_d{limit}'
   rules.append(dict(r,scenario=sc,n_full_pairs=int(r['full_history']),n_candidate_pairs=1,original_years=anchor['original_years'],lost_original=';'.join(map(str,sorted(oldyears-set(ys)))),added_years=';'.join(map(str,sorted(set(ys)-oldyears))),eligible_slope=len(ys)>=3))
   if len(ys)<3:continue
   fz=zone[season];fz=fz[fz.site.isin([r['protected'],r['reference']])&fz.year.isin(ys)&fz.zone.isin(['INNER','MID'])].merge(mm,on='site');assert fz.value.notna().all() and len(fz)==2*2*len(ys)*4
   fz['scenario']=sc;foot.append(fz);rr=fz.groupby(['MPA','status','year','Taxa','Resp'],as_index=False).value.mean();rr['scenario']=sc;roles.append(rr)
   for tx,resp in taxa.values():
    h=ha[(ha.Taxa==tx)&(ha.Resp==resp)];v=rr[(rr.Taxa==tx)&(rr.Resp==resp)].pivot(index='year',columns='status',values='value');q=pd.DataFrame([h.iloc[0].to_dict() for y in v.index]);q['year']=v.index.to_numpy();q['t']=q.year-q.MPA_Start;q['time']=q.t;q['BA']='After';q['baseline_design_verified']=False;q['lnDiff']=np.log((v.mpa.to_numpy()+offsets[(m,tx,resp)])/(v.reference.to_numpy()+offsets[(m,tx,resp)]));q['scenario']=sc;model.append(q)
pd.DataFrame(ledger).to_csv(O/'pair_coverage.csv',index=False);pd.DataFrame(pairyears).to_csv(O/'pair_year_gaps.csv',index=False);pd.DataFrame(rules).to_csv(O/'selected_coverage.csv',index=False);pd.concat(foot).to_csv(O/'selected_footprints.csv',index=False);pd.concat(roles).to_csv(O/'role_values.csv',index=False);pd.concat(model).drop(columns=['mpa','reference','Diff'],errors='ignore').to_csv(O/'local_annual.csv',index=False)
(O/'source_hashes.json').write_text(json.dumps({p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in used},indent=2))
print(pd.DataFrame(rules).query("season in ['AprOct','AprNov'] and zones=='IM'")[['MPA','scenario','protected','reference','years','n_full_pairs','lost_original','added_years']].to_string(index=False))
