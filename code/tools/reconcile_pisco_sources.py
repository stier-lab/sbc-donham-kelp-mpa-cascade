from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,re
O=Path('outputs/pisco_source_reconciliation_20261101');D=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data/PISCO');oldout=Path('outputs/seasonal_zone_audit_20261101');used=[]
def rd(p,**kw):used.append(str(p));return pd.read_csv(p,low_memory=False,**kw)
a=rd(oldout/'annual_input.csv');a=a[a.scenario=='Fixed_all_program_units'].copy();targets=[('Point Vicente SMCA','M. pyrifera','Bio','MACPYRAD'),('Scorpion SMR','M. franciscanus','Den','MESFRAAD'),('South Point SMR','S. pulcher','Bio','SPUL')]
sw=rd(D/'MLPA_kelpforest_swath_2024.csv',encoding='latin1');fi=rd(D/'MLPA_kelpforest_fish_2024.csv',encoding='latin1');sw=sw[(sw.campus=='VRG')|((sw.campus=='UCSB')&~sw.zone.isin(['DEEP','SPECIALO']))];fi=fi[fi.campus.isin(['VRG','UCSB'])&(fi.level!='CAN')]
mp=rd(D/'master_site_table_Emilyedit.csv',encoding='utf-8-sig');mp=mp[mp.BaselineRegion=='SOUTH'].rename(columns={'CA_MPA_Name_Short':'MPA','site_status':'status'})[['site','MPA','status']].drop_duplicates(['site','MPA']);valid=pd.concat([sw[sw.month.between(5,10)][['site','year']],fi[['site','year']]]).drop_duplicates().groupby('site').size();utils=D.parent.parent/'code/R/01_utils.R';used.append(str(utils));exc=re.findall(r'"([^"]+)"',re.search(r'EXCLUDED_REFERENCE_SITES <- c\((.*?)\)',utils.read_text(),re.S).group(1));mp=mp[mp.site.isin(valid[valid>=5].index)&~mp.site.isin(exc)&mp.status.isin(['mpa','reference'])]
raw=rd(Path('data/harmonized/harmonized_raw_responses.csv')).rename(columns={'CA_MPA_Name_Short':'MPA','source':'Source','resp':'Resp'});nm={'Macrocystis pyrifera':'M. pyrifera','Mesocentrotus franciscanus':'M. franciscanus','Semicossyphus pulcher':'S. pulcher'};raw['Taxa']=raw.taxon_name.str.strip().map(nm)
co=rd(Path('data/external/sbc1275/Macrocystis_Biomass_Relationships_20261016.csv'));co=co[co.model==1].set_index('MONTH').slope
checks=[];anomalies=[];zones=[];patches=[];foot=[];availability=[];allrole=[]
for m,tx,resp,code in targets:
 h=a[(a.MPA==m)&(a.Source=='PISCO')&(a.Taxa==tx)&(a.Resp==resp)].sort_values('year').copy();sid=h.series_id.iloc[0];hasbase=h.baseline_design_verified.all() and (h.t<0).sum()>=3;h=h[(h.t!=0)&((h.t>0)|hasbase)];map_=mp[mp.MPA==m];data=(fi if code=='SPUL' else sw[sw.month.between(5,10)]);data=data[data.site.isin(map_.site)&data.year.isin(h.year)].copy();f=data[data.classcode==code].copy()
 if code=='SPUL':
  f['amount']=f['count']*.0144*f.fish_tl**3.04;bad=f[(f['count']>0)&f.amount.isna()].copy();keys=['site','year','zone','transect'];ev=data[keys].drop_duplicates();ct=f.groupby(keys,as_index=False).amount.agg(lambda v: v.sum(skipna=False));grid=ev.merge(ct,on=keys,how='left');grid.amount=grid.amount.fillna(0);grid['value']=grid.amount/60;zz=grid.groupby(['site','year','zone'],as_index=False).value.mean()
 else:
  f['amount']=f['count']*f['size'] if code=='MACPYRAD' else f['count'];bad=f[f.amount.isna()].copy();keys=['site','year','month','zone','transect'];ev=data[keys].drop_duplicates();f.amount=f.amount.fillna(0);ct=f.groupby(keys,as_index=False).amount.sum();grid=ev.merge(ct,on=keys,how='left');grid.amount=grid.amount.fillna(0);grid['value']=grid.amount/60
  if code=='MACPYRAD':grid.value*=grid.month.map(co)*1000
  zz=grid.groupby(['site','year','zone','transect'],as_index=False).value.mean().groupby(['site','year','zone'],as_index=False).value.mean()
 assert len(bad)==1
 # Minimal source evidence excludes observer names, retaining physical keys and missing fields.
 for r in bad.itertuples():anomalies.append(dict(series_id=sid,site=r.site,year=r.year,month=r.month,day=r.day,zone=r.zone,transect=r.transect,count=r.count,size=getattr(r,'size',np.nan),fish_tl=getattr(r,'fish_tl',np.nan),notes=str(r.notes),classcode=code))
 zz=zz.merge(map_,on='site',validate='many_to_one');zz['series_id']=sid;zz['valid']=True
 for r in bad.itertuples():zz.loc[(zz.site==r.site)&(zz.year==r.year)&(zz.zone==r.zone),'valid']=False
 zones.append(zz)
 def ag(z):return z.groupby(['site','status','year'],as_index=False).value.mean().groupby(['status','year'],as_index=False).value.mean()
 orig=ag(zz);rr=raw[(raw.MPA==m)&(raw.Source=='PISCO')&(raw.Taxa==tx)&(raw.Resp==resp)];j=orig.merge(rr,on=['year','status'],suffixes=('_rebuilt','_old'));err=(j.value_rebuilt-j.value_old).abs().max();assert len(j)==2*len(h) and err<1e-7,(sid,err)
 vv=orig.pivot(index='year',columns='status',values='value').reindex(h.year);R=np.exp(h.lnDiff.to_numpy());implied=np.divide(vv.mpa.to_numpy()-R*vv.reference.to_numpy(),R-1,out=np.full(len(R),np.nan),where=abs(R-1)>1e-6);offset=float(np.nanmedian(implied));re=max(abs(np.log((vv.mpa+offset)/(vv.reference+offset)).to_numpy()-h.lnDiff.to_numpy()));assert re<1e-7
 checks.append(dict(series_id=sid,raw_rows=len(j),max_raw_error=err,max_ratio_error=re,offset=offset,affected_year=int(bad.year.iloc[0])))
 for scenario in ['Legacy_reproduced','Complete_zones','Fixed_valid_zones']:
  z=zz.copy() if scenario=='Legacy_reproduced' else zz[zz.valid].copy()
  if scenario=='Fixed_valid_zones':z=z.groupby(['site','zone']).filter(lambda q:set(q.year)==set(h.year))
  v=ag(z).pivot(index='year',columns='status',values='value').reindex(index=h.year,columns=['mpa','reference']);complete=v.notna().all().all();availability.append(dict(series_id=sid,scenario=scenario,complete=complete,station_zones=z[['site','zone']].drop_duplicates().shape[0],sites=z.site.nunique()))
  if not complete:continue
  hh=h.copy();hh['lnDiff']=np.log((v.mpa.to_numpy()+offset)/(v.reference.to_numpy()+offset));hh['scenario']=scenario;patches.append(hh);z['scenario']=scenario;foot.append(z);role=ag(z);role['scenario']=scenario;role['series_id']=sid;allrole.append(role)
pd.DataFrame(checks).to_csv(O/'reconstruction.csv',index=False);pd.DataFrame(anomalies).to_csv(O/'source_anomalies.csv',index=False);pd.concat(zones).to_csv(O/'legacy_zone_values.csv',index=False);pd.concat(foot).to_csv(O/'scenario_zone_values.csv',index=False);pd.concat(allrole).to_csv(O/'role_values.csv',index=False);pd.DataFrame(availability).to_csv(O/'coverage.csv',index=False);p=pd.concat(patches);p.to_csv(O/'patches.csv',index=False)
parts=[]
for s in ['Previous_review','Complete_zones','Fixed_valid_zones']:
 z=a.copy();z.scenario=s
 if s!='Previous_review':
  pp=p[p.scenario=='Complete_zones'].set_index(['series_id','year'])
  if s=='Fixed_valid_zones':
   fixedp=p[p.scenario==s].set_index(['series_id','year']);pp.update(fixedp)
  ix=pd.MultiIndex.from_frame(z[['series_id','year']]);mask=ix.isin(pp.index);z.loc[mask,'lnDiff']=pp.loc[ix[mask],'lnDiff'].to_numpy()
 parts.append(z)
supported=parts[-1].copy();supported.scenario='Supported_corrections'
# Do not promote the mechanically continuous DEEP-only Point Vicente subset.
mask=(supported.MPA=='Point Vicente SMCA')&(supported.Taxa=='M. pyrifera')&(supported.Source=='PISCO')&(supported.Resp=='Bio')
supported.loc[mask,'lnDiff']=parts[1].loc[mask,'lnDiff'].to_numpy();parts.append(supported)
pd.concat(parts).to_csv(O/'annual_input.csv',index=False)
# Independent physical-depth evidence for the invalid DEEP-only interpretation.
depth=sw[(sw.site=='POINT_VICENTE_W')&(sw.year==2013)].groupby(['site','year','month','day','zone']).agg(transects=('transect','nunique'),depth_min=('depth','min'),depth_max=('depth','max')).reset_index();depth['season_retained']=depth.month.between(5,10);depth.to_csv(O/'point_vicente_2013_depth_audit.csv',index=False)

# Verify source copies; they cannot independently recover the absent measurements.
archives=[]
for n in ['MLPA_kelpforest_swath_2024.csv','MLPA_kelpforest_fish_2024.csv']:
 source=D/n;backup=source.resolve().parents[2]/'Github/MPA_Kelp/data/PISCO'/n;used.extend([str(source),str(backup)]);archives.append(dict(file=n,archive=str(backup),identical=hashlib.sha256(source.read_bytes()).digest()==hashlib.sha256(backup.read_bytes()).digest()))
pd.DataFrame(archives).to_csv(O/'archive_comparison.csv',index=False)
(O/'source_hashes.json').write_text(json.dumps({p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in used},indent=2));print(pd.DataFrame(checks).to_string(index=False));print(pd.DataFrame(availability).to_string(index=False))
