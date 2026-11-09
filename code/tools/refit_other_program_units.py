# Survey-unit reconstruction and conservative fixed-unit patch; no canonical writes.
from pathlib import Path
import pandas as pd,numpy as np,re,json,hashlib
O=Path('outputs/seasonal_zone_audit_20261101');D=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data');a=pd.read_csv('outputs/researcher_scope_20261101/annual_input.csv');used=[]
def read(p):
 used.append(str(p));return pd.read_csv(p,low_memory=False)
mm=read(D/'MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv').rename(columns={'site_id':'site','CA_MPA_Name_Short':'MPA'});mm['Source']=mm.data_source.str.upper();text=(D.parent/'code/R/01_utils.R').read_text();exc=[]
for name in ['EXCLUDED_KFM_SITES','EXCLUDED_LTER_SITES','EXCLUDED_LTER_SITE_CODES']:exc+=re.findall(r'"([^\"]+)"',re.search(name+r' <- c\((.*?)\)',text,re.S).group(1))
fishmm=mm[mm.status.isin(['mpa','reference'])][['Source','MPA','site','status']].drop_duplicates()
mm=mm[~mm.site.isin(exc)&mm.status.isin(['mpa','reference'])][['Source','MPA','site','status']].drop_duplicates()
names={'Strongylocentrotus purpuratus':'S. purpuratus','Mesocentrotus franciscanus':'M. franciscanus','Panulirus interruptus':'P. interruptus','Bodianus pulcher':'S. pulcher'}
x=read(D/'MBON/SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv');x['Source']=x.data_source.str.upper();x['Taxa']=x.taxon_name.map(names);x['site']=x.site_id;x['year']=pd.to_datetime(x.date).dt.year;x=x.merge(mm,on=['Source','site']);x=x[x.Taxa.notna()];x['count']=pd.to_numeric(x['count'],errors='coerce');assert x['count'].notna().all();parts=[]
# KFM: sum size-class entries per actual replicate, then average replicate density.
k=x[x.Source=='KFM'];k=k.groupby(['MPA','site','status','year','Taxa','area','replicate_id'],as_index=False)['count'].sum();k['value']=k['count']/k.area;k=k.groupby(['MPA','site','status','year','Taxa','area'],as_index=False).value.mean()
# Area is a method distinction; each focal taxon must have a single area per site-year.
assert not k.duplicated(['MPA','site','status','year','Taxa']).any();k['unit']='station';k['Resp']='Den';k['Source']='KFM';parts.append(k)
l=x[(x.Source=='LTER')&(x.sample_method=='quad')&x.Taxa.isin(['S. purpuratus','M. franciscanus'])].copy();l['value']=l['count']/l.area;l=l.groupby(['MPA','site','status','year','Taxa','transect_id'],as_index=False).value.mean().rename(columns={'transect_id':'unit'});l['Resp']='Den';l['Source']='LTER';parts.append(l)
f=read(D/'MBON/SBCMBON_kelp_forest_integrated_fish_20231022.csv');f=f[(f.data_source=='kfm')&f.proj_taxon_id.isin(['t-k-127','t-k-233'])&(f.sample_method=='visualfish')].copy();f['Source']='KFM';f['site']=f.site_id;f['year']=pd.to_datetime(f.date).dt.year;f['count']=pd.to_numeric(f['count'],errors='coerce');assert f['count'].notna().all();f['value']=f['count']/f.area;f=f.merge(fishmm,on=['Source','site']);f=f.groupby(['MPA','site','status','year','transect_id'],as_index=False).value.mean().rename(columns={'transect_id':'unit'});f['Source']='KFM';f['Taxa']='S. pulcher';f['Resp']='Den';parts.append(f)
f=read(D/'LTER/Annual_fish_comb_20240307.csv').rename(columns={'SITE':'site','YEAR':'year','TRANSECT':'unit'});f=f[(f.SCIENTIFIC_NAME=='Semicossyphus pulcher')&(f.COUNT!=-99999)].copy();f['Source']='LTER';f=f.merge(mm,on=['Source','site']);f['Den']=f.COUNT/80;f['Bio']=f.COUNT*.0144*f.SIZE**3.04/80;f=f.groupby(['MPA','site','status','year','unit'],as_index=False)[['Den','Bio']].sum();f=f.melt(id_vars=['MPA','site','status','year','unit'],value_vars=['Den','Bio'],var_name='Resp',value_name='value');f['Source']='LTER';f['Taxa']='S. pulcher';parts.append(f)
c=read(Path('data/external/sbc1275/Macrocystis_Biomass_Relationships_20261016.csv'));c=c[c.model==1].set_index('MONTH').slope
l=read(D/'LTER/Annual_Kelp_All_Years_20240305.csv').rename(columns={'SITE':'site','YEAR':'year','TRANSECT':'unit'});l=l[l.FRONDS!=-99999].copy();l['Source']='LTER';l=l.merge(mm,on=['Source','site']);l['value']=l.FRONDS/l.AREA*l.MONTH.map(c)*1000;l=l.groupby(['MPA','site','status','year','unit','QUAD','SIDE','AREA'],as_index=False).value.sum().groupby(['MPA','site','status','year','unit'],as_index=False).value.mean();l['Source']='LTER';l['Taxa']='M. pyrifera';l['Resp']='Bio';parts.append(l)
units=pd.concat(parts)[['MPA','site','status','year','unit','Source','Taxa','Resp','value']];units['unit']=units.unit.astype(str);units.to_csv(O/'other_program_unit_values.csv',index=False)
raw=read(Path('data/harmonized/harmonized_raw_responses.csv')).rename(columns={'CA_MPA_Name_Short':'MPA','source':'Source','resp':'Resp'});nm={**names,'Semicossyphus pulcher':'S. pulcher','Macrocystis pyrifera':'M. pyrifera'};raw['Taxa']=raw.taxon_name.str.strip().map(nm)
# Restored histories use original unpaired annual values from the diver review.
for path in [Path('outputs/diver_exclusion_refit_20261101/consistent_kfm/unpaired_annual_values.csv')]:
 if path.exists():
  rr=read(path);raw=pd.concat([raw,rr])
checks=[];patches=[];foot=[]
for sid,h in a[a.Source.isin(['KFM','LTER'])].groupby('series_id'):
 r=h.iloc[0];has_baseline=bool(h.baseline_design_verified.all()) and (h.t<0).sum()>=3;h=h[(h.t!=0)&((h.t>0)|has_baseline)].copy()
 if h.empty:continue
 r=h.iloc[0];z=units[(units.MPA==r.MPA)&(units.Source==r.Source)&(units.Taxa==r.Taxa)&(units.Resp==r.Resp)&units.year.isin(h.year)].copy()
 if z.empty:continue
 def ag(q):return q.groupby(['site','status','year'],as_index=False).value.mean().groupby(['status','year'],as_index=False).value.mean()
 v=ag(z).pivot(index='year',columns='status',values='value').reindex(index=h.year,columns=['mpa','reference']);R=np.exp(h.lnDiff.to_numpy());offsets=np.divide(v.mpa.to_numpy()-R*v.reference.to_numpy(),R-1,out=np.full(len(R),np.nan),where=abs(R-1)>1e-6);good=np.isfinite(offsets)&(offsets>0);offset=float(np.median(offsets[good])) if good.any() else np.nan
 # Ratio equivalence across every year is mandatory; raw agreement is also checked where archived.
 err=np.max(abs(np.log((v.mpa.to_numpy()+offset)/(v.reference.to_numpy()+offset))-h.lnDiff.to_numpy())) if v.notna().all().all() and np.isfinite(offset) else np.inf
 old=raw[(raw.MPA==r.MPA)&(raw.Source==r.Source)&(raw.Taxa==r.Taxa)&(raw.Resp==r.Resp)].drop_duplicates(['year','status'],keep='last');j=ag(z).merge(old,on=['year','status'],suffixes=('_new','_old'));re=float((j.value_new-j.value_old).abs().max()) if len(j) else np.nan
 verified=err<1e-7 and (not np.isfinite(re) or re<1e-7)
 all_years=set(h.year);q=z.groupby(['site','unit']).filter(lambda g:set(g.year)==all_years);v2=ag(q).pivot(index='year',columns='status',values='value').reindex(index=h.year,columns=['mpa','reference']);complete=v2.notna().all().all();status='Eligible' if verified and complete else ('Reconstruction mismatch' if not verified else 'No complete fixed footprint')
 checks.append(dict(series_id=sid,status=status,max_ratio_error=err,max_raw_error=re,n_raw_checks=len(j),offset=offset,units_before=z[['site','unit']].drop_duplicates().shape[0],units_after=q[['site','unit']].drop_duplicates().shape[0]))
 if status!='Eligible':continue
 hh=h.copy();hh['lnDiff']=np.log((v2.mpa.to_numpy()+offset)/(v2.reference.to_numpy()+offset));hh['scenario']='Fixed_all_program_units';hh['corrected']=True;patches.append(hh);q['series_id']=sid;foot.append(q)
pd.DataFrame(checks).to_csv(O/'other_program_reconstruction.csv',index=False);pd.concat(patches).to_csv(O/'other_program_patches.csv',index=False);pd.concat(foot).to_csv(O/'other_program_fixed_units.csv',index=False)
allinput=pd.read_csv(O/'annual_input.csv');allinput=allinput[allinput.scenario!='Fixed_all_program_units'];new=allinput[allinput.scenario=='Fixed_station_zones'].copy();new['scenario']='Fixed_all_program_units';pp=pd.concat(patches).set_index(['series_id','year']);ix=pd.MultiIndex.from_frame(new[['series_id','year']]);mask=ix.isin(pp.index);new.loc[mask,'lnDiff']=pp.loc[ix[mask],'lnDiff'].to_numpy();new.loc[mask,'corrected']=True;pd.concat([allinput,new]).to_csv(O/'annual_input.csv',index=False)
(O/'other_refit_hashes.json').write_text(json.dumps({p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in used},indent=2));print(pd.DataFrame(checks).groupby('status').size())
