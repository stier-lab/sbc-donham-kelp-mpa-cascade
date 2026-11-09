from pathlib import Path
import pandas as pd,numpy as np,json
O=Path('outputs/all_program_station_refits_20261102');a=pd.read_csv('outputs/pisco_fish_rules_20261101/approved_annual_input.csv');tax={'MACPYRAD':'M. pyrifera','MESFRA':'M. franciscanus','MESFRAAD':'M. franciscanus','STRPUR':'S. purpuratus','STRPURAD':'S. purpuratus','PANINT':'P. interruptus','SPUL':'S. pulcher','Macrocystis pyrifera':'M. pyrifera','Strongylocentrotus franciscanus':'M. franciscanus','Mesocentrotus franciscanus':'M. franciscanus','Strongylocentrotus purpuratus':'S. purpuratus','Panulirus interruptus':'P. interruptus'}
parts=[]
def long(x,source):
 x=x.copy();x['Source']=source
 if 'Taxa' not in x:x['Taxa']=x.y.map(tax)
 assert x.Taxa.notna().all()
 x=x[x.status.isin(['mpa','reference'])];x=x.melt(id_vars=['MPA','site','status','Source','year','Taxa'],value_vars=['Bio','Den'],var_name='Resp',value_name='value');return x[~(x.Taxa.eq('M. pyrifera')&x.Resp.eq('Den'))]
x=pd.read_csv(O/'pisco_sites.csv').rename(columns={'CA_MPA_Name_Short':'MPA','site_status':'status'});x['Bio']=np.where(x.y.eq('MACPYRAD'),x.biomass,x.biomass/60);x['Den']=x.density;x=long(x,'PISCO');x=x[x.Taxa.ne('S. pulcher')]
up=pd.read_csv('outputs/upstream_sites_20261101/pisco_selected_station_values.csv');up=up[up.MPA.isin(['Point Dume SMR','Point Dume SMCA'])&up.Taxa.ne('S. pulcher')];x=x[~x.MPA.isin(up.MPA)];parts.extend([x,up[x.columns]])
f=pd.read_csv('outputs/pisco_fish_rules_20261101/station_values.csv');f['Taxa']='S. pulcher';f['Source']='PISCO';parts.append(f[['MPA','site','status','Source','year','Taxa','Resp','value']])
x=pd.read_csv(O/'kfm_sites.csv').rename(columns={'CA_MPA_Name_Short':'MPA','site_id':'site'});x=x[(x.y!='Macrocystis pyrifera')|(x.area==10)].copy();x['Bio']=x.biomass;x['Den']=np.where(x.y.eq('Macrocystis pyrifera'),x['den.y'],x.den);parts.append(long(x,'KFM'))
x=pd.read_csv(O/'kfm_fish_sites.csv').rename(columns={'CA_MPA_Name_Short':'MPA','site_id':'site'});x=x[x.sample_method=='visualfish'].copy();x['Taxa']='S. pulcher';x['Bio']=np.nan;x['Den']=x.den;parts.append(long(x,'KFM'))
x=pd.read_csv(O/'Swath.LTER.csv').rename(columns={'CA_MPA_Name_Short':'MPA','site_status':'status'});x['Bio']=x.biomass;x['Den']=x.Density;parts.append(long(x,'LTER'))
for name,tx,area,den in [('lter.lob.site.ave','P. interruptus',300,'COUNT'),('lter.macro.site.ave','M. pyrifera',1,'den'),('lter.fish.ave','S. pulcher',80,'COUNT')]:
 x=pd.read_csv(O/(name+'.csv')).rename(columns={'CA_MPA_Name_Short':'MPA','SITE':'site','YEAR':'year'});x['Taxa']=tx;x['Bio']=x.biomass/area;x['Den']=x[den]/area;parts.append(long(x,'LTER'))
D=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data');x=pd.read_csv(D/'LANDSAT/MPA_Runs_new.csv');bad=pd.read_csv(D/'landsat_empty_reference_corrections.csv');x=x.merge(bad[['MPA','status','rep']].assign(exclude=True),on=['MPA','status','rep'],how='left');x=x[x.exclude.isna()];x['site']=x.MPA+'|'+x.status+'|'+x.rep.astype(str);x=x.melt(id_vars=['MPA','site','status'],value_vars=[c for c in x if c.isdigit()],var_name='year',value_name='value');x['year']=x.year.astype(int);x['Source']='Landsat';x['Taxa']='M. pyrifera';x['Resp']='Bio';parts.append(x)
v=pd.concat(parts,ignore_index=True);v['series_id']=v[['MPA','Taxa','Source','Resp']].agg('|'.join,axis=1);v=v.merge(a[['series_id','year','t','MPA_Start']],on=['series_id','year'],validate='many_to_one');keys=['series_id','site','status','year'];dups=v[v.duplicated(keys,keep=False)];dups.to_csv(O/'duplicate_station_records.csv',index=False);assert (v.groupby(keys).value.nunique()<=1).all();v=v.groupby(['series_id','MPA','Taxa','Source','Resp','site','status','year','t','MPA_Start'],as_index=False).value.mean()
v.to_csv(O/'candidate_station_values.csv',index=False)
checks=[];verified=[];offsets=[]
for sid,g in a.groupby('series_id'):
 z=v[v.series_id.eq(sid)];means=z.groupby(['year','status']).value.mean().unstack().reindex(columns=['mpa','reference']);j=g.merge(means,on='year',how='left');q=np.exp(j.lnDiff);c=(j.mpa-q*j.reference)/(q-1);cc=c[np.isfinite(c)&(abs(q-1)>1e-8)];off=float(cc.median()) if len(cc) else np.nan
 if off>=-1e-10:off=max(0.,off)
 err=np.max(abs(np.log((j.mpa+off)/(j.reference+off))-j.lnDiff)) if np.isfinite(off) and off>=0 else np.nan
 ok=bool(len(j)==len(g) and j[['mpa','reference']].notna().all().all() and np.isfinite(err) and err<1e-7)
 checks.append(dict(series_id=sid,Source=g.Source.iloc[0],MPA=g.MPA.iloc[0],Taxa=g.Taxa.iloc[0],Resp=g.Resp.iloc[0],n_years=len(g),n_reconstructed=int(j[['mpa','reference']].notna().all(axis=1).sum()),offset=off,max_ratio_error=err,verified=ok))
 if ok:verified.append(z);offsets.append(dict(series_id=sid,offset=off))
pd.DataFrame(checks).to_csv(O/'reconstruction_checks.csv',index=False);pd.concat(verified).to_csv(O/'verified_station_values.csv',index=False);pd.DataFrame(offsets).to_csv(O/'offsets.csv',index=False);print(pd.DataFrame(checks).groupby(['Source','verified']).size());print(pd.DataFrame(checks).query('not verified').to_string(index=False))
