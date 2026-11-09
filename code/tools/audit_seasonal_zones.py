from pathlib import Path
import pandas as pd,numpy as np,re,json,hashlib,subprocess
O=Path('outputs/seasonal_zone_audit_20261101');D=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026')
a=pd.read_csv('outputs/researcher_scope_20261101/annual_input.csv'); ps=a[a.Source=='PISCO']; used=[]
def read(p,**kw):
 used.append(str(p));return pd.read_csv(p,low_memory=False,**kw)
sw=read(D/'data/PISCO/MLPA_kelpforest_swath_2024.csv',encoding='latin1');sw=sw[(sw.campus=='VRG')|((sw.campus=='UCSB')&~sw.zone.isin(['DEEP','SPECIALO']))].copy()
fi=read(D/'data/PISCO/MLPA_kelpforest_fish_2024.csv',encoding='latin1');fi=fi[fi.campus.isin(['VRG','UCSB'])&(fi.level!='CAN')].copy()
master=read(D/'data/PISCO/master_site_table_Emilyedit.csv',encoding='utf-8-sig');master=master[master.BaselineRegion=='SOUTH'].rename(columns={'CA_MPA_Name_Short':'MPA','site_status':'status'})[['site','MPA','status']].drop_duplicates(['site','MPA'])
text=(D/'code/R/01_utils.R').read_text();excluded=re.findall(r'"([^"]+)"',re.search(r'EXCLUDED_REFERENCE_SITES <- c\((.*?)\)',text,re.S).group(1))
valid=pd.concat([sw[sw.month.between(5,10)][['site','year']],fi[['site','year']]]).drop_duplicates().groupby('site').size();master=master[master.site.isin(valid[valid>=5].index)&~master.site.isin(excluded)&master.MPA.isin(ps.MPA)&master.status.isin(['mpa','reference'])]
sw=sw[sw.site.isin(master.site)].copy();fi=fi[fi.site.isin(master.site)].copy()
fixed=read(Path('tables/table_s_pisco_fixed_coverage.csv'));fixed=fixed[(fixed.scenario=='strict_min3')&fixed.eligible]
raw=read(Path('data/harmonized/harmonized_raw_responses.csv')).rename(columns={'CA_MPA_Name_Short':'MPA','source':'Source','resp':'Resp'})
tax={'MESFRAAD':'M. franciscanus','STRPURAD':'S. purpuratus','PANINT':'P. interruptus','MACPYRAD':'M. pyrifera'}
name={'Mesocentrotus franciscanus':'M. franciscanus','Strongylocentrotus franciscanus':'M. franciscanus','Strongylocentrotus purpuratus':'S. purpuratus','Panulirus interruptus':'P. interruptus','Macrocystis pyrifera':'M. pyrifera','Semicossyphus pulcher':'S. pulcher'}
raw['Taxa']=raw.taxon_name.str.strip().map(name)
newraw=read(Path('outputs/upstream_sites_20261101/unpaired_annual_values.csv'));raw=pd.concat([raw[['MPA','Source','Taxa','Resp','year','status','value']],newraw[['MPA','Source','Taxa','Resp','year','status','value']]])
raw=raw.drop_duplicates(['MPA','Source','Taxa','Resp','year','status'],keep='last')
# Audit all mapped PISCO survey events before and after the inherited month filter.
cover=[]
for branch,x in [('benthic',sw),('fish',fi)]:
 c=x.groupby(['site','year','month','zone'],dropna=False).agg(transects=('transect','nunique'),depth_min=('depth','min'),depth_max=('depth','max')).reset_index();c=c.merge(master,on='site');c['branch']=branch;c['season_retained']=c.month.between(5,10) if branch=='benthic' else True;cover.append(c)
cov=pd.concat(cover);cov.to_csv(O/'pisco_survey_coverage.csv',index=False)
# Actual events define zeros; never manufacture a missing survey.
ev=['site','year','month','zone','transect'];z=sw[sw.month.between(5,10)].copy();event=z[ev].drop_duplicates();grid=event.merge(pd.DataFrame({'classcode':list(tax)}),how='cross')
f=z[z.classcode.isin(tax)].copy();assert (f['count'].dropna()>=0).all();badcounts=f[f['count'].isna()][['site','year','classcode']].drop_duplicates();badcounts.to_csv(O/'missing_source_counts.csv',index=False);f['amount']=f['count'];m=f.classcode=='MACPYRAD';f.loc[m,'amount']=f.loc[m,'count']*f.loc[m,'size'];badcounts=f[f.amount.isna()][['site','year','classcode']].drop_duplicates();badcounts.to_csv(O/'missing_source_counts_or_sizes.csv',index=False)
counts=f.groupby(ev+['classcode'],as_index=False).amount.sum();grid=grid.merge(counts,on=ev+['classcode'],how='left');grid.amount=grid.amount.fillna(0);grid['Taxa']=grid.classcode.map(tax);grid['Resp']=np.where(grid.classcode=='MACPYRAD','Bio','Den')
coef=read(Path('data/external/sbc1275/Macrocystis_Biomass_Relationships_20261016.csv'),encoding='utf-8-sig');coef=coef[coef.model==1].set_index('MONTH').slope
assert len(coef)==12;grid['value']=grid.amount/60;mk=grid.Resp=='Bio';grid.loc[mk,'value']*=grid.loc[mk,'month'].map(coef)*1000
zone=grid.groupby(['site','year','zone','transect','Taxa','Resp'],as_index=False).value.mean().groupby(['site','year','zone','Taxa','Resp'],as_index=False).value.mean();zone['branch']='benthic'
for v in badcounts.itertuples():zone.loc[(zone.site==v.site)&(zone.year==v.year)&(zone.Taxa==tax[v.classcode]),'value']=np.nan
# Fish branch reproduces the existing across-level sums. Repeated visits are audited separately.
f=fi[fi.classcode=='SPUL'].copy();f['Den']=f['count'];f['Bio']=f['count']*.0144*f.fish_tl**3.04
bad=f[(f['count']>0)&(~np.isfinite(f.fish_tl)|(f.fish_tl<=0))][['site','year']].drop_duplicates();bad['invalid_size']=True
keys=['site','year','zone','transect'];fg=fi[keys].drop_duplicates().merge(f.groupby(keys,as_index=False)[['Den','Bio']].sum(),on=keys,how='left').fillna({'Den':0,'Bio':0});fg=fg.merge(bad,on=['site','year'],how='left');fg.loc[fg.invalid_size==True,'Bio']=np.nan
fg=fg.melt(id_vars=keys,value_vars=['Den','Bio'],var_name='Resp',value_name='value');fg.value/=60;fz=fg.groupby(['site','year','zone','Resp'],as_index=False).value.mean();fz['Taxa']='S. pulcher';fz['branch']='fish';zone=pd.concat([zone,fz],ignore_index=True)
fi.groupby(keys).month.nunique().reset_index(name='n_months').query('n_months>1').to_csv(O/'fish_repeated_visit_keys.csv',index=False)
# Reuse frozen simulated masses only where cached zone totals reproduce archived
# raw annual biomass; this tests survey support, not new size distributions.
subprocess.run(['Rscript','-e','write.csv(readRDS("data/cache/pisco_urchin_bootstrap.rds"),"outputs/seasonal_zone_audit_20261101/cached_urchin_zones.csv",row.names=FALSE);write.csv(readRDS("data/cache/vrg_panint_bootstrap.rds"),"outputs/seasonal_zone_audit_20261101/cached_lobster_zones.csv",row.names=FALSE)'],check=True)
used.extend(['data/cache/pisco_urchin_bootstrap.rds','data/cache/vrg_panint_bootstrap.rds'])
cacheparts=[]
for file,tx in [('cached_urchin_zones.csv',None),('cached_lobster_zones.csv','P. interruptus')]:
 c=read(O/file);c['Taxa']=c.y.map(tax) if tx is None else tx;c['Resp']='Bio';c['branch']='benthic'
 c.loc[(c['count']==0)&c.biomass.isna(),'biomass']=0
 # Keep explicit unavailable positive-count mass missing at the whole site/year.
 bad=c[c.biomass.isna()][['site','year','Taxa']].drop_duplicates()
 c['cache_MPA']=c.CA_MPA_Name_Short if tx is None else ''
 c['value']=c.biomass/c.transect/60
 c=c[['site','year','zone','Taxa','Resp','branch','value','cache_MPA']].drop_duplicates()
 assert not c.duplicated(['site','year','zone','Taxa','cache_MPA']).any()
 for v in bad.itertuples():c.loc[(c.site==v.site)&(c.year==v.year)&(c.Taxa==v.Taxa),'value']=np.nan
 cacheparts.append(c)
zone=pd.concat([zone]+cacheparts,ignore_index=True)
zone.to_csv(O/'source_zone_values.csv',index=False)
validation=[];ledger=[];annuals=[];patches=[];footprints=[]
for sid,h in ps.groupby('series_id'):
 r=h.iloc[0];has_baseline=bool(h.baseline_design_verified.all()) and (h.t<0).sum()>=3;h=h[(h.t!=0)&((h.t>0)|has_baseline)].copy()
 if h.empty:continue
 r=h.iloc[0];branch='fish' if r.Taxa=='S. pulcher' else 'benthic';mp=master[master.MPA==r.MPA].copy()
 if branch=='fish' and r.MPA in set(fixed.MPA):
  ff=fixed[fixed.MPA==r.MPA].iloc[0];ss=(ff.protected+';'+ff.reference).split(';');mp=mp[mp.site.isin(ss)]
 zz=zone[(zone.Taxa==r.Taxa)&(zone.Resp==r.Resp)].merge(mp,on='site');zz=zz[zz.year.isin(h.year)&np.isfinite(zz.value)&(zz.cache_MPA.isna()|(zz.cache_MPA=='')|(zz.cache_MPA==r.MPA))]
 if zz.empty:
  ledger.append(dict(series_id=sid,status='Not reconstructed: size-based biomass requires zone-level size/count linkage'));continue
 def ag(q):return q.groupby(['site','status','year'],as_index=False).value.mean().groupby(['status','year'],as_index=False).value.mean()
 old=raw[(raw.MPA==r.MPA)&(raw.Source=='PISCO')&(raw.Taxa==r.Taxa)&(raw.Resp==r.Resp)]
 b=ag(zz);j=b.merge(old,on=['year','status'],suffixes=('_rebuilt','_old'));err=(j.value_rebuilt-j.value_old).abs().max() if len(j) else np.inf
 offset=old.loc[old.value>0,'value'].min()/2
 pivot=b.pivot(index='year',columns='status',values='value').reindex(h.year)
 # Recover the frozen additive constant from archived paired raw values and ratios;
 # this accommodates older pre-filter domains, and must reproduce every year.
 rr=np.exp(h.lnDiff.to_numpy());inferred=np.divide(pivot.mpa.to_numpy()-rr*pivot.reference.to_numpy(),rr-1,out=np.full(len(rr),np.nan),where=abs(rr-1)>1e-6);good=np.isfinite(inferred)&(abs(rr-1)>1e-6)
 if good.any() and np.nanmedian(inferred[good])>0:offset=float(np.nanmedian(inferred[good]))
 ln=np.log((pivot.mpa+offset)/(pivot.reference+offset));raterr=np.max(abs(ln.to_numpy()-h.lnDiff.to_numpy())) if ln.notna().all() else np.inf
 verified=len(j)==2*len(h) and err<1e-7 and raterr<1e-7
 validation.append(dict(series_id=sid,n=len(j),expected=2*len(h),max_raw_error=err,max_ratio_error=raterr,offset=offset,verified=verified))
 core={s:set.intersection(*(set(v.zone) for _,v in q.groupby('year'))) for s,q in zz.groupby('site')}
 fixed_sites=[s for s,q in zz.groupby('site') if set(q.year)==set(h.year) and len(core[s])]
 for scenario in ['Within_station_zones','Fixed_station_zones']:
  q=zz[[v.zone in core[v.site] for v in zz.itertuples()]].copy()
  if scenario=='Fixed_station_zones':q=q[q.site.isin(fixed_sites)]
  aa=ag(q);pp=aa.pivot(index='year',columns='status',values='value').reindex(index=h.year,columns=['mpa','reference']);complete=pp.notna().all().all()
  ledger.append(dict(series_id=sid,scenario=scenario,status='Eligible' if verified and complete else ('Reconstruction mismatch' if not verified else 'No complete fixed footprint'),n_sites_before=zz.site.nunique(),n_sites_after=q.site.nunique(),years=len(h)))
  fq=q[['site','status','year','zone']].drop_duplicates();fq['series_id']=sid;fq['scenario']=scenario;footprints.append(fq)
  if not (verified and complete):continue
  hh=h.copy();hh['lnDiff']=np.log((pp.mpa.to_numpy()+offset)/(pp.reference.to_numpy()+offset));hh['scenario']=scenario;patches.append(hh)
  aa['series_id']=sid;aa['scenario']=scenario;annuals.append(aa)
pd.DataFrame(validation).to_csv(O/'reconstruction_checks.csv',index=False);pd.DataFrame(ledger).to_csv(O/'eligibility.csv',index=False);pd.concat(footprints).to_csv(O/'fixed_footprints.csv',index=False);pd.concat(annuals).to_csv(O/'corrected_role_values.csv',index=False);pd.concat(patches).to_csv(O/'patches.csv',index=False)
# Full-cohort counterfactual: blocked histories are held unchanged, explicitly marked.
parts=[]
for scenario in ['Original','Within_station_zones','Fixed_station_zones']:
 q=a.copy();q['scenario']=scenario;q['corrected']=False
 if scenario!='Original':
  p=pd.concat(patches);p=p[p.scenario==scenario].set_index(['series_id','year']);idx=pd.MultiIndex.from_frame(q[['series_id','year']]);mask=idx.isin(p.index);q.loc[mask,'lnDiff']=p.loc[idx[mask],'lnDiff'].to_numpy();q.loc[mask,'corrected']=True
 parts.append(q)
pd.concat(parts).to_csv(O/'annual_input.csv',index=False)
(O/'input_hashes.json').write_text(json.dumps({p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in used},indent=2))
print(pd.DataFrame(validation)[['verified']].value_counts());print(pd.DataFrame(ledger).groupby(['scenario','status']).size().to_string())
