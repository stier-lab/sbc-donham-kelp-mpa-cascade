from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,itertools,shutil,os
root=Path.cwd();out=root/'outputs/upstream_sites_20261101';base=root/'outputs/diver_exclusion_refit_20261101/consistent_kfm';dr=Path(os.environ.get('MPA_DATA_REPO',str(root.parent/'Donham-Stier-CA-MPA-Data-2026')))
targets=['Anacapa Island SMCA','Painted Cave SMCA','Point Dume SMR','Point Dume SMCA'];old=pd.read_csv(base/'annual_input.csv')
tax={'MACPYRAD':'M. pyrifera','MESFRA':'M. franciscanus','MESFRAAD':'M. franciscanus','STRPUR':'S. purpuratus','STRPURAD':'S. purpuratus','PANINT':'P. interruptus','SPUL':'S. pulcher','Macrocystis pyrifera':'M. pyrifera','Strongylocentrotus franciscanus':'M. franciscanus','Mesocentrotus franciscanus':'M. franciscanus','Strongylocentrotus purpuratus':'S. purpuratus','Panulirus interruptus':'P. interruptus'}
ps=pd.read_csv(out/'pisco_site_values.csv');ps=ps.rename(columns={'CA_MPA_Name_Short':'MPA','site_status':'status'});ps['Taxa']=ps.y.map(tax);assert ps.Taxa.notna().all();ps['Source']='PISCO';ps['Bio']=np.where(ps.Taxa=='M. pyrifera',ps.biomass,ps.biomass/60);ps['Den']=ps.density
# Refresh the formerly excluded PISCO urchin branches with the restored size pool.
rebuild=pd.read_csv(out/'pisco_refreshed_urchin_biomass.csv').rename(columns={'CA_MPA_Name_Short':'MPA','biomass':'refreshed_biomass'})
ps=ps.merge(rebuild[['MPA','site','year','y','refreshed_biomass']],on=['MPA','site','year','y'],how='left',validate='many_to_one')
mask=ps.MPA.isin(targets)&ps.Taxa.isin(['M. franciscanus','S. purpuratus'])
ps.loc[mask,'Bio']=ps.loc[mask,'refreshed_biomass']/60
bad=pd.read_csv(out/'invalid_sheephead_lengths.csv');badkeys=set(zip(bad.site,bad.year));invalid=(ps.Taxa=='S. pulcher')&pd.Series(list(zip(ps.site,ps.year)),index=ps.index).isin(badkeys);ps.loc[invalid,'Bio']=np.nan
# Source lobster tables carry duplicate branch rows with identical density and
# at most one finite biomass estimate; reconcile by identity, never row order.
ps=ps[ps.status.isin(['mpa','reference'])].copy();keys=['MPA','site','status','Source','year','Taxa']
dups=ps[ps.duplicated(keys,keep=False)];dups.to_csv(out/'pisco_duplicate_source_rows.csv',index=False)
assert (ps.groupby(keys)[['Bio','Den']].nunique()<=1).all().all()
ps=ps.groupby(keys,as_index=False)[['Bio','Den']].mean()
p=ps.melt(id_vars=['MPA','site','status','Source','year','Taxa'],value_vars=['Bio','Den'],var_name='Resp',value_name='value');p=p[~((p.Taxa=='M. pyrifera')&(p.Resp=='Den'))];assert not p.duplicated(['MPA','site','status','Source','year','Taxa','Resp']).any();p['branch']=np.where(p.Taxa=='S. pulcher','fish','benthic')
# Verify independent unit conversion/aggregation against current retained control sites.
raw=pd.read_csv(root/'data/harmonized/harmonized_raw_responses.csv').rename(columns={'CA_MPA_Name_Short':'MPA','source':'Source','resp':'Resp'});fullnames={str(x).strip():y for x,y in old[['y','Taxa']].drop_duplicates().itertuples(index=False,name=None)};raw['Taxa']=raw.taxon_name.str.strip().map(fullnames)
agg=p.groupby(['MPA','Source','Taxa','Resp','year','status'],as_index=False).value.mean()
control=agg[agg.MPA.isin(['Abalone Cove SMCA','Matlahuayl SMR'])].merge(raw,on=['MPA','Source','Taxa','Resp','year','status'],suffixes=('_rebuilt','_current'))
control['difference']=control.value_rebuilt-control.value_current;assert len(control)>100;assert max(abs(control.difference))<1e-7;control.to_csv(out/'control_reconstruction_check.csv',index=False)
p=p[p.MPA.isin(targets)].copy();p.to_csv(out/'pisco_candidate_station_values.csv',index=False)
# Fixed branch coverage is based only on valid observations, not their magnitudes.
meta=pd.read_csv(root/'data/harmonized/harmonized_site_metadata.csv').rename(columns={'CA_MPA_Name_Short':'MPA'}).drop_duplicates('MPA');starts=meta.set_index('MPA').MPA_Start.to_dict();choices=[];feasible=[];annualparts=[];membership=[]
def subsets(v):
 v=sorted(v)
 return itertools.chain.from_iterable(itertools.combinations(v,i) for i in range(1,len(v)+1))
for (m,branch),z in p.groupby(['MPA','branch']):
 required=2 if branch=='fish' else 7
 counts=z[np.isfinite(z.value)].groupby(['site','year']).size();complete=counts[counts==required].reset_index()[['site','year']]
 statuses=z[['site','status']].drop_duplicates();assert not statuses.site.duplicated().any();pa=statuses.loc[statuses.status=='mpa','site'];re=statuses.loc[statuses.status=='reference','site'];options=[]
 for pp in subsets(pa):
  for rr in subsets(re):
   ss=pp+rr;yrs=sorted(set.intersection(*(set(complete.loc[complete.site==s,'year']) for s in ss)))
   nb=sum(y<starts[m] for y in yrs);na=sum(starts[m]<y<=starts[m]+9 for y in yrs)
   if nb>=3 and na>=3:options.append(dict(MPA=m,branch=branch,protected=';'.join(pp),reference=';'.join(rr),n_stations=len(ss),n_years=len(yrs),nb=nb,na=na,years=';'.join(map(str,yrs))))
 if options:
  opt=sorted(options,key=lambda a:(-a['n_stations'],-a['n_years'],a['protected'],a['reference']))[0];ss=(opt['protected']+';'+opt['reference']).split(';');yy=list(map(int,opt['years'].split(';')));zz=z[z.site.isin(ss)&z.year.isin(yy)].copy();verified=True;feasible.extend(options);choices.append(dict(**opt,baseline_design_verified=True))
 else:zz=z.copy();verified=False;choices.append(dict(MPA=m,branch=branch,baseline_design_verified=False,reason='No fixed station group with all branch responses and >=3 strict Before and >=3 years1-9 After'))
 membership.append(zz);aa=zz.groupby(['MPA','Source','Taxa','Resp','year','status'],as_index=False).value.mean();aa['baseline_design_verified']=verified;annualparts.append(aa)
pd.DataFrame(choices).to_csv(out/'pisco_baseline_design.csv',index=False);pd.DataFrame(feasible).to_csv(out/'pisco_feasible_baselines.csv',index=False)
pd.concat(membership).to_csv(out/'pisco_selected_station_values.csv',index=False)
# KFM annual values: same source averaging, preserving missing biomass rather than treating it as zero.
k=pd.read_csv(out/'kfm_site_values.csv').rename(columns={'CA_MPA_Name_Short':'MPA','site_id':'site'});k=k[k.MPA.isin(targets)].copy();k['Taxa']=k.y.map(tax);assert k.Taxa.notna().all();k=k[(k.Taxa!='M. pyrifera')|(k.area==10)].copy();k['Source']='KFM';k['Bio']=k.biomass;k['Den']=np.where(k.Taxa=='M. pyrifera',k['den.y'],k.den)
kl=k.melt(id_vars=['MPA','site','status','Source','year','Taxa','area'],value_vars=['Bio','Den'],var_name='Resp',value_name='value');kl=kl[~((kl.Taxa=='M. pyrifera')&(kl.Resp=='Den'))];kl=kl[~((kl.Taxa=='P. interruptus')&(kl.Resp=='Bio'))];assert not kl.duplicated(['MPA','site','year','Taxa','Resp']).any()
ka=kl.groupby(['MPA','Source','Taxa','Resp','year','status'],as_index=False).value.agg(lambda x: x.mean(skipna=False));ka['baseline_design_verified']=True;annualparts.append(ka)
f=pd.read_csv(out/'kfm_fish_site_values.csv').rename(columns={'CA_MPA_Name_Short':'MPA','site_id':'site'});f=f[f.MPA.isin(targets)&(f.sample_method=='visualfish')].copy();f['Source']='KFM';f['Taxa']='S. pulcher';f['Resp']='Den';f['value']=f.den
fa=f.groupby(['MPA','Source','Taxa','Resp','year','status'],as_index=False).value.agg(lambda x:x.mean(skipna=False));fa['baseline_design_verified']=True;annualparts.append(fa)
kstation=pd.concat([kl,f[['MPA','site','status','Source','year','Taxa','Resp','value']]],ignore_index=True);kstation.to_csv(out/'kfm_selected_station_values.csv',index=False)
annual=pd.concat(annualparts,ignore_index=True);annual=annual[annual.status.isin(['mpa','reference'])];assert (annual.value.dropna()>=0).all();annual.to_csv(out/'unpaired_annual_values.csv',index=False)
parts=[];ledger=[]
for key,z in annual.groupby(['MPA','Source','Taxa','Resp']):
 m,source,taxa,resp=key;positive=z.value[z.value>0];floor=positive.min()/2 if len(positive) else np.nan
 assert not z.duplicated(['year','status']).any();v=z.pivot(index='year',columns='status',values='value').reindex(columns=['mpa','reference']).dropna().reset_index();pre=int((v.year<starts[m]).sum());post=int(v.year.between(starts[m]+1,starts[m]+9).sum());verified=bool(z.baseline_design_verified.all())
 ledger.append(dict(MPA=m,Source=source,Taxa=taxa,Resp=resp,n_paired=len(v),before=pre,post_years1_9=post,years=';'.join(map(str,v.year)),zero_correction=floor,baseline_design_verified=verified,coverage_eligible=post>=3,reason='No positive observations for a corrected ratio' if not np.isfinite(floor) else ('Fewer than3paired post years' if post<3 else 'Coverage passes; model variance checked separately')))
 if len(v)==0 or not np.isfinite(floor):continue
 v['lnDiff']=np.log((v.mpa+floor)/(v.reference+floor));mx=z.value.max();v['mpa']=(v.mpa+floor)/mx;v['reference']=(v.reference+floor)/mx;v['Diff']=np.exp(v.lnDiff);v['MPA']=m;v['Source']=source;v['Taxa']=taxa;v['Resp']=resp;v['MPA_Start']=starts[m];v['baseline_design_verified']=verified;parts.append(v)
new=pd.concat(parts,ignore_index=True);new=new.merge(meta[['MPA','type','Location']],on='MPA',validate='many_to_one');new['t']=new.year-new.MPA_Start;new['time']=new.t.clip(lower=0);new['BA']=np.where(new.t<=0,'Before','After');new['series_id']=new[['MPA','Taxa','Source','Resp']].agg('|'.join,axis=1);names=old[['Taxa','y']].drop_duplicates().set_index('Taxa').y.to_dict();new['y']=new.Taxa.map(names);new['Hectares']=np.nan
assert not set(new.series_id)&set(old.series_id);assert not new.duplicated(['series_id','year']).any();pd.DataFrame(ledger).to_csv(out/'candidate_eligibility.csv',index=False)
# Shared-reference connected groups: preserve all old edges and add new source footprints.
oldclusters=pd.read_csv(base/'reference_clusters.csv');parent={m:m for m in set(old.MPA)|set(new.MPA)}
def find(x):
 while parent[x]!=x:x=parent[x]
 return x
def union(a,b):
 a,b=find(a),find(b)
 if a!=b:parent[max(a,b)]=min(a,b)
for _,z in oldclusters.groupby('cluster'):
 for m in z.MPA:union(z.MPA.iloc[0],m)
refs=pd.read_csv(root/'data/harmonized/reference_sites.csv');refs=refs[refs.included&~refs.excluded];edges=refs[['MPA','reference_id']].drop_duplicates().copy()
selected=pd.concat([pd.concat(membership)[['MPA','Source','site','status','value']],kstation[['MPA','Source','site','status','value']]],ignore_index=True);selected=selected[(selected.status=='reference')&np.isfinite(selected.value)].drop_duplicates(['MPA','Source','site'])
selected['reference_id']='site|'+selected.site.astype(str)
edges=pd.concat([edges,selected[['MPA','reference_id']]],ignore_index=True)
# Match aliases already validated in the registry by reference name; preserve their physical IDs.
for _,r in selected.iterrows():
 aliases=refs[(refs.Source==r.Source)&((refs.site_id==r.site)|(refs.reference_name==r.site))].reference_id.unique()
 for alias in aliases:edges.loc[len(edges)]=[r.MPA,alias]
for _,z in edges.groupby('reference_id'):
 ms=[m for m in z.MPA.unique() if m in parent]
 for m in ms[1:]:union(ms[0],m)
cl=pd.DataFrame({'MPA':sorted(parent)});cl['cluster']=cl.MPA.map(find);cl.to_csv(out/'reference_clusters.csv',index=False);edges.to_csv(out/'reference_connections.csv',index=False);selected.to_csv(out/'selected_reference_sites.csv',index=False)
full=pd.concat([old.drop(columns='cluster'),new],ignore_index=True).merge(cl,on='MPA',validate='many_to_one');full=full.reindex(columns=old.columns);full.to_csv(out/'annual_input.csv',index=False);new.to_csv(out/'new_annual_histories.csv',index=False)
shutil.copy(base/'protected_hashes.json',out/'protected_hashes.json')
plan=(base/'PLAN.md').read_text()+'\nUpstream reconstruction extension: add any estimable focal responses from Anacapa SMCA, Painted Cave, Point Dume SMR and Point Dume SMCA using current source methods/calibrations. PISCO reference exclusions and source minimum-five-year station filter retained; KFM Keyhole restored for review because the inherited benthic exclusion has no recovered methodological justification. Modern KFM kelp protocol only. PISCO baseline admission additionally requires a fixed station group with every branch response observed in >=3strict Before and >=3primary After years; rank eligible groups by station count, then total common years, then alphabetical station names. Never use effect magnitudes. Correction constants use all valid unpaired annual values before matching. Missing fish lengths invalidate biomass at that station/year, not its density. Keep all prior annual outcomes unchanged. Shared reference footprints update robust clusters. Anacapa reference size memberships are explicitly restored from source station assignments rather than limited to existing PISCO-table MPA labels. Recompute candidate PISCO urchin biomass with that complete pool. Additional sensitivity omits Anacapa SMCA KFM to assess dependence on restoring Keyhole.\n'
(out/'PLAN.md').write_text(plan)
files=[out/n for n in ['PLAN.md','annual_input.csv','reference_clusters.csv','candidate_eligibility.csv','pisco_baseline_design.csv','unpaired_annual_values.csv']]
(out/'freeze_manifest.json').write_text(json.dumps({str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest() for f in files},indent=2))
print(pd.DataFrame(ledger).to_string(index=False));print('Candidate roster',full.series_id.nunique())
