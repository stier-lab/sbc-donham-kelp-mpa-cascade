from pathlib import Path
import pandas as pd,numpy as np,re,json,hashlib
O=Path('outputs/pisco_fish_effort_20261101');D=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026');I=Path('outputs/dataset_inventory_20261101')
a=pd.read_csv('outputs/researcher_scope_20261101/annual_input.csv');ps=a[(a.Source=='PISCO')&(a.Taxa=='S. pulcher')].copy()
f=pd.read_csv(D/'data/PISCO/MLPA_kelpforest_fish_2024.csv',encoding='latin1',low_memory=False);f=f[f.campus.isin(['UCSB','VRG']) & f.level.ne('CAN')].copy()
master=pd.read_csv(D/'data/PISCO/master_site_table_Emilyedit.csv',encoding='utf-8-sig');master=master[master.BaselineRegion.eq('SOUTH')].rename(columns={'CA_MPA_Name_Short':'MPA','site_status':'status'})[['site','MPA','status']].drop_duplicates(['site','MPA'])
text=(D/'code/R/01_utils.R').read_text();excluded=re.findall(r'"([^"]+)"',re.search(r'EXCLUDED_REFERENCE_SITES <- c\((.*?)\)',text,re.S).group(1))
sw=pd.read_csv(D/'data/PISCO/MLPA_kelpforest_swath_2024.csv',encoding='latin1',low_memory=False);sw=sw[(sw.campus=='VRG')|((sw.campus=='UCSB')&~sw.zone.isin(['DEEP','SPECIALO']))]
v=pd.concat([sw[sw.month.between(5,10)][['site','year']],f[['site','year']]]).drop_duplicates().groupby('site').size();master=master[master.site.isin(v[v>=5].index)&~master.site.isin(excluded)&master.MPA.isin(ps.MPA)&master.status.isin(['mpa','reference'])]
fixed=pd.read_csv('tables/table_s_pisco_fixed_coverage.csv');fixed=fixed[(fixed.scenario=='strict_min3')&fixed.eligible]
f=f[f.site.isin(master.site)].copy();f['date']=pd.to_datetime(f[['year','month','day']],errors='coerce');assert f.date.notna().all()
f['standard']=f.method.isin(['SBTL_FISH_PISCO','SBTL_FISH_VRG']) & f.LOOKEDFOR.eq('yes')
ev=['site','year','date','method','zone','transect'];unit=['site','year','zone','transect']
# Retain all-taxon survey evidence, even when sheephead was absent.
e=f.groupby(ev,dropna=False).agg(levels=('level',lambda x:';'.join(sorted(x.unique()))),standard=('standard','all')).reset_index();e['bot_upper']=e.levels.str.split(';').map(lambda x:'BOT' in x and bool({'MID','CNMD'}&set(x)))
e['bot_mid']=e.levels.str.split(';').map(lambda x: {'BOT','MID'}<=set(x))
fish=f[f.classcode.eq('SPUL')].copy();fish['Den']=fish['count'];fish['Bio']=fish['count']*.0144*fish.fish_tl**3.04
bad=fish[(fish['count']>0)&(fish.fish_tl.isna()|fish.fish_tl.le(0))];bad[['site','year','date','zone','transect','count','fish_tl']].to_csv(O/'missing_lengths.csv',index=False)
def summed(x,keys):
 g=x.groupby(keys,dropna=False)[['Den','Bio']].sum(min_count=1);badg=x.groupby(keys,dropna=False).Bio.apply(lambda s:s.isna().any());g.loc[badg,'Bio']=np.nan;return g.reset_index()
# Exact legacy missing-value behavior retained only in the reconstruction comparison.
l=f[unit].drop_duplicates().merge(summed(fish,unit),on=unit,how='left').fillna({'Den':0,'Bio':0});l['scenario']='Legacy'
legacy=l.groupby(['site','year','zone'],as_index=False)[['Den','Bio']].mean().groupby(['site','year'],as_index=False)[['Den','Bio']].mean();legacy['scenario']='Legacy'
tables=[legacy];events=[]
for scenario,sel,levels in [('Visit_mean',e,None),('Standard_visits',e[e.standard],None),('Complete_noncan',e[e.standard&e.bot_upper],None),('Complete_BOT',e[e.standard&e.bot_upper],['BOT']),('Matched_noncan',e[e.standard&e.bot_mid],None),('Matched_BOT_MID',e[e.standard&e.bot_mid],['BOT','MID']),('Matched_BOT',e[e.standard&e.bot_mid],['BOT'])]:
 ff=fish if levels is None else fish[fish.level.isin(levels)]
 vals=summed(ff,ev);vals['observed_species']=True
 g=sel.merge(vals,on=ev,how='left');absent=g.observed_species.isna();g.loc[absent,['Den','Bio']]=0
 # Keep site-year mass unknown if a positive-count missing length belongs to its footprint.
 invalid=g[g.Bio.isna()][['site','year']].drop_duplicates()
 g['scenario']=scenario;events.append(g)
 # Preserve original hierarchy: average repeat dates within unit, then equal units/zones.
 # Do not give days with one zone the same weight as full-site survey days.
 z=g.groupby(unit,as_index=False)[['Den','Bio']].mean().groupby(['site','year','zone'],as_index=False)[['Den','Bio']].mean().groupby(['site','year'],as_index=False)[['Den','Bio']].mean()
 for b in invalid.itertuples():z.loc[z.site.eq(b.site)&z.year.eq(b.year),'Bio']=np.nan
 z['scenario']=scenario;tables.append(z)
sv=pd.concat(tables);sv[['Den','Bio']]/=60;sv.to_csv(O/'site_year_values.csv',index=False)
evdata=pd.concat(events);evdata.to_csv(O/'event_values.csv',index=False)
e.to_csv(O/'survey_units.csv',index=False)
# Source flags joined to actual mapped members and annual analysis windows.
flags=[];validation=[];role=[];paired=[];coverage=[]
oldcheck=pd.read_csv('outputs/seasonal_zone_audit_20261101/reconstruction_checks.csv').set_index('series_id')
for sid,h in ps.groupby('series_id'):
 h=h.sort_values('year');r=h.iloc[0];mp=master[master.MPA.eq(r.MPA)].copy()
 if r.MPA in set(fixed.MPA):
  ff=fixed[fixed.MPA.eq(r.MPA)].iloc[0];mp=mp[mp.site.isin((ff.protected+';'+ff.reference).split(';'))]
 vals=sv.merge(mp,on='site');vals=vals[vals.year.isin(h.year)].rename(columns={r.Resp:'value'})
 q=vals.groupby(['scenario','year','status'],as_index=False).value.mean();q['series_id']=sid;role.append(q)
 pp=q.pivot(index=['scenario','year'],columns='status',values='value');offset=oldcheck.loc[sid,'offset']
 orig=pp.loc['Legacy'].reindex(h.year);ln=np.log((orig.mpa+offset)/(orig.reference+offset));err=np.max(abs(ln.to_numpy()-h.lnDiff.to_numpy()))
 validation.append(dict(series_id=sid,offset=offset,max_ratio_error=err,verified=bool(np.isfinite(err) and err<1e-7)))
 for sc in sv.scenario.unique():
  ss=vals[vals.scenario.eq(sc)];rows=[]
  # Same role/site/year membership as Legacy, never silently drop a station for a missing response.
  baseline=vals[vals.scenario.eq('Legacy')][['site','year','status']]
  z=baseline.merge(ss[['site','year','status','value']],on=['site','year','status'],how='left',validate='one_to_one')
  ok=z.groupby('year').value.apply(lambda x:x.notna().all());good=set(ok[ok].index)
  p=pp.loc[sc].reindex(h.year)
  hh=h.copy();hh['scenario']=sc;hh['value_mpa']=p.mpa.to_numpy();hh['value_reference']=p.reference.to_numpy();hh['lnDiff']=np.log((hh.value_mpa+offset)/(hh.value_reference+offset));hh['complete_support']=hh.year.isin(good);hh['offset']=offset
  paired.append(hh)
  coverage.append(dict(series_id=sid,scenario=sc,original_years=len(h),complete_years=int(hh.complete_support.sum()),lost_years=';'.join(map(str,hh.loc[~hh.complete_support,'year'])),max_ratio_change=float((hh.lnDiff-h.lnDiff.to_numpy()).abs().max())))
 ee=e.merge(mp,on='site');ee=ee[ee.year.isin(h.year)]
 for (site,year,status),g in ee.groupby(['site','year','status']):
  nrep=int((g.groupby(['zone','transect']).date.nunique()>1).sum());comb=' | '.join(sorted(g.levels.unique()));nl=g.levels.str.split(';').map(len)
  flags.append(dict(series_id=sid,MPA=r.MPA,Resp=r.Resp,site=site,year=year,status=status,dates=g.date.nunique(),revisited_units=nrep,level_combinations=comb,unequal_levels=nl.nunique()>1,units=len(g),units_without_bot_mid=int((~g.bot_mid).sum()),nonstandard_units=int((~g.standard).sum())))
pd.DataFrame(validation).to_csv(O/'reconstruction_checks.csv',index=False);pd.DataFrame(coverage).to_csv(O/'coverage_changes.csv',index=False);pd.DataFrame(flags).to_csv(O/'affected_site_years.csv',index=False);pd.concat(role).to_csv(O/'role_values.csv',index=False)
d=pd.concat(paired,ignore_index=True);d.to_csv(O/'paired_values.csv',index=False)
print(pd.DataFrame(validation).to_string(index=False));print(pd.DataFrame(coverage).groupby('scenario')[['original_years','complete_years']].sum())

# Contrasts isolate each decision; preserve exact paired years for both arms.
contrasts=[('Repeat_visits','Legacy','Visit_mean'),('Survey_method','Visit_mean','Standard_visits'),('Complete_levels','Standard_visits','Complete_noncan'),('Strict_MID_footprint','Complete_noncan','Matched_noncan'),('Bottom_only','Complete_noncan','Complete_BOT'),('Combined','Legacy','Complete_noncan')]
model=[];pairchecks=[]
other=a[(a.Taxa=='S. pulcher') & (a.Source!='PISCO')]
for contrast,left,right in contrasts:
 for sid,h in d.groupby('series_id'):
  l=h[h.scenario.eq(left) & h.complete_support].copy();r=h[h.scenario.eq(right) & h.complete_support].copy();ys=set(l.year)&set(r.year);l=l[l.year.isin(ys)];r=r[r.year.isin(ys)]
  for arm,x in [('Before',l),('After',r)]:
   x=x.copy();x['contrast']=contrast;x['arm']=arm;x['scenario']=contrast+'__'+arm;model.append(x)
  pairchecks.append(dict(contrast=contrast,series_id=sid,n_years=len(ys),same_years=True))
 for arm in ['Before','After']:
  x=other.copy();x['contrast']=contrast;x['arm']=arm;x['scenario']=contrast+'__'+arm;model.append(x)
pd.concat(model,ignore_index=True).to_csv(O/'model_input.csv',index=False);pd.DataFrame(pairchecks).to_csv(O/'matched_year_checks.csv',index=False)
# Record calendar/provider-campaign context rather than interpreting every date as a repeat visit.
cc=f.groupby(['site','year','survey_year','month','day','method','zone','level'],dropna=False).transect.nunique().reset_index(name='n_transects');cc.to_csv(O/'campaign_date_ledger.csv',index=False)
inputs=[D/'data/PISCO/MLPA_kelpforest_fish_2024.csv',D/'data/PISCO/MLPA_kelpforest_swath_2024.csv',D/'data/PISCO/master_site_table_Emilyedit.csv',D/'code/R/04_pisco_processing.R',D/'code/R/01_utils.R',Path('outputs/researcher_scope_20261101/annual_input.csv'),Path('tables/table_s_pisco_fixed_coverage.csv'),Path('outputs/seasonal_zone_audit_20261101/reconstruction_checks.csv')]
(O/'input_hashes.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},indent=2)+'\n')
