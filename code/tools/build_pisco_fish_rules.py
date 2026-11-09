"""Authorized PISCO measurement policy, isolated from manuscript/canonical outputs."""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
O=Path('outputs/pisco_fish_rules_20261101');D=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data');A=Path('outputs/pisco_fish_effort_20261101')
policy=dict(decision_date='2026-11-01',authorization='Adrian: keep CRANE, average repeat surveys, leave incomplete passes out of combined response',keep_crane=True,complete_levels='BOT plus MID or CNMD; no standalone CAN',repeat_aggregation='Mean repeat dates within site/calendar-year/zone/transect number; then equal transect, zone and site means',missing_lengths='Retain missing flag; do not zero-fill positive-count biomass',unchanged='Researcher exclusions, site/year membership, allometry, ratio offsets, calendar clock, baseline eligibility and uncertainty',paper_updated=False)
(O/'frozen_rules.json').write_text(json.dumps(policy,indent=2)+'\n')
x=pd.read_csv(D/'PISCO/MLPA_kelpforest_fish_2024.csv',encoding='latin1',low_memory=False);x=x[x.campus.isin(['UCSB','VRG'])&x.level.ne('CAN')].copy();x['date']=pd.to_datetime(x[['year','month','day']]);ev=['site','year','date','method','zone','transect'];unit=['site','year','zone','transect']
u=x.groupby(ev,dropna=False)['level'].apply(lambda s:';'.join(sorted(s.unique()))).reset_index(name='levels');u['complete']=u.levels.str.split(';').map(lambda s:'BOT' in s and bool({'MID','CNMD'}&set(s)))
f=x[x.classcode=='SPUL'].copy();f['Den']=f['count'];f['Bio']=f['count']*.0144*f.fish_tl**3.04;f.loc[(f['count']>0)&(f.fish_tl.isna()|f.fish_tl.le(0)),'Bio']=np.nan
v=f.groupby(ev,dropna=False)[['Den','Bio']].sum(min_count=1);bad=f.groupby(ev,dropna=False).Bio.apply(lambda s:s.isna().any());v.loc[bad,'Bio']=np.nan;v=v.reset_index();v['observed']=True
g=u[u.complete].merge(v,on=ev,how='left',validate='one_to_one');g.loc[g.observed.isna(),['Den','Bio']]=0
bad=g[g.Bio.isna()][['site','year']].drop_duplicates();s=g.groupby(unit,as_index=False)[['Den','Bio']].mean().groupby(['site','year','zone'],as_index=False)[['Den','Bio']].mean().groupby(['site','year'],as_index=False)[['Den','Bio']].mean();s[['Den','Bio']]/=60
for r in bad.itertuples():s.loc[s.site.eq(r.site)&s.year.eq(r.year),'Bio']=np.nan
links=pd.read_csv(A/'affected_site_years.csv')[['series_id','MPA','Resp','site','year','status']].drop_duplicates();v=links.merge(s,on=['site','year'],how='left',validate='many_to_one');v['value']=np.where(v.Resp.eq('Bio'),v.Bio,v.Den);v.to_csv(O/'station_values.csv',index=False)
valid=v.groupby(['series_id','year']).value.apply(lambda z:z.notna().all());r=v.groupby(['series_id','year','status']).value.mean().unstack();r['complete']=valid;r=r.reset_index();offset=pd.read_csv(A/'reconstruction_checks.csv').set_index('series_id').offset;r['offset']=r.series_id.map(offset);r['lnDiff']=np.log((r.mpa+r.offset)/(r.reference+r.offset));r.to_csv(O/'paired_values.csv',index=False)
# Verify previous proposal, independently reconstructed from raw fish here.
old=pd.read_csv('outputs/pisco_pass_crane_20261101/proposed_paired_values.csv');j=r.merge(old[['series_id','year','lnDiff','complete_support']],on=['series_id','year'],suffixes=('_new','_proposal'),validate='one_to_one');assert len(j)==420 and j.complete.eq(j.complete_support).all();err=(j.loc[j.complete,'lnDiff_new']-j.loc[j.complete,'lnDiff_proposal']).abs().max();assert err<1e-10
cols=['MPA','year','Taxa','Resp','Source','MPA_Start','series_id','t','baseline_design_verified','cluster','lnDiff']
a=pd.read_csv('outputs/researcher_scope_20261101/annual_input.csv')[cols];a['scenario']='Previous';b=a.copy();b['scenario']='Approved';patch=r.set_index(['series_id','year']);removed=[]
for idx,z in b.iterrows():
 k=(z.series_id,z.year)
 if k in patch.index:
  q=patch.loc[k]
  if not q.complete:removed.append(idx)
  else:b.loc[idx,'lnDiff']=q.lnDiff
b=b.drop(index=removed);a.loc[removed].to_csv(O/'removed_history_years.csv',index=False)
keys=set(map(tuple,b[['series_id','year']].to_numpy()));matched=a[[k in keys for k in zip(a.series_id,a.year)]].copy();parts=[]
for contrast,l,rr in [('Total',a,b),('Same_years',matched,b),('Coverage_only',a,matched)]:
 for arm,q in [('Before',l),('After',rr)]:
  q=q.copy();q['contrast']=contrast;q['arm']=arm;q['scenario']=contrast+'__'+arm;parts.append(q)
pd.concat(parts).to_csv(O/'model_input.csv',index=False);b.to_csv(O/'approved_annual_input.csv',index=False)
# Scope the event ledger to existing station/year memberships, preserving original source identity.
scope=links[['site','year']].drop_duplicates();ug=u.merge(scope,on=['site','year']);ug.to_csv(O/'survey_unit_decisions.csv',index=False)
inputs=[D/'PISCO/MLPA_kelpforest_fish_2024.csv',A/'affected_site_years.csv',A/'reconstruction_checks.csv',Path('outputs/researcher_scope_20261101/annual_input.csv'),Path('outputs/pisco_pass_crane_20261101/proposed_paired_values.csv')];(O/'input_hashes.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},indent=2)+'\n')
print('Prior proposal reproduced:',err,'; histories',a.series_id.nunique(),b.series_id.nunique(),'; annual rows',len(a),len(b));print(a.loc[removed,['MPA','Taxa','Resp','year']].to_string(index=False))
