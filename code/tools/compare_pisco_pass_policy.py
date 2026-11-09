from pathlib import Path
import pandas as pd,numpy as np
O=Path('outputs/pisco_pass_crane_20261101');A=Path('outputs/pisco_fish_effort_20261101')
e=pd.read_csv(A/'event_values.csv');e=e[e.scenario.eq('Visit_mean')].copy();e=e[e.levels.str.split(';').map(lambda l:'BOT' in l and bool(set(l)&{'MID','CNMD'}))]
unit=['site','year','zone','transect'];bad=e[e.Bio.isna()][['site','year']].drop_duplicates()
s=e.groupby(unit,as_index=False)[['Den','Bio']].mean().groupby(['site','year','zone'],as_index=False)[['Den','Bio']].mean().groupby(['site','year'],as_index=False)[['Den','Bio']].mean();s[['Den','Bio']]/=60
for r in bad.itertuples():s.loc[s.site.eq(r.site)&s.year.eq(r.year),'Bio']=np.nan
links=pd.read_csv(A/'affected_site_years.csv')[['series_id','MPA','Resp','site','year','status']].drop_duplicates();v=links.merge(s,on=['site','year'],how='left',validate='many_to_one');v['value']=np.where(v.Resp.eq('Bio'),v.Bio,v.Den)
ok=v.groupby(['series_id','year']).value.apply(lambda x:x.notna().all()).rename('valid').reset_index();role=v.groupby(['series_id','year','status'],as_index=False).value.mean();p=role.pivot(index=['series_id','year'],columns='status',values='value').reset_index();p=p.merge(ok,on=['series_id','year'])
base=pd.read_csv(A/'paired_values.csv');z=base[base.scenario.eq('Visit_mean')].merge(p,on=['series_id','year'],validate='one_to_one');z['complete_support']=z.valid;z['lnDiff']=np.log((z.mpa_y+z.offset)/(z.reference_y+z.offset)) if 'mpa_y' in z else np.log((z.mpa+z.offset)/(z.reference+z.offset));z['scenario']='Keep_CRANE_complete';z.to_csv(O/'proposed_paired_values.csv',index=False)
# Add one isolated pass-completeness contrast on identical years, and one explicit coverage contrast for CRANE.
annual=pd.read_csv('outputs/researcher_scope_20261101/annual_input.csv');other=annual[annual.Taxa.eq('S. pulcher')&annual.Source.ne('PISCO')];parts=[]
for contrast,left,right,matched in [('Incomplete_passes',base[base.scenario.eq('Visit_mean')],z,True),('CRANE_coverage',base[base.scenario.eq('Complete_noncan')],z,False)]:
 for sid,r in right.groupby('series_id'):
  l=left[left.series_id.eq(sid)&left.complete_support].copy();r=r[r.complete_support].copy()
  if matched:
   ys=set(l.year)&set(r.year);l=l[l.year.isin(ys)];r=r[r.year.isin(ys)]
  for arm,g in [('Before',l),('After',r)]:
   g=g.copy();g['contrast']=contrast;g['arm']=arm;g['scenario']=contrast+'__'+arm;parts.append(g)
 for arm in ['Before','After']:
  g=other.copy();g['contrast']=contrast;g['arm']=arm;g['scenario']=contrast+'__'+arm;parts.append(g)
pd.concat(parts,ignore_index=True).to_csv(O/'model_input.csv',index=False)
print('proposed complete history-years',int(z.complete_support.sum()),'of',len(z));print(z[~z.complete_support][['MPA','Resp','year']].to_string(index=False))
# Do not recode or remove source rows. This is a comparison-level recommendation ledger.
single=pd.read_csv(O/'single_pass_comparison_links.csv');single['classification']=np.where(single.notes.fillna('').str.contains('NO BOTTOMS DUE TO POOR VIS'),'Confirmed omitted bottom pass','Opposite pass unrecorded; reason unresolved');single['proposed_use']='Retain measured pass; omit unit from combined bottom+upper response until supported. No synthetic zeros.';single.to_csv(O/'single_pass_comparison_links.csv',index=False)
