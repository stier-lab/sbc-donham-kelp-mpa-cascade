from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
O=Path('outputs/pisco_pass_crane_20261101');D=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data');A=Path('outputs/pisco_fish_effort_20261101')
x=pd.read_csv(D/'PISCO/MLPA_kelpforest_fish_2024.csv',encoding='latin1',low_memory=False);x['date']=pd.to_datetime(x[['year','month','day']]);keys=['site','year','date','method','zone','transect'];x['level_clean']=x.level.astype(str).str.strip().str.upper()
u=pd.read_csv(A/'survey_units.csv');u.date=pd.to_datetime(u.date);single=u[~u.levels.str.contains(';')].copy();allgroups={k:g for k,g in x.groupby(keys,dropna=False)}
rows=[]
for r in single.itertuples():
 k=tuple(getattr(r,c) for c in keys);g=allgroups[k];other=x[x.site.eq(r.site)&x.year.eq(r.year)&x.zone.eq(r.zone)&x.transect.eq(r.transect)]
 desired=['MID','CNMD'] if r.levels=='BOT' else ['BOT'];alt=other[other.level_clean.isin(desired)]
 exact=alt[alt.date.eq(r.date)];near=alt[~alt.date.eq(r.date)];gaps=abs((near.date-r.date).dt.days)
 rows.append(dict(**dict(zip(keys,k)),level_present=r.levels,rows=len(g),observed_species=g.classcode.nunique(),no_org_rows=int(g.classcode.eq('NO_ORG').sum()),sheephead_rows=int(g.classcode.eq('SPUL').sum()),depth_min=g.depth.min(),depth_max=g.depth.max(),visibility_min=g.vis.min(),campaigns=';'.join(map(str,sorted(g.survey_year.unique()))),notes=' | '.join(sorted(g.notes.dropna().astype(str).unique())),observers=';'.join(sorted(g.observer.dropna().astype(str).unique())),opposite_same_date_rows=len(exact),opposite_other_dates=';'.join(sorted(near.date.dt.strftime('%Y-%m-%d').unique())),nearest_opposite_gap_days=gaps.min() if len(gaps) else np.nan))
s=pd.DataFrame(rows);s.to_csv(O/'single_pass_ledger.csv',index=False)
# Exact legacy paired support; do not count every mapped site as current inclusion.
flags=pd.read_csv(A/'affected_site_years.csv');links=flags[['series_id','MPA','Resp','site','year','status']].drop_duplicates();sc=s.merge(links,on=['site','year']);sc.to_csv(O/'single_pass_comparison_links.csv',index=False)
c=x[x.method.str.contains('CRANE',na=False)].copy();c.to_csv(O/'crane_source_records.csv',index=False)
cu=c.groupby(keys,dropna=False).agg(levels=('level_clean',lambda v:';'.join(sorted(v.unique()))),rows=('classcode','size'),sheephead_rows=('classcode',lambda v:int(v.eq('SPUL').sum())),no_org_rows=('classcode',lambda v:int(v.eq('NO_ORG').sum())),depth_min=('depth','min'),depth_max=('depth','max'),notes=('notes',lambda v:' | '.join(sorted(v.dropna().astype(str).unique())))).reset_index();cu.to_csv(O/'crane_units.csv',index=False);cu.merge(links,on=['site','year']).to_csv(O/'crane_comparison_links.csv',index=False)
print('single-pass levels',s.level_present.value_counts().to_dict(),'rows',len(s),'same-date other method',s.opposite_same_date_rows.gt(0).sum(),'other dates',s.nearest_opposite_gap_days.notna().sum());print('Detailed notes retained in single_pass_ledger.csv')
print('CRANE',len(c),'units',len(cu),'sites',c.site.nunique(),'SBI VRG2004',len(x[x.site.str.startswith('SBI')&x.campus.eq('VRG')&x.year.eq(2004)]));print(cu.groupby(['site','levels']).size().to_string())
(O/'input_hashes.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [D/'PISCO/MLPA_kelpforest_fish_2024.csv',A/'survey_units.csv',A/'affected_site_years.csv',O/'upstream_merge_fish.R']},indent=2)+'\n')
