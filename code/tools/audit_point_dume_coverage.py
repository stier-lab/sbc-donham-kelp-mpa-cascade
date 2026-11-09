from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
out=Path('outputs/point_dume_coverage_20261101');source=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data/PISCO/MLPA_kelpforest_swath_2024.csv')
x=pd.read_csv(source,encoding='latin1',low_memory=False);sites=['LECHUZA','NICHOLAS_CANYON_W'];x=x[x.site.isin(sites)].copy();x=x[(x.campus=='VRG')|((x.campus=='UCSB')&~x.zone.isin(['DEEP','SPECIALO']))]
# Keep raw rows only at their canonical source; export bounded derived summaries.
events=['site','year','month','zone','transect','campus'];taxa={'MESFRAAD':'M. franciscanus','STRPURAD':'S. purpuratus','PANINT':'P. interruptus'}
# Existing aggregation identifies events without day; audit any multiple dates rather than silently merging them.
daycheck=x.groupby(events).day.nunique();assert (daycheck<=1).all()
ev=x[events].drop_duplicates();grid=ev.merge(pd.DataFrame({'classcode':list(taxa)}),how='cross')
f=x[x.classcode.isin(taxa)].copy();assert not f['count'].isna().any();assert (f['count']>=0).all()
ct=f.groupby(events+['classcode'],dropna=False)['count'].sum().reset_index();g=grid.merge(ct,on=events+['classcode'],how='left',validate='one_to_one');g['explicit_zero']=g['count'].isna();g['count']=g['count'].fillna(0)
# Zero fill only focal taxa absent from an otherwise observed swath event, matching pipeline convention.
g['Taxa']=g.classcode.map(taxa);g['density']=g['count']/60;g.to_csv(out/'event_counts.csv',index=False)
survey=x.groupby(['site','year','month','zone','campus','method'],dropna=False).agg(rows=('count','size'),transects=('transect','nunique'),depth_min=('depth','min'),depth_max=('depth','max')).reset_index();survey['season_retained']=survey.month.between(5,10);survey.to_csv(out/'survey_coverage.csv',index=False)
def aggregate(z):
 a=z.groupby(['site','year','zone','transect','Taxa'],as_index=False).density.mean();b=a.groupby(['site','year','zone','Taxa'],as_index=False).density.mean();return b.groupby(['site','year','Taxa'],as_index=False).density.mean()
season=g[g.month.between(5,10)];annual=aggregate(season);annual.to_csv(out/'reconstructed_site_density.csv',index=False)
old=pd.read_csv('outputs/upstream_sites_20261101/pisco_selected_station_values.csv');old=old[(old.MPA=='Point Dume SMCA')&old.site.isin(sites)&(old.Resp=='Den')&old.Taxa.isin(taxa.values())]
j=annual.merge(old,on=['site','year','Taxa'],validate='one_to_one');err=max(abs(j.density-j.value));assert err<1e-10
post=season[season.year.between(2021,2023)];fixed={s:set.intersection(*[set(d.zone) for _,d in post[post.site==s].groupby('year')]) for s in sites}
allpost=g[g.year.between(2021,2023)];common=set.intersection(*[set(d.zone) for _,d in allpost.groupby(['site','year'])]);parts=[]
for name,z in [('Seasonal_available',post),('Seasonal_fixed_within_site',post[[r.zone in fixed[r.site] for _,r in post.iterrows()]]),('All_months_common_zones',allpost[allpost.zone.isin(common)])]:
 a=aggregate(z);a['scenario']=name;parts.append(a)
r=pd.concat(parts);r.to_csv(out/'footprint_site_densities.csv',index=False)
# Freeze original taxon-specific additive correction across every footprint diagnostic, including positive values.
w=r.pivot(index=['scenario','year','Taxa'],columns='site',values='density').reset_index();ledger=pd.read_csv('outputs/upstream_sites_20261101/candidate_eligibility.csv');floors=ledger[(ledger.MPA=='Point Dume SMCA')&(ledger.Resp=='Den')].set_index('Taxa').zero_correction;w['offset']=w.Taxa.map(floors);w['ln_ratio']=np.log((w.LECHUZA+w.offset)/(w.NICHOLAS_CANYON_W+w.offset));w.to_csv(out/'footprint_ratios.csv',index=False)
existing=pd.read_csv('outputs/researcher_scope_20261101/annual_input.csv');existing=existing[(existing.MPA=='Point Dume SMCA')&(existing.Source=='PISCO')&(existing.Resp=='Den')&existing.year.between(2021,2023)&existing.Taxa.isin(taxa.values())]
jj=w[w.scenario=='Seasonal_available'].merge(existing,on=['year','Taxa']);assert max(abs(jj.ln_ratio-jj.lnDiff))<1e-10
# Coverage definitions use no outcomes.
e=pd.read_csv('outputs/full_followup_recovery_20261101/effects.csv');e=e[e.scenario=='Full'].copy();e['span']=e.last_after-e.first_after
coverage=e[e.metric=='Post_slope'][['series_id','MPA','Source','Taxa','Resp','na','span','first_age','last_age','after_years']].copy();coverage['category']=np.where((coverage.na>=5)&(coverage.span>=5),np.where(coverage.first_age<=3,'Early-start multiyear trajectory','Later-start multiyear trajectory'),'Short-period trajectory');coverage.to_csv(out/'coverage_categories.csv',index=False)
parts=[];membership=[]
for n,span in [(4,3),(5,5),(7,7)]:
 for early in [None,3,5]:
  name=f'n{n}_span{span}_early'+('any' if early is None else str(early));mask=(e.na>=n)&(e.span>=span)&(True if early is None else e.first_age<=early);z=e[mask].copy();z['scenario']=name;parts.append(z)
  membership.append(coverage.assign(rule=name,included=(coverage.na>=n)&(coverage.span>=span)&(True if early is None else coverage.first_age<=early)))
base=e.copy();base['scenario']='All_three_year';parts.insert(0,base)
late=e[(e.na>=5)&(e.span>=5)&(e.first_age>3)].copy();late['scenario']='Late_only_n5_span5';parts.append(late)
short=e[(e.na<5)|(e.span<5)].copy();short['scenario']='Short_only';parts.append(short)
pd.concat(parts).to_csv(out/'effects.csv',index=False);pd.concat(membership).to_csv(out/'rule_membership.csv',index=False)
result={'raw_source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'max_site_density_error':err,'matched_site_year_taxa':len(j),'fixed_seasonal_zones':{k:sorted(v) for k,v in fixed.items()},'all_month_common_zones':sorted(common),'missing_raw_focal_counts':int(f['count'].isna().sum()),'event_dates_unique':True,'seasonal_ratio_reproduced':True}
(out/'source_verification.json').write_text(json.dumps(result,indent=2));print(result);print(coverage.category.value_counts().to_string())
