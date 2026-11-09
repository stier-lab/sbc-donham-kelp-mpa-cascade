"""Read-only source reconciliation; current count roster, no biomass promotion."""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
from collections import Counter
A=Path(__file__).resolve().parents[2]; D=A.parent/'sbc-donham-kelp-mpa-data'; P=A.parent/'sbc-donham-kelp-mpa-data/data/source_archive/fish_lter_reconciliation_20261102'; O=A/'outputs/fish_lter_reconciliation_20261102'
p=P/'chis_kelp_forest_fish_size.csv'; d=pd.read_csv(p); s=d[d.scientificName.str.contains('pulcher',na=False)].copy();s['source_row']=s.index+2
mapping=pd.read_csv(D/'config/kfm_physical_sites.csv').set_index('source_name').site
s['site_id']=s.siteName.map(mapping); assert s.site_id.notna().all()
s['year']=s.startDate.str[:4].astype(int);s['stage']=s.speciesRemarks.str.lower();s['length']=s.totalLengthInCentimeters
assert (s.speciesCount>0).all() and (s.length>0).all() and (s.speciesCount%1==0).all()
assert (s.startDate==s.endDate).all()
s['mass']=.0144*s.length**3.04
s['lower_mass']=.0144*s.minimumLengthInCentimeters.fillna(s.length)**3.04
s['upper_mass']=.0144*s.maximumLengthInCentimeters.fillna(s.length)**3.04
s['range_present']=s.minimumLengthInCentimeters.notna()&s.maximumLengthInCentimeters.notna()
# Legacy multiset comparison preserves repeated equal-sized fish.
old=pd.read_csv(D/'data/source_archive/kfm_size_feasibility_2026-10-29/legacy/Fish size frequency data.csv',encoding='cp1252'); old=old[old.ScientificName.str.contains('pulcher',na=False)].copy()
old['site_id']=old.Site.map(lambda x:f'a-k-{int(x):02}');old['date']=pd.to_datetime(old.Date,format='%d%b%Y').dt.strftime('%Y-%m-%d');old['stage']=old.ScientificName.str.split(',').str[-1].str.strip().str.lower()
a=Counter(zip(old.site_id,old.date,old.stage,old.Count,old.TotalLength))
b=Counter(zip(s.site_id,s.startDate,s.stage,s.speciesCount,s.length));dif=[]
for key in sorted(a.keys()|{k for k in b if int(k[1][:4])<=2011}):
 if a[key]!=b[key]:dif.append(dict(site_id=key[0],date=key[1],stage=key[2],count=key[3],length_cm=key[4],old_rows=a[key],new_rows=b[key]))
pd.DataFrame(dif,columns=['site_id','date','stage','count','length_cm','old_rows','new_rows']).to_csv(O/'fish_legacy_differences.csv',index=False)
pools=[]
for (site,year),g in s.groupby(['site_id','year']):
 n=g.speciesCount.sum();pools.append(dict(site_id=site,year=year,new_n_sized=n,new_n_rows=len(g),new_mean_mass_g=np.average(g.mass,weights=g.speciesCount),new_lower_mass_g=np.average(g.lower_mass,weights=g.speciesCount),new_upper_mass_g=np.average(g.upper_mass,weights=g.speciesCount),size_dates=';'.join(sorted(g.startDate.unique())),n_range_records=g.range_present.sum(),n_juveniles=g.loc[g.stage.eq('juvenile'),'speciesCount'].sum()))
pools=pd.DataFrame(pools);pools.to_csv(O/'fish_size_pools.csv',index=False)
legacy=pd.read_csv(A/'tables/table_s_kfm_size_pools.csv'); comp=legacy.merge(pools,on=['site_id','year'],how='outer',indicator=True);comp['mass_difference_g']=comp.new_mean_mass_g-comp.mean_mass_g;comp.to_csv(O/'fish_pool_comparison.csv',index=False)
c=pd.read_csv(A/'tables/table_s_kfm_size_site_year_matches.csv');c=c[c.species=='S. pulcher'].copy(); c=c.merge(pools,on=['site_id','year'],how='left',validate='many_to_one');c['new_has_sizes']=c.new_n_sized.notna()
c['all_size_dates_match']=c.apply(lambda r:set(str(r.size_dates).split(';'))<=set(r.count_dates.split(';')) if r.new_has_sizes else False,axis=1)
c['new_supported_point']=c.new_has_sizes|c.count_density.eq(0);c['old_supported_point']=c.has_size_pool|c.count_density.eq(0)
c['provisional_biomass_g_m2']=np.where(c.count_density.eq(0),0,c.count_density*c.new_mean_mass_g)
c.to_csv(O/'fish_count_size_matches.csv',index=False)
ann=c.groupby(['MPA','year']).agg(n_sites=('site_id','size'),size_sites=('new_has_sizes','sum'),supported_sites=('new_supported_point','sum'),old_supported_sites=('old_supported_point','sum'),roles=('status','nunique')).reset_index();ann['complete']=ann.supported_sites.eq(ann.n_sites)&ann.roles.eq(2);ann['old_complete']=ann.old_supported_sites.eq(ann.n_sites)&ann.roles.eq(2);ann.to_csv(O/'fish_comparison_year_support.csv',index=False)
summary=c.groupby('MPA').agg(count_site_years=('site_id','size'),old_size_site_years=('has_size_pool','sum'),new_size_site_years=('new_has_sizes','sum'),same_date_site_years=('all_size_dates_match','sum'),new_sized_fish=('new_n_sized','sum'));summary['complete_comparison_years']=ann.groupby('MPA').complete.sum();summary['old_complete_comparison_years']=ann.groupby('MPA').old_complete.sum();summary.to_csv(O/'fish_coverage_summary.csv')
result=dict(file_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),all_species_rows=len(d),sheephead_rows=len(s),sheephead_fish=int(s.speciesCount.sum()),first_date=s.startDate.min(),last_date=s.startDate.max(),legacy_rows=len(old),legacy_difference_keys=len(dif),current_count_site_years=len(c),old_matched_site_years=int(c.has_size_pool.sum()),new_matched_site_years=int(c.new_has_sizes.sum()),new_all_dates_match=int(c.all_size_dates_match.sum()),old_complete_pair_years=int(ann.old_complete.sum()),new_complete_pair_years=int(ann.complete.sum()),positive_count_without_size=int((c.count_density.gt(0)&~c.new_has_sizes).sum()),matched_pools_under10=int((c.new_n_sized<10).sum()),pre2003_size_records=int((s.year<2003).sum()),notes='Provisional site-year mass only; roving size versus visual count footprint remains unverified. No primary biomass or new exclusions adopted.')
(O/'fish_summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));print(summary.to_string())
# Separate exact-date pools from annual pooling; repeat visits receive equal weight.
visit=[]
for _,r in c.iterrows():
 g=s[(s.site_id==r.site_id)&(s.year==r.year)]; target=set(r.count_dates.split(';'));g=g[g.startDate.isin(target)]
 masses=[np.average(v.mass,weights=v.speciesCount) for _,v in g.groupby('startDate')]
 visit.append(dict(MPA=r.MPA,status=r.status,site_id=r.site_id,year=int(r.year),count_dates=r.count_dates,n_exact_date_fish=int(g.speciesCount.sum()),n_size_visits=len(masses),exact_visit_mean_mass_g=float(np.mean(masses)) if masses else np.nan,coverage_status='count_zero_no_sizes_needed_for_point' if r.count_density==0 else 'exact_date_sizes' if masses else 'year_only_sizes' if r.new_has_sizes else 'no_sizes',biomass_point_g_m2=0 if r.count_density==0 else r.count_density*np.mean(masses) if masses else np.nan))
pd.DataFrame(visit).to_csv(O/'fish_exact_visit_support.csv',index=False)
