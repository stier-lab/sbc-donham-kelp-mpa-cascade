from pathlib import Path
import pandas as pd,numpy as np,ast,json
from scipy.stats import t as student
O=Path('outputs/biomass_protocols_20261102');D=Path('/Users/adrianstier/Donham-Stier-CA-MPA-Data-2026/data')
tree=ast.parse(Path('code/tools/audit_station_coverage_refits.py').read_text());exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['fit','slope']],type_ignores=[]),'fit','exec'))
a=pd.read_csv('outputs/biomass_resolution_20261102/review_annual_input.csv');i=pd.read_csv('outputs/biomass_support_20261102/animal_pool_inventory.csv');off=pd.read_csv('outputs/all_program_station_refits_20261102/offsets.csv').set_index('series_id').offset
l=pd.read_csv(D/'source_archive/kfm_size_feasibility_2026-10-29/legacy/Invertebrate Size Frequency.csv');p=pd.read_csv('outputs/biomass_support_20261102/private/pool_members.csv',low_memory=False);h=pd.read_csv(O/'kfm_label_hypotheses.csv');h=h[h.candidate.eq('Southeast Reef')];assert len(h)==14 and h.exact_match.all()
rows=[];fits=[];corrected=a.copy()
for tax,code,power in [('M. franciscanus',11005,2.917),('S. purpuratus',11006,2.870)]:
 sid='Santa Barbara Island SMR|'+tax+'|KFM|Bio';z=a[a.series_id.eq(sid)&a.t.gt(0)].copy();z['new']=z.lnDiff
 for year in range(2005,2012):
  q=l[l.Site.eq(14)&l.Year.eq(year)&l.Species.eq(code)];q37=l[l.Site.eq(37)&l.Year.eq(year)&l.Species.eq(code)];old=p[p.series_id.eq(sid)&p.year.eq(year)&p.status.eq('mpa')]
  assert len(q)>0 and set(old.site)=={'a-k-14'}
  oldmass=np.average(.00059*old['size']**power,weights=old['count']);newmass=np.mean(.00059*(q.Size*10)**power);factor=newmass/oldmass
  vv=i[i.series_id.eq(sid)&i.year.eq(year)].set_index('status').value;new=np.log((vv['mpa']*factor+off[sid])/(vv['reference']+off[sid]));z.loc[z.year.eq(year),'new']=new;corrected.loc[corrected.series_id.eq(sid)&corrected.year.eq(year),'lnDiff']=new
  rows.append(dict(series_id=sid,year=year,old_size_site='Southeast Reef (legacy37)',count_and_replacement_size_site='SE Sea Lion Rookery (legacy14)',old_n=int(old['count'].sum()),replacement_n=len(q),old_mean_mass=oldmass,replacement_mean_mass=newmass,biomass_factor=factor,old_lnDiff=float(z.loc[z.year.eq(year),'lnDiff'].iloc[0]),new_lnDiff=new))
 for scope,zz in [('all_observed_post_years',z),('2005_2011_same_years',z[z.year.between(2005,2011)])]:
  for label,col in [('current','lnDiff'),('original_count_site_sizes','new')]:fits.append(dict(series_id=sid,scope=scope,scenario=label,**slope(zz,col)))
pd.DataFrame(rows).to_csv(O/'sbi_size_replacements.csv',index=False);pd.DataFrame(fits).to_csv(O/'sbi_reconciliation_effects.csv',index=False);corrected.to_csv(O/'sbi_reconciled_review_variant.csv',index=False)
# Verify count-site identity independently against original legacy densities.
c=pd.read_csv(D/'MBON/SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv',low_memory=False);c=c[c.site_id.eq('a-k-14')&c.taxon_name.isin(['Mesocentrotus franciscanus','Strongylocentrotus purpuratus'])].copy();c['year']=pd.to_datetime(c.date).dt.year;c['den']=pd.to_numeric(c['count'],errors='coerce')/c.area
ld=pd.read_csv(D/'source_archive/kfm_size_feasibility_2026-10-29/legacy/Benthic Density Data.csv');checks=[]
for name,code in [('Mesocentrotus franciscanus',11005),('Strongylocentrotus purpuratus',11006)]:
 for year in range(2005,2012):
  actual=c[c.taxon_name.eq(name)&c.year.eq(year)].den.mean();v=ld[ld.Site.eq(14)&ld.Year.eq(year)&ld.Species.eq(code)].DensityMean.iloc[0];assert abs(actual-v)<1e-9;checks.append(dict(species=name,year=year,mbon_mean_density=actual,legacy14_density=v))
pd.DataFrame(checks).to_csv(O/'sbi_count_identity_checks.csv',index=False)
print(pd.DataFrame(fits).to_string(index=False))
