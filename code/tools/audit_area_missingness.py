from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
D=Path('/Users/adrianstier/sbc-donham-kelp-mpa-data');O=Path('outputs/required_scope_area_20261102');O.mkdir(exist_ok=True)
p=D/'data/MBON/SBCMBON_kelp_forest_integrated_quad_and_swath_20231022.csv';raw=pd.read_csv(p,low_memory=False);x=raw[raw.data_source.eq('lter')&raw.taxon_name.isin(['Strongylocentrotus purpuratus','Mesocentrotus franciscanus'])].copy();x['count']=pd.to_numeric(x['count'],errors='coerce');x['area']=pd.to_numeric(x.area,errors='coerce');x['year']=pd.to_datetime(x.date).dt.year
x.groupby(['sample_method','area'],dropna=False).agg(records=('count','size'),missing_count=('count',lambda x:x.isna().sum())).reset_index().to_csv(O/'lter_source_area.csv',index=False)
assert x.area.eq(1).all() and x['count'].notna().all() and x['count'].ge(0).all()
c=pd.read_csv(O/'lter_legacy_bootstrap.csv');missing=c[c.biomass.isna()].copy();missing.to_csv(O/'lter_missing_biomass.csv',index=False)
# Reproduce old zero-then-remove rule and compare explicit missing support without imputing zero.
keys=['site','year','y'];old=c.copy();old=old.fillna({k:0 for k in old.select_dtypes('number').columns});bad=old[old['count'].gt(0)&old.biomass.eq(0)][keys].drop_duplicates();keep=old.merge(bad.assign(drop=True),on=keys,how='left');keep=keep[keep['drop'].isna()].drop(columns='drop')
new=c.copy();bad2=new[new['count'].gt(0)&(new.biomass.isna()|new.biomass.le(0))][keys].drop_duplicates();keep2=new.merge(bad2.assign(drop=True),on=keys,how='left');keep2=keep2[keep2['drop'].isna()].drop(columns='drop');pd.testing.assert_frame_equal(keep.drop(columns=['site_designation','BaselineRegion']).reset_index(drop=True),keep2.drop(columns=['site_designation','BaselineRegion']).reset_index(drop=True))
keep2['old_g_per_quad']=keep2.biomass/keep2.quad;keep2['explicit_g_m2']=keep2.biomass/(keep2.quad*keep2.area);assert np.allclose(keep2.old_g_per_quad,keep2.explicit_g_m2)
keep2[['site','year','y','area','quad','old_g_per_quad','explicit_g_m2']].to_csv(O/'lter_area_comparison.csv',index=False)
# Check dates and quadrat identities within physical samples before metadata joins.
physical=['site_id','transect_id','taxon_name','year'];visits=x.groupby(physical).agg(dates=('date','nunique'),rows=('count','size'),replicates=('replicate_id','nunique')).reset_index();visits.to_csv(O/'lter_visit_counts.csv',index=False)
geo=pd.read_csv(D/'data/MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv');g=geo[geo.CA_MPA_Name_Short.isin(['Campus Point SMCA','Naples SMCA'])&~geo.site_id.eq('a-l-02')];dup=g.groupby('site_id').CA_MPA_Name_Short.nunique();dup.to_csv(O/'lter_comparison_memberships.csv')
report=dict(lter_source_urchin_rows=len(x),all_areas_one_m2=True,missing_source_counts=int(x['count'].isna().sum()),cached_transect_rows=len(c),missing_biomass_rows=len(missing),affected_site_year_species=len(bad2),legacy_and_explicit_missing_selection_identical=True,maximum_area_conversion_change=float((keep2.old_g_per_quad-keep2.explicit_g_m2).abs().max()),multi_date_transect_taxon_years=int(visits.dates.gt(1).sum()),physical_sites_in_both_lter_comparisons=dup[dup.gt(1)].index.tolist(),source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),scope='LTER source areas and cached aggregation; not a fresh biomass/model refit')
(O/'area_missing_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
# Wider source inventory: keep this distinct from the numerically checked LTER branch.
focal=['Strongylocentrotus purpuratus','Mesocentrotus franciscanus','Macrocystis pyrifera','Panulirus interruptus']
z=raw[raw.taxon_name.isin(focal)].copy();z['count_num']=pd.to_numeric(z['count'],errors='coerce');z['area_num']=pd.to_numeric(z.area,errors='coerce');z['invalid_area']=z.area_num.isna()|z.area_num.le(0)
z.groupby(['data_source','sample_method','taxon_name']).agg(rows=('count_num','size'),missing_count=('count_num',lambda v:v.isna().sum()),invalid_area=('invalid_area','sum'),area_values=('area_num',lambda v:';'.join(map(str,sorted(v.dropna().unique()))))).reset_index().to_csv(O/'mbon_area_missing_inventory.csv',index=False)
z[z.data_source.eq('lter')&z.taxon_name.eq('Macrocystis pyrifera')&z.sample_method.eq('swath')&z.count_num.isna()].to_csv(O/'lter_kelp_missing_counts.csv',index=False)
y=pd.read_csv(D/'data/PISCO/MLPA_kelpforest_swath_2024.csv',encoding='cp1252',low_memory=False);f=y[y.classcode.isin(['MESFRAAD','STRPURAD','PANINT','MACPYRAD'])]
f.groupby(['campus','classcode']).agg(rows=('count','size'),missing_count=('count',lambda v:v.isna().sum()),missing_size=('size',lambda v:v.isna().sum())).reset_index().to_csv(O/'pisco_swath_missing_inventory.csv',index=False)
f[f.campus.isin(['UCSB','VRG'])&((f.classcode.eq('MACPYRAD')&f['size'].isna())|(f.classcode.eq('MESFRAAD')&f['count'].isna()))].to_csv(O/'pisco_missing_focal_values.csv',index=False)
# Pelican original-record inspection; retain both versions without choosing on effect direction.
pel=D/'data/source_archive/kfm_size_feasibility_2026-10-29/legacy/Invertebrate Size Frequency.csv';z=pd.read_csv(pel);z['source_row']=range(1,len(z)+1);z=z[z.Site.eq(8)&z.Year.eq(1995)&z.Species.eq(11005)];z.to_csv(O/'pelican_original_rows.csv',index=False)
c=pd.read_csv('data/ALL_sizefreq_2024.csv',low_memory=False);c=c[c.campus.eq('KFM')&c.site.eq('a-k-08')&c.year.eq(1995)&c.classcode.eq('Strongylocentrotus franciscanus')];a=z.groupby('Size').size();freq=pd.DataFrame({'original_n':a.set_axis((a.index*10).round(6)),'combined_n':c.groupby('size')['count'].sum()}).fillna(0);freq['difference']=freq.original_n-freq.combined_n;assert freq.difference.ge(0).all() and freq.difference.sum()==5;freq.to_csv(O/'pelican_size_frequencies.csv')
inputs=[p,pel,D/'data/cache/lter_urchin_bootstrap.rds',Path('data/ALL_sizefreq_2024.csv'),D/'config/deferred_kfm_source_groups.csv']
(O/'input_hashes.json').write_text(json.dumps({str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},indent=2)+'\n')
