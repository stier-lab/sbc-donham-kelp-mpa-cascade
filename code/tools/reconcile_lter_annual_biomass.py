"""Compare current raw kelp reconstruction to pinned SBC50.21; no promotion."""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
A=Path(__file__).resolve().parents[2];D=A.parent/'sbc-donham-kelp-mpa-data';P=A.parent/'sbc-donham-kelp-mpa-data/data/source_archive/fish_lter_reconciliation_20261102';O=A/'outputs/fish_lter_reconciliation_20261102'
registry={s['source_id']:s for s in json.loads((D/'config/source_manifest.json').read_text())['sources']}
raw=pd.read_csv(D/registry['lter_kelp']['path']); slopes=pd.read_csv(D/registry['kelp_calibration']['path']);slopes=slopes[slopes.model.eq(1)].set_index('MONTH').slope
r=raw[raw.FRONDS.ne(-99999)].copy();assert r.AREA.gt(0).all();r['current_biomass']=r.FRONDS/r.AREA*r.MONTH.map(slopes)*1000
q=r.groupby(['YEAR','MONTH','DATE','SITE','TRANSECT','QUAD','SIDE','AREA'],dropna=False).agg(fronds=('FRONDS','sum'),current_biomass=('current_biomass','sum')).reset_index();q['current_density']=q.fronds/q.AREA
keys=['YEAR','SITE','TRANSECT']; a=q.groupby(keys).agg(current_density=('current_density','mean'),current_biomass=('current_biomass','mean')).reset_index()
p=P/'Annual_All_Species_Biomass_at_transect_20261012.csv';d=pd.read_csv(p);d=d[d.SP_CODE.eq('MAPY')].copy();assert not d.duplicated(keys).any();assert not a.duplicated(keys).any()
d[['DRY_GM2','DENSITY']]=d[['DRY_GM2','DENSITY']].replace(-99999,np.nan)
b=d.rename(columns={'DRY_GM2':'provider_dry_g_m2','DENSITY':'provider_density'})
z=a.merge(b[keys+['MONTH','DATE','provider_dry_g_m2','provider_density']],on=keys,how='outer',indicator=True,validate='one_to_one');z['density_delta']=z.provider_density-z.current_density;z['biomass_delta']=z.provider_dry_g_m2-z.current_biomass;z.to_csv(O/'lter_transect_comparison.csv',index=False)
oldpath=Path('/Volumes/Stier_Lab/sbc-lter-data-catalogue/raw/packages/knb-lter-sbc/50/20/entities/24d18d9ebe4f6e8b94e222840096963c/Annual_All_Species_Biomass_at_transect_20260928.csv');old=pd.read_csv(oldpath);old=old[old.SP_CODE.eq('MAPY')].copy();old['DRY_GM2']=old.DRY_GM2.replace(-99999,np.nan);v=old[keys+['DRY_GM2']].merge(b[keys+['provider_dry_g_m2']],on=keys,how='outer',indicator=True);v['difference']=v.provider_dry_g_m2-v.DRY_GM2;v.to_csv(O/'lter_provider_revision_comparison.csv',index=False)
# Current source mapping and ecological scope, maintaining identical years/transects.
m=pd.read_csv(A/'data/MBON/SBCMBON_kelp_forest_site_geolocation_20210120_KFM_LTER.csv');m=m[m.data_source.eq('lter')];print(m[['site_id','CA_MPA_Name_Short','status']].to_string(index=False))
j=z[z._merge.eq('both')].merge(m[['site_id','CA_MPA_Name_Short','status']],left_on='SITE',right_on='site_id',validate='many_to_one'); j=j[j.CA_MPA_Name_Short.isin(['Campus Point SMCA','Naples SMCA']) & ~j.SITE.eq('AHND')]
j.to_csv(O/'lter_mapped_comparison.csv',index=False)
site=j.groupby(['CA_MPA_Name_Short','status','YEAR','SITE'])[['current_biomass','provider_dry_g_m2']].mean().reset_index(); role=site.groupby(['CA_MPA_Name_Short','status','YEAR'])[['current_biomass','provider_dry_g_m2']].mean().reset_index(); role.to_csv(O/'lter_role_year_comparison.csv',index=False)
summary=dict(package='knb-lter-sbc.50.21',raw_rows=len(raw),current_transects=len(a),provider_mapy_transects=len(d),matched_transects=int(z._merge.eq('both').sum()),current_only=int(z._merge.eq('left_only').sum()),provider_only=int(z._merge.eq('right_only').sum()),max_abs_density_delta=float(z.density_delta.abs().max()),max_abs_biomass_delta=float(z.biomass_delta.abs().max()),biomass_changed_above_1e_minus6=int(z.biomass_delta.abs().gt(1e-6).sum()),mapped_matched_transects=len(j),mapped_max_abs_biomass_delta=float(j.biomass_delta.abs().max()),older_release_max_abs_difference=float(v.difference.abs().max()),older_release_changed=int(v.difference.abs().gt(1e-6).sum()),new_file_sha256=hashlib.sha256(p.read_bytes()).hexdigest())
(O/'lter_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

# Compare focal taxa across provider releases independently of our reconstruction.
old_all=pd.read_csv(oldpath); new_all=pd.read_csv(p); revisions=[]
for species,response in [('MAPY','DRY_GM2'),('PAIN','DRY_GM2'),('SFL','DRY_GM2'),('SPL','DRY_GM2'),('SPUL','WM_GM2')]:
    left=old_all.loc[old_all.SP_CODE.eq(species),keys+[response]].replace(-99999,np.nan)
    right=new_all.loc[new_all.SP_CODE.eq(species),keys+[response]].replace(-99999,np.nan)
    joined=left.merge(right,on=keys,how='outer',suffixes=('_old','_new'),indicator=True,validate='one_to_one')
    delta=joined[response+'_new']-joined[response+'_old']
    revisions.append(dict(species_code=species,response=response,matched_rows=int(joined._merge.eq('both').sum()),old_only=int(joined._merge.eq('left_only').sum()),new_only=int(joined._merge.eq('right_only').sum()),changed_finite_rows=int(delta.abs().gt(1e-6).sum()),max_absolute_change=float(delta.abs().max()),scope='provider21 versus provider20; not a replacement for this papers raw-based estimates'))
pd.DataFrame(revisions).to_csv(O/'lter_provider_focal_taxa_revisions.csv',index=False)
