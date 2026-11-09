"""Review-only common-cohort reconstruction of the archived PISCO urchin size pool.

This compares point inputs. It does not claim to reproduce Emily's original full model.
"""
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import linregress
R=Path(__file__).resolve().parents[2]; B=R/'outputs/integrated_review_20261103';O=R/'outputs/emily_biomass_bridge_20261103';O.mkdir(exist_ok=True)
x=pd.read_parquet(B/'private/supported_measurement_uses.parquet')
x=x[(x.Source=='PISCO')&(x.Resp=='Bio')&x.Taxa.isin(['M. franciscanus','S. purpuratus'])].copy()
# Replay the archived size-file and physical-site crosswalk used in Emily's script.
A=Path('/Users/adrianstier/Library/CloudStorage/GoogleDrive-<account>/My Drive/Stier Lab/People/Emily Donham/Projects/Kelp MPA/Github/MPA_Kelp/data')
raw=pd.read_csv(A/'ALL_sizefreq_2024.csv',low_memory=False)
site=pd.read_csv(A/'PISCO/master_site_table_Emilyedit.csv',low_memory=False)
site=site.drop_duplicates(['site','CA_MPA_Name_Short'])
excluded_kfm={'a-k-21','a-k-05','a-k-07','a-k-26','a-k-27','a-k-28','a-k-29','a-k-30','a-k-36','a-k-37','a-k-12','a-k-13','a-k-31','a-k-35'}
excluded_ref={'ANACAPA_BLACK_SEA_BASS','SMI_BAY_POINT','SCAI_SHIP_ROCK','SMI_CROOK_POINT_E','SMI_CROOK_POINT_W','SBI_SUTIL','SCI_YELLOWBANKS_CEN','SCI_YELLOWBANKS_W','SRI_FORD_POINT','SMI_PRINCE_ISLAND_CEN','SMI_PRINCE_ISLAND_N','\tSMI_HARE_ROCK','SMI_TYLER_BIGHT_E','SMI_TYLER_BIGHT_W','SBI_GRAVEYARD_CANYON','SBI_GRAVEYARD_CANYON_N'}
raw=raw[~raw.site.isin(excluded_kfm|excluded_ref)].copy();raw['source_record_id']='emily_archive:'+raw.index.astype(str)
raw['Taxa']=raw.classcode.map({'MESFRA':'M. franciscanus','Strongylocentrotus franciscanus':'M. franciscanus','SFL':'M. franciscanus','SFS':'M. franciscanus','STRPUR':'S. purpuratus','Strongylocentrotus purpuratus':'S. purpuratus','SPL':'S. purpuratus','SPS':'S. purpuratus'})
raw=raw[raw.Taxa.notna()&raw['size'].ge(25)&raw['count'].gt(0)].copy()
p=raw.merge(site[['site','site_new','campus','CA_MPA_Name_Short','site_status']],on=['site','site_new','campus'],how='inner',validate='many_to_many')
p=p[p.CA_MPA_Name_Short.isin(x.series_id.str.split('|').str[0].unique())].copy()
p['mass']=np.where(p.Taxa.eq('M. franciscanus'),.00059*p['size']**2.917,.00059*p['size']**2.870)
p['weighted_mass']=p['mass']*p['count'];p['series_id']=p.CA_MPA_Name_Short+'|'+p.Taxa+'|PISCO|Bio';p=p.rename(columns={'site_status':'status'});
p['pool_id']='archive_pool|'+p.series_id+'|'+p.status+'|'+p.year.astype(str)
p[['pool_id','source_record_id','mass','count']].rename(columns={'count':'frequency'}).to_csv(O/'archived_donor_members.csv',index=False)
pool=p.groupby(['series_id','status','year','Taxa']).agg(total_mass=('weighted_mass','sum'),n_sized=('count','sum'),n_sites=('site','nunique'),n_pools=('site_new','nunique'),source_rows=('source_record_id','size')).reset_index()
pool['legacy_mean_mass']=pool.total_mass/pool.n_sized
pool.to_csv(O/'emily_style_pools.csv',index=False)
x=x.merge(pool[['series_id','status','year','Taxa','legacy_mean_mass']],on=['series_id','status','year','Taxa'],how='left',validate='many_to_one')
x['legacy_biomass']=x.density*x.legacy_mean_mass
x.loc[x.density.eq(0),'legacy_biomass']=0.
# Site -> zone -> date/method -> unit, same weights as current integrated estimate.
def agg(g,col):
    if g[col].isna().any():return np.nan
    return float(g[col].mean())
def annual(g,col):
    vals=[]
    for _,site in g.groupby('site'):
        zones=[]
        for _,zone in site.groupby('zone'):
            visits=[]
            z=zone.assign(method=zone.unit.str.split('|').str[0])
            for _,event in z.groupby(['date','method']):visits.append(agg(event,col))
            zones.append(np.mean(visits))
        vals.append(np.mean(zones))
    return np.mean(vals)
rows=[]
for (sid,yr),g in x.groupby(['series_id','year']):
    rec=dict(series_id=sid,year=int(yr),Taxa=g.Taxa.iloc[0],MPA_Start=int(g.MPA_Start.iloc[0]))
    for status,sg in g.groupby('status'):
        for col,name in [('biomass','current'),('legacy_biomass','emily_style')]:rec[status+'_'+name]=annual(sg,col)
    rows.append(rec)
out=pd.DataFrame(rows)
original=pd.read_csv(B/'annual_values.csv');original=original[original.series_id.isin(x.series_id.unique())][['series_id','year','protected','reference','offset','response_supported']]
out=out.merge(original,on=['series_id','year'],validate='one_to_one')
assert np.nanmax(abs(out.mpa_current-out.protected))<1e-8
assert np.nanmax(abs(out.reference_current-out.reference))<1e-8
out['matched']=out[['mpa_current','reference_current','mpa_emily_style','reference_emily_style']].notna().all(axis=1)
for v in ['current','emily_style']:
 out['log_ratio_'+v]=np.log((out['mpa_'+v]+out.offset)/(out['reference_'+v]+out.offset))
out['size_matching_delta_log_ratio']=out.log_ratio_emily_style-out.log_ratio_current
out.to_csv(O/'paired_annual_values.csv',index=False)
q=out[out.matched].copy();assert len(q)
res=[]
for sid,g in q.groupby('series_id'):
 g=g[g.year>=g.MPA_Start.iloc[0]].sort_values('year')
 if len(g)<3:continue
 r=dict(series_id=sid,taxon=g.Taxa.iloc[0],matched_post_years=len(g),first_year=int(g.year.min()),last_year=int(g.year.max()))
 for v in ['current','emily_style']:
  fit=linregress(g.year,g['log_ratio_'+v]);r[v+'_annual_pct']=100*(np.exp(fit.slope)-1);r[v+'_slope_log']=fit.slope
 r['size_matching_slope_delta_pct_points']=r['emily_style_annual_pct']-r['current_annual_pct'];res.append(r)
pd.DataFrame(res).to_csv(O/'same_year_descriptive_slopes.csv',index=False)
print('paired years',len(q),'of',len(out),'histories',q.series_id.nunique(),'post slopes',len(res));print(pd.DataFrame(res).groupby('taxon')[['current_annual_pct','emily_style_annual_pct','size_matching_slope_delta_pct_points']].agg(['mean','median']).to_string())
# A transparent 2 x 2 diagnostic: fixed matched years and annual ratios,
# with independent-year versus lag-1 HAC uncertainty. This is NOT Emily's
# archived regression or the integrated joint covariance model.
import statsmodels.api as sm
unc=[]
for sid,g in q.groupby('series_id'):
 g=g[g.year>=g.MPA_Start.iloc[0]].sort_values('year')
 if len(g)<4:continue
 X=sm.add_constant(g.year.to_numpy(float)-g.year.min())
 for v in ['current','emily_style']:
  y=g['log_ratio_'+v].to_numpy(float)
  for method in ['independent_years','HAC_lag1']:
   f=sm.OLS(y,X).fit(cov_type='HAC',cov_kwds={'maxlags':1}) if method=='HAC_lag1' else sm.OLS(y,X).fit()
   lo,hi=f.conf_int()[1];unc.append(dict(series_id=sid,taxon=g.Taxa.iloc[0],n_years=len(g),size_pool=v,uncertainty=method,annual_pct=100*np.expm1(f.params[1]),low_pct=100*np.expm1(lo),high_pct=100*np.expm1(hi),slope_se=f.bse[1]))
pd.DataFrame(unc).to_csv(O/'factorial_diagnostic.csv',index=False)
print('factorial fits',len(unc),'on',len(set(v['series_id'] for v in unc)),'series')
import hashlib,json
(O/'archive_input_hashes.json').write_text(json.dumps({str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in [A/'ALL_sizefreq_2024.csv',A/'PISCO/master_site_table_Emilyedit.csv']},indent=2))
