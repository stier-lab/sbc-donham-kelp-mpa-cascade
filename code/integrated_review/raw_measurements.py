"""Source rows -> physical survey units; unknown counts never become zeros.
No frozen analysis values enter these calculations. Site inclusion is applied
later through the approved comparison roster, and physical exclusions here.
"""
from common import *
frames=[];decisions=[]
cal=read_source('kelp_calibration');cal=cal[(cal.SP_CODE=='MAPY')&(cal.model==1)].set_index('MONTH').slope

def finish(x,source,source_id):
 x=x.copy();x['Source']=source;x['source_id']=source_id
 x=x.loc[~pd.Series([blocked(source,s,t) for s,t in zip(x.site,x.Taxa)],index=x.index,dtype=bool)]
 if x.empty:return
 x['year']=pd.to_datetime(x.date,errors='raise').dt.year
 x=x[x.year<=2023].copy() # registered manuscript source horizon; later provider years not adopted
 assert x.area.gt(0).all();assert not x.density.dropna().lt(0).any();assert not x.biomass.dropna().lt(0).any()
 for col,default in [('zone','not_recorded'),('size_basis','unsupported'),('donor_pool',''),('bio_mean_size_proxy',False),('imputed_density',0.),('imputed_pool','')]:
  if col not in x:x[col]=default
 x['unit_id']=[sid(source,s,d,z,u) for s,d,z,u in x[['site','date','zone','unit']].itertuples(index=False,name=None)]
 # Units may contain different species but cannot repeat within a species.
 assert not x.duplicated(['unit_id','Taxa']).any(),x[x.duplicated(['unit_id','Taxa'],False)].head().to_dict('records')
 frames.append(x[['Source','source_id','source_rows','site','year','date','zone','unit','unit_id','Taxa','area','density','biomass','size_basis','donor_pool','bio_mean_size_proxy','imputed_density','imputed_pool']])

def pisco():
 p=DATA/'data/processed/source_build/pisco_swath_reconciled.csv';sw=read_tracked(p,encoding='utf-8-sig');sw['_row']=np.arange(1,len(sw)+1)
 sw=sw[sw.campus.eq('VRG')|(sw.campus.eq('UCSB')&~sw.zone.isin(['DEEP','SPECIALO']))].copy()
 sw['date']=pd.to_datetime(sw[['year','month','day']]).dt.strftime('%Y-%m-%d');sw['unit']=sw.method.astype(str)+'|'+sw.transect.astype(str)
 key=['site','date','zone','unit'];roster=sw[key].drop_duplicates()
 for code in ['MESFRAAD','STRPURAD','PANINT','MACPYRAD']:
  f=sw[sw.classcode.eq(code)].copy();f['density']=pd.to_numeric(f['count'],errors='coerce')/60
  f['biomass']=np.nan;f['size_basis']='unsupported'
  if code=='MACPYRAD':
   f['biomass']=f['count']*f['size']*f.month.map(cal)*1000/60;f['size_basis']='same_record_stipes'
  elif code=='PANINT':
   # UCSB swaths contain individual/group measured lengths; VRG uses roving sizes.
   ok=f.campus.eq('UCSB')&f['size'].gt(0);f.loc[ok,'biomass']=f.loc[ok,'count']*mass(TAXA[code],f.loc[ok,'size'])/60;f.loc[ok,'size_basis']='same_record_length'
  f.loc[f['count'].eq(0),'biomass']=0
  g=f.groupby(key,dropna=False).agg(density=('density',strict_sum),biomass=('biomass',strict_sum),source_rows=('_row',lambda s:';'.join(map(str,s))),size_basis=('size_basis',lambda s:';'.join(sorted(set(s))))).reset_index()
  x=roster.merge(g,on=key,how='left',indicator=True,validate='one_to_one');absent=x._merge.eq('left_only');x.loc[absent,['density','biomass']]=0;x.loc[absent,'size_basis']='surveyed_zero';x['Taxa']=TAXA[code];x['area']=60
  finish(x,'PISCO','pisco_swath_reconciled')
 fi=read_source('pisco_fish',encoding='latin1');fi['_row']=np.arange(1,len(fi)+1);fi=fi[fi.campus.isin(['UCSB','VRG'])&fi.level.ne('CAN')].copy();fi['date']=pd.to_datetime(fi[['year','month','day']]).dt.strftime('%Y-%m-%d');fi['unit']=fi.method.astype(str)+'|'+fi.transect.astype(str)
 r=fi.groupby(key,dropna=False)['level'].agg(lambda s:set(s)).reset_index();r['complete']=r.level.map(lambda s:'BOT'in s and bool({'MID','CNMD'}&s));r.drop(columns='level').to_csv(OUT/'pisco_pass_decisions.csv',index=False)
 f=fi[fi.classcode=='SPUL'].copy();f['density']=f['count']/60;f['biomass']=f['count']*mass('S. pulcher',f.fish_tl)/60;f.loc[f.fish_tl.le(0)|f.fish_tl.isna(),'biomass']=np.nan;f.loc[f['count'].eq(0),'biomass']=0
 g=f.groupby(key,dropna=False).agg(density=('density',strict_sum),biomass=('biomass',strict_sum),source_rows=('_row',lambda s:';'.join(map(str,s)))).reset_index();x=r[r.complete].merge(g,on=key,how='left',indicator=True,validate='one_to_one');x.loc[x._merge.eq('left_only'),['density','biomass']]=0;x['size_basis']='same_record_fish_length';x['Taxa']='S. pulcher';x['area']=60;finish(x,'PISCO','pisco_fish')

def mbon():
 x=read_source('mbon_counts',na_values=['.']);x['_row']=np.arange(1,len(x)+1);x['Taxa']=x.taxon_name.str.strip().map(TAXA);x=x[x.data_source.eq('kfm')&x.Taxa.notna()].copy();x=x[~x.proj_taxon_id.eq('t-k-002')];x=x[x.Taxa.ne('M. pyrifera')|x.area.eq(10)]
 x['site']=x.site_id;x['unit']=x.sample_method.astype(str)+'|'+x.transect_id.astype(str)+'|'+x.replicate_id.astype(str);x['date']=pd.to_datetime(x.date).dt.strftime('%Y-%m-%d');x['density']=x['count']/x.area;x['biomass']=np.nan;x.loc[x['count'].eq(0),'biomass']=0
 # Distinct taxon IDs represent additive size categories, but never add repeat dates.
 key=['site','date','unit','Taxa','area'];x=x.groupby(key,dropna=False).agg(density=('density',strict_sum),biomass=('biomass',strict_sum),source_rows=('_row',lambda s:';'.join(map(str,s)))).reset_index();finish(x,'KFM','mbon_counts')
 f=read_source('mbon_fish',na_values=['.']);f['_row']=np.arange(1,len(f)+1);f=f[f.data_source.eq('kfm')&f.sample_method.eq('visualfish')&f.taxon_name.str.strip().isin(['Semicossyphus pulcher','Bodianus pulcher'])].copy();f['site']=f.site_id;f['Taxa']='S. pulcher';f['unit']=f.transect_id.astype(str)+'|'+f.replicate_id.astype(str);f['date']=pd.to_datetime(f.date).dt.strftime('%Y-%m-%d');f['density']=pd.to_numeric(f['count'],errors='coerce')/f.area;f['biomass']=np.nan;f.loc[f.density.eq(0),'biomass']=0;f['source_rows']=f['_row'].astype(str)
 # Several supplied fish records may share a transect; retain each original observation.
 f['unit']=f['unit']+'|row'+f['_row'].astype(str);finish(f,'KFM','mbon_fish')

def lter():
 # Original quadrat group measurements: group means are not individual lengths.
 q=read_source('lter_benthic_19_33');q['_row']=np.arange(1,len(q)+1);q=q[q.SP_CODE.isin(['SFL','SFS','SPL','SPS'])].copy();q['Taxa']=q.SP_CODE.map(TAXA);q['site']=q.SITE;q['date']=pd.to_datetime(q.DATE).dt.strftime('%Y-%m-%d');q['unit']=q.TRANSECT.astype(str)+'|'+q.QUAD.astype(str)+'|'+q.SIDE.astype(str)
 q['density']=q.COUNT.where(q.COUNT.ge(0))/q.AREA;q['biomass']=[c*mass(t,z*10)/a if c>=0 and z>0 else 0 if c==0 else np.nan for c,t,z,a in zip(q.COUNT,q.Taxa,q.SIZE,q.AREA)]
 q=q.rename(columns={'AREA':'area'});key=['site','date','unit','Taxa','area'];q=q.groupby(key,dropna=False).agg(density=('density',strict_sum),biomass=('biomass',strict_sum),source_rows=('_row',lambda s:';'.join(map(str,s)))).reset_index();q['size_basis']='same_quadrat_size_class_mean';q['bio_mean_size_proxy']=True;finish(q,'LTER','lter_benthic_19_33')
 for source_id in ['lter_kelp','lter_fish','lter_lobster']:
  x=read_source(source_id);x['_row']=np.arange(1,len(x)+1);x['site']=x.SITE;x['date']=pd.to_datetime(x.DATE).dt.strftime('%Y-%m-%d');x['area']=x.AREA
  if source_id=='lter_kelp':
   x=x[x.SP_CODE.eq('MAPY')].copy();x['Taxa']='M. pyrifera';x['unit']=x.TRANSECT.astype(str)+'|'+x.QUAD.astype(str)+'|'+x.SIDE.astype(str);x['density']=x.FRONDS.where(x.FRONDS.ge(0))/x.area;x['biomass']=x.density*x.MONTH.map(cal)*1000;x['size_basis']='same_record_stipes'
  else:
   if source_id=='lter_fish':
    x=x[x.SP_CODE.eq('SPUL')].copy();x['Taxa']='S. pulcher';x['unit']=x.TRANSECT.astype(str);length=x.SIZE;x['size_basis']='same_record_fish_length'
   else:
    x['Taxa']='P. interruptus';x['unit']=x.TRANSECT.astype(str)+'|'+x.REPLICATE.astype(str);length=x.SIZE_MM;x['size_basis']='same_record_lobster_length'
   x['imputed_density']=0.;x['imputed_pool']='';x['density']=x.COUNT.where(x.COUNT.ge(0))/x.area;x['individual_mass']=mass(x.Taxa.iloc[0],length.where(length.gt(0)));x['biomass']=x.density*x.individual_mass
   if source_id=='lter_lobster':
    for (site,date),g in x.groupby(['site','date']):
     donor=g[g.COUNT.gt(0)&g.individual_mass.notna()];missing=g[g.COUNT.gt(0)&g.individual_mass.isna()]
     if len(donor) and len(missing):
      mu=np.average(donor.individual_mass,weights=donor.COUNT);x.loc[missing.index,'biomass']=x.loc[missing.index,'density']*mu;x.loc[missing.index,'size_basis']='imputed_same_site_date_mean_mass'
      pool=sid('lter_lobster_imputation',site,date);x.loc[missing.index,'imputed_density']=x.loc[missing.index,'density'];x.loc[missing.index,'imputed_pool']=pool
      for _,r in donor.iterrows():decisions.append(dict(pool_id=pool,source_record_id=source_id+':'+SOURCES[source_id]['sha256']+':'+str(r['_row']),mass=r.individual_mass,frequency=r.COUNT,site=site,date=date))
   x.loc[x.COUNT.eq(0),'biomass']=0
  if 'imputed_density' not in x:x['imputed_density']=0.;x['imputed_pool']=''
  key=['site','date','unit','Taxa','area'];x=x.groupby(key,dropna=False).agg(density=('density',strict_sum),biomass=('biomass',strict_sum),source_rows=('_row',lambda s:';'.join(map(str,s))),size_basis=('size_basis',lambda s:';'.join(sorted(set(s)))),imputed_density=('imputed_density',strict_sum),imputed_pool=('imputed_pool',lambda s:next((v for v in s if v),''))).reset_index();finish(x,'LTER',source_id)

def satellite():
 x=read_source('satellite_canopy');bad=read_source('empty_reference_decisions');x=x.merge(bad[['MPA','status','rep']].assign(_exclude=True),on=['MPA','status','rep'],how='left');x=x[x._exclude.isna()];x['site']=x.MPA+'|'+x.status+'|'+x.rep.astype(str);x=x.melt(id_vars=['site'],value_vars=[c for c in x if c.isdigit()],var_name='year',value_name='biomass');x['date']=x.year+'-01-01';x['unit']='polygon';x['Taxa']='M. pyrifera';x['area']=1.;x['density']=np.nan;x['source_rows']='registered_polygon';x['size_basis']='provider_canopy_biomass';finish(x,'Landsat','satellite_canopy')

if __name__=='__main__':
 pisco();mbon();lter();satellite();pd.DataFrame(decisions).to_csv(PRIVATE/'lobster_imputation_members.csv',index=False);x=pd.concat(frames,ignore_index=True);x.to_parquet(PRIVATE/'measurement_units.parquet',index=False);x.groupby(['Source','Taxa']).agg(units=('unit_id','size'),sites=('site','nunique'),first_year=('year','min'),last_year=('year','max'),missing_counts=('density',lambda s:s.isna().sum()),missing_biomass=('biomass',lambda s:s.isna().sum())).to_csv(OUT/'measurement_summary.csv');write_inputs();print('Taxon-by-survey-unit records',len(x),'unique physical units',x.unit_id.nunique())
