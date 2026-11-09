"""Integration invariants, checked against source keys and frozen decisions."""
from common import *
u=pd.read_parquet(PRIVATE/'measurement_units.parquet');x=pd.read_parquet(PRIVATE/'supported_measurement_uses.parquet');a=read_tracked(OUT/'annual_values.csv');ro=read_tracked(OUT/'comparison_year_roster.csv');assign=read_tracked(OUT/'station_assignments.csv');vis=read_tracked(OUT/'visit_biomass_support.csv');p=read_tracked(OUT/'size_pool_inventory.csv');checks=[]
def check(name,ok,detail=''):
 checks.append(dict(check=name,passed=bool(ok),detail=detail))
 if not ok:raise AssertionError(name+': '+detail)
check('unique source units within species',not u.duplicated(['unit_id','Taxa']).any())
check('positive areas',u.area.gt(0).all())
check('no negative density',u.density.dropna().ge(0).all())
check('no negative biomass',x.biomass.dropna().ge(0).all())
check('excluded physical sites absent from admitted uses',not any(blocked(r.Source,r.site,r.Taxa) for r in x[['Source','site','Taxa']].drop_duplicates().itertuples()))
check('benthic exclusions not silently imposed on fish',set(['a-k-21','a-k-27','a-k-28']).issubset(set(x.loc[x.Source.eq('KFM')&x.Taxa.eq('S. pulcher'),'site'])))
check('one annual row per comparison/year',not a.duplicated(['series_id','year']).any())
check('actual calendar clock',(a.t==a.year-a.MPA_Start).all())
check('both roles in every admitted comparison/year',x.groupby(['series_id','year']).status.agg(set).map(lambda s:s=={'mpa','reference'}).all())
check('no silently discarded unknowns',all(not np.isfinite(r.ln_ratio) for r in a[(a.protected.isna()|a.reference.isna())].itertuples()))
check('every admitted row has explicit assignment',len(x.merge(assign[assign.included][['series_id','site','role']].rename(columns={'role':'status'}),on=['series_id','site','status'],how='inner',validate='many_to_one'))==len(x))
check('shared zones within comparison/year',x.groupby(['series_id','year']).apply(lambda z:z.groupby('site').zone.agg(lambda v:';'.join(sorted(set(v)))).nunique()==1,include_groups=False).all())
check('unknown density survives zero handling',u.density.isna().any())
# Source hashes and manuscript hold cover data separately from Git state.
hashes=json.loads((OUT/'input_hashes.json').read_text());changed=[k for k,v in hashes.items() if sha(k)!=v];check('inputs match recorded hashes',not changed,str(changed))
protected=json.loads((ROOT/'outputs/researcher_scope_20261101/protected_hashes.json').read_text());changed=[k for k,v in protected.items() if sha(k)!=v];check('canonical and manuscript hold',not changed,str(len(protected))+' protected files')
j=read_tracked(OUT/'joint_model_diagnostics.csv');check('joint coefficient covariance nonnegative',j.covariance_min_eigenvalue.ge(-1e-8).all());check('positive residual degrees of freedom',j.df.gt(0).all())
# Density coverage never takes its years from biomass availability.
d=a[a.Resp.eq('Den')];b=a[a.Resp.eq('Bio')].copy();b['series_id']=b.series_id.str.replace('|Bio','|Den',regex=False);z=d.merge(b[['series_id','year','response_supported']],on=['series_id','year'],suffixes=('_density','_biomass'));check('unsupported biomass retains supported density',((z.response_supported_density)&(~z.response_supported_biomass)).any())
pd.DataFrame(checks).to_csv(OUT/'validation_checks.csv',index=False);(OUT/'verification.json').write_text(json.dumps(dict(checks=len(checks),passed=True,protected_files=len(protected),annual_rows=len(a),supported_annual_rows=int(a.response_supported.sum()),candidate_histories=ro.series_id.nunique(),fitted_histories=int(pd.read_csv(OUT/'model_results.csv').series_id.nunique()),uncertainty='Conditional; see limitations ledger'),indent=2));print('PASS',len(checks),'invariants')

# The exact cross-year measurement-zero rule must remain valid as sources change.
u=pd.read_parquet(PRIVATE/'supported_measurement_uses.parquet')
assert u[u.Source.isin(['KFM','PISCO'])].groupby('unit_id').year.nunique().max()==1
members=pd.read_parquet(PRIVATE/'size_pool_members.parquet');pools=pd.read_csv(OUT/'size_pool_inventory.csv')
z=members.merge(pools[['pool_id','Source','year']],on='pool_id')
assert z[z.Source.eq('KFM')].groupby('source_record_id').year.nunique().max()==1
from common import measurement_covariance
rng=np.random.default_rng(9);d=rng.normal(size=(4,100));years=np.array([2000,2000,2001,2001])
c=measurement_covariance(d,years,'KFM','Den','P. interruptus')
assert np.all(c[:2,2:]==0) and np.linalg.eigvalsh(c).min()>0
assert c[0,1]==np.cov(d)[0,1]
print('Independent-year measurement covariance and retained shared-year covariance verified.')
