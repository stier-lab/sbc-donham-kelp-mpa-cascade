from pathlib import Path
import pandas as pd,numpy as np,json,ast,hashlib
from scipy.stats import t as student
O=Path('outputs/kfm_remaining_sources_20261102')
p=pd.read_csv('outputs/biomass_support_20261102/private/pool_members.csv',low_memory=False);inv=pd.read_csv('outputs/biomass_support_20261102/animal_pool_inventory.csv').set_index(['series_id','year','status']);off=pd.read_csv('outputs/all_program_station_refits_20261102/offsets.csv').set_index('series_id').offset
l=pd.read_csv('/Users/adrianstier/sbc-donham-kelp-mpa-data/data/source_archive/kfm_size_feasibility_2026-10-29/legacy/Invertebrate Size Frequency.csv');a=pd.read_csv('outputs/kfm_scott_reconciliation_20261102/review_annual_input.csv');v=a.copy();records=[];sizes=[]
for sid,site,y,code,power in [('Harris Point SMR|S. purpuratus|KFM|Bio','a-k-02',1987,11006,2.870),('Scorpion SMR|M. franciscanus|KFM|Bio','a-k-08',1995,11005,2.917)]:
 role=p.loc[p.series_id.eq(sid)&p.year.eq(y)&p.site.eq(site)&p.campus.eq('KFM'),'status'].unique();assert len(role)==1;role=role[0]
 pool=p[p.series_id.eq(sid)&p.year.eq(y)&p.status.eq(role)];target=pool[pool.site.eq(site)&pool.campus.eq('KFM')];oldmean=np.average(.00059*pool['size']**power,weights=pool['count']);z=l[l.Site.eq(int(site[-2:]))&l.Year.eq(y)&l.Species.eq(code)].copy();z['mm']=(z.Size*10).round(6);keep=pool[~(pool.site.eq(site)&pool.campus.eq('KFM'))];newmean=((.00059*keep['size']**power*keep['count']).sum()+(.00059*z.mm**power).sum())/(keep['count'].sum()+len(z));factor=newmean/oldmean
 oldhist=target.groupby('size')['count'].sum();newhist=z.groupby('mm').size()
 for mm in sorted(set(oldhist.index)|set(newhist.index)):
  if oldhist.get(mm,0)!=newhist.get(mm,0):sizes.append(dict(site=site,year=y,size_mm=mm,combined_count=int(oldhist.get(mm,0)),legacy_count=int(newhist.get(mm,0))))
 values={st:inv.loc[(sid,y,st),'value'] for st in ['mpa','reference']};values[role]*=factor;newratio=np.log((values['mpa']+off[sid])/(values['reference']+off[sid]));ix=v.series_id.eq(sid)&v.year.eq(y);assert ix.sum()==1;oldratio=v.loc[ix,'lnDiff'].iloc[0];v.loc[ix,'lnDiff']=newratio
 records.append(dict(series_id=sid,site=site,year=y,affected_role=role,combined_n=int(target['count'].sum()),legacy_n=len(z),old_pool_mean_mass=oldmean,new_pool_mean_mass=newmean,affected_mass_factor=factor,old_review_lnratio=oldratio,legacy_review_lnratio=newratio,delta_review_lnratio=newratio-oldratio,original_date=';'.join(z.Date.unique())))
c=pd.DataFrame(records);c.to_csv(O/'older_mass_comparison.csv',index=False);pd.DataFrame(sizes).to_csv(O/'exact_size_disagreements.csv',index=False);v.to_csv(O/'legacy_sensitivity_review_annual.csv',index=False)
# Same step contrast and conditional lm uncertainty as the frozen primary.
frozen=pd.read_csv('data/harmonized/harmonized_response_ratios.csv');raw=pd.read_csv('data/harmonized/harmonized_raw_responses.csv');f=frozen[frozen.source.eq('KFM')&frozen.CA_MPA_Name_Short.eq('Scorpion SMR')&frozen.y.eq('Mesocentrotus franciscanus')&frozen.resp.eq('Bio')].copy();r=raw[raw.source.eq('KFM')&raw.CA_MPA_Name_Short.eq('Scorpion SMR')&raw.taxon_name.eq('Mesocentrotus franciscanus')&raw.resp.eq('Bio')&raw.year.eq(1995)];base=f.loc[f.year.eq(1995)].iloc[0];offset_raw=off['Scorpion SMR|M. franciscanus|KFM|Bio'];scales=[base[st]/(r.loc[r.status.eq(st),'value'].iloc[0]+offset_raw) for st in ['mpa','reference']];assert np.isclose(scales[0],scales[1],atol=1e-12);offset=offset_raw*scales[0]
# Verify scaling and offset on every primary role/year, not just the modified year.
rr=raw[raw.source.eq('KFM')&raw.CA_MPA_Name_Short.eq('Scorpion SMR')&raw.taxon_name.eq('Mesocentrotus franciscanus')&raw.resp.eq('Bio')].pivot(index='year',columns='status',values='value')
for st in ['mpa','reference']:assert np.allclose(f[st].values,(rr.loc[f.year,st].values+offset_raw)*scales[0],atol=1e-12)

factor=c.loc[c.site.eq('a-k-08'),'affected_mass_factor'].iloc[0];f['legacy']=f.lnDiff;f.loc[f.year.eq(1995),'legacy']=np.log(base.mpa/((base.reference-offset)*factor+offset))
e=[]
for label,col in [('current_primary','lnDiff'),('legacy_source_sensitivity','legacy')]:
 X=np.column_stack([np.ones(len(f)),f.BA.eq('After').astype(int)]);beta=np.linalg.lstsq(X,f[col],rcond=None)[0];res=f[col]-X@beta;df=len(f)-2;se=np.sqrt((res@res/df)*np.linalg.inv(X.T@X)[1,1]);ci=student.ppf(.975,df)*se;e.append(dict(scenario=label,estimate=beta[1],se=se,low=beta[1]-ci,high=beta[1]+ci,n=len(f),before_n=f.BA.eq('Before').sum(),after_n=f.BA.eq('After').sum(),offset=offset))
e=pd.DataFrame(e);s=pd.read_csv('data/sumstats_final.csv');s=s[s.MPA.eq('Scorpion SMR')&s.Source.eq('KFM')&s.Taxa.eq('M. franciscanus')&s.Resp.eq('Bio')].iloc[0];assert abs(e.iloc[0].estimate-s.Mean)<1e-10 and abs(e.iloc[0].se-s.SE)<1e-10;e.to_csv(O/'scorpion_primary_step_sensitivity.csv',index=False);f.to_csv(O/'scorpion_primary_annual_sensitivity.csv',index=False)
# Post-protection slopes: exactly unchanged since the changed years precede protection.
post=a.t.gt(0);pd.testing.assert_frame_equal(a.loc[post],v.loc[post]);tree=ast.parse(Path('code/tools/audit_station_coverage_refits.py').read_text());exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['fit','slope']],type_ignores=[]),'fit','exec'));rows=[]
for sid in c.series_id:
 for label,x in [('current_review',a),('legacy_source_sensitivity',v)]:rows.append(dict(series_id=sid,scenario=label,**slope(x[x.series_id.eq(sid)&x.t.gt(0)],'lnDiff')))
pd.DataFrame(rows).to_csv(O/'post_protection_slopes.csv',index=False)
keys=pd.read_csv('data/harmonized/reference_primary_keys.csv');uses=[]
for site in ['a-k-09','a-k-15']:
 for tax in ['Strongylocentrotus franciscanus','Strongylocentrotus purpuratus']:
  z=p[p.campus.eq('KFM')&p.site.eq(site)&p.year.eq(2023)&p.classcode.eq(tax)];mpa='Scorpion SMR' if site=='a-k-09' else 'Santa Barbara Island SMR';short='M. franciscanus' if 'franciscanus' in tax else 'S. purpuratus';kk=keys[keys.MPA.eq(mpa)&keys.Taxa.eq(short)&keys.Resp.eq('Bio')];uses.append(dict(site=site,year=2023,taxon=tax,review_uses=len(z),primary_histories_using_2023=int(kk.observed_years.fillna('').str.split(',').apply(lambda y:'2023' in y).sum()),status='Outside current matched biomass histories; source identity still unverified'))
pd.DataFrame(uses).to_csv(O/'2023_analysis_use.csv',index=False)
h=json.loads(Path('outputs/researcher_scope_20261101/protected_hashes.json').read_text());assert all(hashlib.sha256(Path(k).read_bytes()).hexdigest()==vv for k,vv in h.items());print(c.to_string(index=False));print(e.to_string(index=False));print('All observed post-protection rows unchanged; 487 protected files unchanged.')
