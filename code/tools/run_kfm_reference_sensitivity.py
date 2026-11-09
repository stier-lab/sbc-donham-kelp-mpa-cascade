"""KFM reference sensitivity: same pipeline, explicit omissions, no primary promotion."""
from pathlib import Path
import os,sys,json,subprocess,hashlib,concurrent.futures
import pandas as pd
import numpy as np
R=Path(__file__).resolve().parents[2];B=R/'outputs/integrated_review_20261103';O=R/'outputs/kfm_reference_sensitivity_20261103';O.mkdir(exist_ok=True)
SCENARIOS={'without_webster':['a-k-35'],'without_miracle':['a-k-21'],'without_wyckoff':['a-k-01'],'without_flagged':['a-k-01','a-k-21','a-k-35','a-k-05']}
x=pd.read_parquet(B/'private/supported_measurement_uses.parquet');r=pd.read_csv(B/'comparison_year_roster.csv')
inputs=[Path(__file__),R/'code/integrated_review/linear_mixed.py',R/'outputs/all_program_station_refits_20261102/offsets.csv',B/'station_annual_values.csv',R.parent/'sbc-donham-kelp-mpa-data/code/R/01_utils.R']+list(B.glob('*.csv'))+[B/'completed.json']+[B/'private/supported_measurement_uses.parquet',B/'private/size_pool_members.parquet',B/'comparison_year_roster.csv',B/'annual_values.csv',B/'private/annual_measurement_draws.npz',B/'station_assignments.csv',B/'visit_biomass_support.csv',B/'size_pool_inventory.csv',R/'config/integrated_review_policy.json']+list((R/'code/integrated_review').glob('*.py'))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
completion=json.loads((B/'completed.json').read_text())
assert sha(R/'config/integrated_review_policy.json')==completion['policy_sha256'], 'Review policy changed; rebuild main review first'
assert all(sha(R/f)==h for f,h in completion['source_code_hashes'].items()), 'Main review code changed; rebuild before sensitivity'
original={os.path.relpath(p,R):sha(p) for p in inputs};(O/'input_manifest.json').write_text(json.dumps(original,indent=2))
settings={'scenarios':SCENARIOS,'primary_changed':False,'rules':'Original supported units, zones, offsets, draw keys and aggregation; only flagged KFM reference rows removed. No new reference added. Role-incomplete years omitted. Matched-main counterparts hold supported series and years fixed; diagnostics, not adopted exclusions.'}
(O/'design.json').write_text(json.dumps(settings,indent=2))
def folder(name):return B/'sensitivities'/('kfm_reference_'+name)
def link(src,dst):
 dst.parent.mkdir(parents=True,exist_ok=True)
 if not dst.exists():dst.symlink_to(src)
def launch(name,script):
 target=folder(name);env=os.environ.copy();env.update(INTEGRATED_REVIEW_VARIANT='kfm_reference_'+name,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
 with (target/(script+'.log')).open('w') as f:subprocess.run([sys.executable,str(R/'code/integrated_review'/f'{script}.py')],cwd=R,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
def prepare(name,sites):
 target=folder(name);(target/'private').mkdir(parents=True,exist_ok=True)
 omitted=x.Source.eq('KFM')&x.site.isin(sites);assert x.loc[omitted,'status'].eq('reference').all()
 z=x[~omitted].copy();roles=z.groupby(['series_id','year']).status.agg(set);valid=set(roles[roles.map(lambda s:s=={'mpa','reference'})].index)
 z=z[[k in valid for k in zip(z.series_id,z.year)]];z.to_parquet(target/'private/supported_measurement_uses.parquet',index=False)
 roster=r.copy()
 for (sid,yr),g in z.groupby(['series_id','year']):
  mask=roster.series_id.eq(sid)&roster.year.eq(yr)
  for col,vals in [('stations',g.site),('protected_stations',g.loc[g.status.eq('mpa'),'site']),('reference_stations',g.loc[g.status.eq('reference'),'site']),('shared_zones',g.zone)]:roster.loc[mask,col]=';'.join(sorted(vals.unique()))
  roster.loc[mask,'visits']=len(g[['site','date','zone']].drop_duplicates());roster.loc[mask,'first_survey_date']=g.date.min();roster.loc[mask,'last_survey_date']=g.date.max()
 missing=roster.survey_eligible & ~pd.Series([k in valid for k in zip(roster.series_id,roster.year)],index=roster.index)
 roster.loc[missing,'survey_eligible']=False;roster.loc[missing,'reason']='No protected/reference pair after named reference omission'
 roster['coverage_metadata_scope']=['Retained sensitivity units' if k in valid else 'Original candidate metadata; no retained paired units' for k in zip(roster.series_id,roster.year)]
 roster.to_csv(target/'comparison_year_roster.csv',index=False)
 assignments=pd.read_csv(B/'station_assignments.csv');bad=assignments.source.eq('KFM')&assignments.site.isin(sites);assignments['sensitivity_flagged']=bad;newly=bad & assignments.included;assignments.loc[newly,'included']=False;assignments.loc[newly,'reason']='Sensitivity only: named reference omission';assignments.to_csv(target/'station_assignments.csv',index=False)
 for f in ['visit_biomass_support.csv','size_pool_inventory.csv']:link(B/f,target/f)
 link(B/'private/size_pool_members.parquet',target/'private/size_pool_members.parquet')
 pd.DataFrame([{'scenario':name,'site':site,'omitted_unit_uses':int((omitted&x.site.eq(site)).sum())} for site in sites]).to_csv(target/'omissions.csv',index=False)
 print('Prepared',name,len(z),'unit uses',flush=True)
 launch(name,'estimate');print('Estimated',name,flush=True);launch(name,'joint_models');print('Fitted',name,flush=True)
def matched_main(name):
 # Copy the original annual point values/draws on the alternative's supported keys.
 target=folder('matched_main_'+name);(target/'private').mkdir(parents=True,exist_ok=True)
 alt=pd.read_csv(folder(name)/'annual_values.csv');alt=alt[alt.response_supported]
 base=pd.read_csv(B/'annual_values.csv');base=base[base.response_supported]
 keys=set(zip(alt.series_id,alt.year));g=base[[k in keys for k in zip(base.series_id,base.year)]].copy()
 # These cases must not gain previously unsupported years: otherwise stage a matched alternative too.
 assert keys==set(zip(g.series_id,g.year)), 'New supported years require paired intersection analysis'
 d=np.load(B/'private/annual_measurement_draws.npz')['draws'][g.annual_index.to_numpy(int)]
 g['annual_index']=np.arange(len(g));g.to_csv(target/'annual_values.csv',index=False);np.savez_compressed(target/'private/annual_measurement_draws.npz',draws=d)
 st=pd.read_csv(B/'station_annual_values.csv');st=st[[k in keys for k in zip(st.series_id,st.year)]];st.to_csv(target/'station_annual_values.csv',index=False)
 elig=[]
 for sid,h in g.groupby('series_id'):
  pre=int(h.t.lt(0).sum());post=int(h.t.ge(0).sum());elig.append(dict(series_id=sid,before_years=pre,after_years=post,eligible_model='before_after_and_post_trend' if pre>=3 and post>=3 else 'post_trend_only' if post>=3 else 'insufficient_post_years',eligible_four_before=pre>=4,first_year=h.year.min(),last_year=h.year.max(),years=';'.join(map(str,sorted(h.year))),missing_response_years=''))
 pd.DataFrame(elig).to_csv(target/'model_eligibility.csv',index=False);launch('matched_main_'+name,'joint_models');print('Fitted matched main',name,flush=True)
if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
  jobs=[ex.submit(prepare,n,s) for n,s in SCENARIOS.items()]
  for job in jobs:job.result()
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
  jobs=[ex.submit(matched_main,n) for n in ['without_wyckoff','without_flagged']]
  for job in jobs:job.result()
 assert all(sha(R/p)==h for p,h in original.items()),'Primary input/code changed during sensitivity'
 (O/'completed.json').write_text(json.dumps({'scenarios':list(SCENARIOS),'matched_main':['without_wyckoff','without_flagged'],'primary_unchanged':True},indent=2))
 print('All reference scenarios completed',flush=True)
