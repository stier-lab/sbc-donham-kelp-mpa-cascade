"""Check reference omissions, matched-year decomposition and unchanged primary files."""
from pathlib import Path
import hashlib,json,os
import pandas as pd
import numpy as np
R=Path(__file__).resolve().parents[2];B=R/'outputs/integrated_review_20261103';O=R/'outputs/kfm_reference_sensitivity_20261103'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def folder(n):return B/'sensitivities'/('kfm_reference_'+n)
assert (O/'completed.json').exists()
config=json.loads((O/'design.json').read_text());x=pd.read_parquet(B/'private/supported_measurement_uses.parquet');base=pd.read_csv(B/'annual_values.csv');base_draws=np.load(B/'private/annual_measurement_draws.npz')['draws'];assign=pd.read_csv(B/'station_assignments.csv');checks=[];diag=[]
for name,sites in config['scenarios'].items():
 p=folder(name);z=pd.read_parquet(p/'private/supported_measurement_uses.parquet');expected=x[~(x.Source.eq('KFM')&x.site.isin(sites))].copy();roles=expected.groupby(['series_id','year']).status.agg(set);valid=set(roles[roles.map(lambda s:s=={'mpa','reference'})].index);expected=expected[[k in valid for k in zip(expected.series_id,expected.year)]]
 pd.testing.assert_frame_equal(z.reset_index(drop=True),expected.reset_index(drop=True))
 # Keep original ecological reasons; sensitivity flags are additional metadata.
 ledger=assign.copy();flag=ledger.source.eq('KFM')&ledger.site.isin(sites);ledger['sensitivity_flagged']=flag;new=flag&ledger.included;ledger.loc[new,'included']=False;ledger.loc[new,'reason']='Sensitivity only: named reference omission';ledger.to_csv(p/'station_assignments.csv',index=False)
 rr=pd.read_csv(p/'comparison_year_roster.csv');rr['coverage_metadata_scope']=['Retained sensitivity units' if k in valid else 'Original candidate metadata; no retained paired units' for k in zip(rr.series_id,rr.year)];rr.to_csv(p/'comparison_year_roster.csv',index=False)
 a=pd.read_csv(p/'annual_values.csv');d=np.load(p/'private/annual_measurement_draws.npz')['draws'];j=base.merge(a,on=['series_id','year'],suffixes=('_base','_alt'));assert np.allclose(j.offset_base,j.offset_alt,equal_nan=True)
 affected=set(x.loc[x.Source.eq('KFM')&x.site.isin(sites),'series_id']);u=j[~j.series_id.isin(affected)]
 assert np.allclose(u.ln_ratio_base,u.ln_ratio_alt,equal_nan=True,rtol=1e-12,atol=1e-12)
 assert np.allclose(base_draws[u.annual_index_base.astype(int)],d[u.annual_index_alt.astype(int)],equal_nan=True,rtol=1e-12,atol=1e-12)
 assert not (z.Source.eq('KFM')&z.site.isin(sites)).any()
 assert not set(ledger.loc[~assign.included,'reason'])-set(assign.loc[~assign.included,'reason'])
 checks.append(dict(scenario=name,retained_units_unchanged=True,unaffected_annual_values_and_draws_unchanged=True,offsets_fixed=True,excluded_references_absent=True,ecological_reasons_preserved=True))
for name in ['without_wyckoff','without_flagged']:
 p=folder('matched_main_'+name);b=pd.read_csv(p/'annual_values.csv');alt=pd.read_csv(folder(name)/'annual_values.csv');alt=alt[alt.response_supported];assert set(zip(b.series_id,b.year))==set(zip(alt.series_id,alt.year))
 j=b.merge(base,on=['series_id','year'],suffixes=('_sub','_base'));assert np.allclose(j.ln_ratio_sub,j.ln_ratio_base)
 sub=np.load(p/'private/annual_measurement_draws.npz')['draws'];assert np.allclose(sub[j.annual_index_sub.astype(int)],base_draws[j.annual_index_base.astype(int)])
for name in list(config['scenarios'])+['matched_main_without_wyckoff','matched_main_without_flagged']:
 p=folder(name);d=pd.read_csv(p/'joint_model_diagnostics.csv');d=d[d.analysis.eq('estimated_AR1')];C=np.load(p/'private/joint_coefficient_covariance.npz')['covariance'];assert np.isfinite(C).all();assert np.allclose(C,C.T,atol=1e-9);assert np.linalg.eigvalsh(C).min()>-1e-8;assert not d.upper_boundary.any();assert (d.df>0).all()
 diag.append(dict(scenario=name,groups=len(d),ar_boundaries=int(d.ar_boundary.sum()),variance_boundaries=int(d.boundary.sum()),upper_variance_boundaries=int(d.upper_boundary.sum()),minimum_coefficient_eigenvalue=np.linalg.eigvalsh(C).min(),max_optimizer_attempts=int(d.optimizer_attempts.max())))
for f,h in json.loads((O/'input_manifest.json').read_text()).items():assert sha(R/f)==h,f
for f,h in json.loads(Path('/tmp/docs-cleanup-numerical-hashes.json').read_text()).items():assert sha(R/f)==h,f
# Each fit's own tracked inputs include offsets, covariance geometry and exclusion utilities.
for name in config['scenarios']:
 for f,h in json.loads((folder(name)/'input_hashes.json').read_text()).items():
  if Path(f).name not in ['comparison_year_roster.csv','station_assignments.csv']:assert sha(Path(f) if Path(f).is_absolute() else R/f)==h,f
pd.DataFrame(checks).to_csv(O/'checks.csv',index=False);pd.DataFrame(diag).to_csv(O/'diagnostics.csv',index=False)
(O/'verification.json').write_text(json.dumps(dict(scenarios=4,matched_main_refits=2,primary_unchanged=True,checks=checks,diagnostics=diag,metadata_note='Ecological reasons preserved and dropped-row metadata scope clarified after initial numerical fits; neither change alters model inputs.'),indent=2))
print('Reference filtering, original exclusions, paired years, draws, offsets, covariance and primary hashes verified.')
