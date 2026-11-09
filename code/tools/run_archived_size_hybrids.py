"""Two matched-cohort joint refits separating size-pool point and sampling covariance."""
from pathlib import Path
import shutil,subprocess,os
import pandas as pd,numpy as np
R=Path(__file__).resolve().parents[2];M=R/'outputs/integrated_review_20261103';C=M/'sensitivities/emily_matched_current';A=M/'sensitivities/emily_matched_archived'
c=pd.read_csv(C/'annual_values.csv');a=pd.read_csv(A/'annual_values.csv');assert len(c)==len(a) and (c[['series_id','year']]==a[['series_id','year']]).all().all();assert c.response_supported.equals(a.response_supported)
cd=np.load(C/'private/annual_measurement_draws.npz')['draws'];ad=np.load(A/'private/annual_measurement_draws.npz')['draws'];assert cd.shape==ad.shape
mask=c.Source.eq('PISCO')&c.Resp.eq('Bio')&c.Taxa.isin(['M. franciscanus','S. purpuratus'])&c.response_supported
ix=np.flatnonzero(mask.to_numpy());assert len(ix)==30
for name,point,draws,donorpoint in [('emily_archived_point_current_cov',A,cd,c),('emily_current_point_archived_cov',C,ad,a)]:
 V=M/'sensitivities'/name;P=V/'private';P.mkdir(parents=True,exist_ok=True)
 for n in ['annual_values.csv','station_annual_values.csv','model_eligibility.csv','input_hashes.json']:shutil.copy2(point/n,V/n)
 target=pd.read_csv(point/'annual_values.csv');base=draws.copy();delta=target.loc[ix,'ln_ratio'].to_numpy()-donorpoint.loc[ix,'ln_ratio'].to_numpy();base[ix]+=delta[:,None]
 np.savez_compressed(P/'annual_measurement_draws.npz',draws=base)
 env=dict(os.environ,INTEGRATED_REVIEW_VARIANT=name)
 subprocess.run(['python3',str(R/'code/integrated_review/joint_models.py')],cwd=R,env=env,check=True)
 print(name,'complete')
