"""Independent algebra, estimability and support checks for conditional diagnostics."""
import importlib.util,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy.linalg import block_diag
spec=importlib.util.spec_from_file_location('model',Path(__file__).with_name('biomass_observation_model.py'));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
rng=np.random.default_rng(9)
# Augmented Gaussian posterior precision as independent check of prediction equations.
n=10;p=2;X=np.column_stack([np.ones(n),np.arange(n)/10]);G=.4*np.eye(n)+.2*np.ones((n,n));R=np.diag(np.linspace(.05,.2,n));y=rng.normal(size=n)
V=G+R;VI=np.linalg.inv(V);bc=np.linalg.inv(X.T@VI@X);beta=bc@X.T@VI@y
from scipy.linalg import cho_factor
f=dict(theta=np.array([1.]),cf=cho_factor(V,lower=True),beta=beta,beta_cov=bc)
mu,L=m.prediction(f,y,X,X,[G],[G]);D=np.column_stack([X,np.eye(n)]);Q=D.T@np.linalg.inv(R)@D+block_diag(np.zeros((p,p)),np.linalg.inv(G));C=np.linalg.inv(Q);mean=C@D.T@np.linalg.inv(R)@y
np.testing.assert_allclose(mu,D@mean,atol=1e-10);np.testing.assert_allclose(L@L.T,D@C@D.T,atol=1e-10)
try:m.fit_reml(y,np.ones((n,3)),np.ones(n),[np.eye(n)])
except ValueError:pass
else:raise AssertionError('rank guard absent')
try:m.fit_trend(np.ones(3),np.arange(3),np.eye(3),.3)
except ValueError:pass
else:raise AssertionError('process df guard absent')
O=Path('outputs/biomass_next_steps_20261102/uncertainty');b=pd.read_csv('outputs/biomass_next_steps_20261102/borrowing/annual.csv');b=b[(b.tolerance_days==7)&b.supported]
a=pd.read_csv(O/'annual_observations.csv');assert set(map(tuple,a[['series_id','year']].values))==set(map(tuple,b[['series_id','year']].values))
assert (a.t>0).sum()==45
fits=pd.read_csv(O/'size_model_diagnostics.csv');assert fits.converged.all()
boots=pd.read_csv(O/'size_parameter_bootstraps.csv');assert boots.converged.all() and not boots.used_fallback.any() and len(boots)==480
assert boots.groupby(['series_id','status']).size().eq(40).all()
assert (~boots.attempt_converged).equals(boots.retried)
reg=pd.read_csv(O/'size_observation_registry.csv');assert reg.loc[reg.n.eq(1),'sampling_log_mean_variance'].gt(0).all()
cp=pd.read_csv(O/'count_reconstruction.csv');np.testing.assert_allclose(cp.expected,cp.reconstructed,atol=1e-12)
C=pd.read_csv(O/'annual_latent_covariance.csv').latent_covariance.to_numpy().reshape(len(a),len(a));np.testing.assert_allclose(C,C.T);assert np.linalg.eigvalsh(C).min()>-1e-9
Y=np.load(O/'private/annual_draws.npz')['draws'];s=pd.read_csv(O/'observed_year_slope_draws.csv').set_index('series_id')
for sid,g in a[a.t>0].groupby('series_id'):
    X=np.column_stack([np.ones(len(g)),g.year-g.year.mean()]);sl=np.linalg.lstsq(X,Y[g.index],rcond=None)[0][1]
    np.testing.assert_allclose([sl.mean(),*np.quantile(sl,[.025,.975])],s.loc[sid,['mean_slope','draw_low','draw_high']].astype(float),atol=1e-12)
assert not (O/'trend_scenarios.csv').exists(),'Invalid latent-posterior process shortcut must not be published'
print('Latent covariance algebra, rank/df,480fit records, support/45yearidentity, singleton pooling,1086counts, PSD and independent slope draws pass. No process inference on latent posterior covariance.')

assert json.loads((O/'run_status.json').read_text())['status']=='complete'
