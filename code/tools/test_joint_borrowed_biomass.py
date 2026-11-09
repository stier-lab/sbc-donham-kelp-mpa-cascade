#!/usr/bin/env python3
"""Read-only independent checks. Run: python3 -B test_joint.py [--repo PATH].
No model runner is called. Guard tests fail if protected invalid inputs are accepted.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse,sys,importlib.util,hashlib,unittest,json
from pathlib import Path
sys.dont_write_bytecode=True
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,default=Path('/Users/adrianstier/sbc-donham-kelp-mpa-cascade'));a,rest=p.parse_known_args();ROOT=a.repo.resolve()
sys.path.insert(0,str(ROOT/'code/tools'))
import numpy as np
import pandas as pd
import joint_borrowed_biomass as m
from scipy.optimize import OptimizeResult
from unittest.mock import patch

def fixture():
    rows=[];counts=[]
    for role in ['mpa','reference']:
      for site in [role+'A',role+'B']:
       for year in [2013,2014,2016,2019,2021]:
        for zone in ['INNER','MID']:
         date=f'{year}-06-10';idx=len(rows);den=1+(idx%7)/3
         rows.append(dict(status=role,site=site,year=year,zone=zone,date=date,density=den,mass=10+idx%9,mean_mass=10+idx%9,relative_ss=1.5+.1*(idx%4),df=4.,n=5.,method='measured_same_visit',sample_id_physical=str(idx)))
         for tr,fac in enumerate([.4,1.,1.6]): counts.append(dict(site=site,year=year,zone=zone,date=date,method='swath',density=den*fac,physical_unit_id=f'{idx}:{tr}'))
    s=pd.DataFrame(rows);v=s.copy();return s,v,np.array([2013,2014,2016,2019,2021]),.7,pd.DataFrame(counts)

def build(fix=None,seed=42,rho=.3):
    return m.build(*(fixture() if fix is None else fix),np.random.default_rng(seed),rho)

def independent(y,X,C):
    # Whitening plus least squares: no production GLS or inverse formula reused.
    L=np.linalg.cholesky(C);wX=np.linalg.solve(L,X);wy=np.linalg.solve(L,y)
    b=np.linalg.lstsq(wX,wy,rcond=None)[0];Q,R=np.linalg.qr(wX,mode='reduced');ri=np.linalg.solve(R,np.eye(R.shape[0]));bc=ri@ri.T
    res=wy-wX@b;obj=2*np.log(np.diag(L)).sum()+2*np.log(abs(np.diag(R))).sum()+res@res
    return b,bc,obj

class JointTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
    cls.source_hash=hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest()
    cls.fix=fixture();cls.d=build(cls.fix);cls.f=m.fit(cls.d['y'],cls.d['X'],cls.d['S'],cls.d['K'])
 @classmethod
 def tearDownClass(cls):
    if hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest()!=cls.source_hash:raise AssertionError('Implementation changed during tests; rerun')
 def completed_run(self):
    main=ROOT/'outputs/joint_borrowed_biomass_20261102';marker=main/'completed.json'
    if not marker.exists():self.skipTest('Parent output is in progress; completed.json absent')
    stamp=marker.read_bytes();paths=list(main.glob('*.csv'))+list((main/'fits').glob('*.npz'))+[main/'input_hashes.json',main/'policy.json']
    state={p:hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    def check_stable():
      if not marker.exists() or marker.read_bytes()!=stamp or any(hashlib.sha256(p.read_bytes()).hexdigest()!=h for p,h in state.items()):raise AssertionError('Parent outputs changed during review; rerun on completed outputs')
    self.addCleanup(check_stable)
    return main,json.loads((main/'policy.json').read_text())
 def test_covariance_signed_incidence(self):
    s,v,years,off,counts=self.fix;d=self.d
    # Independently create one column per physical random-effect level across s/v.
    for j,label in enumerate(d['labels'][:6]):
      role,kind=label.split('_');cols={'site':['site'],'siteyear':['site','year'],'event':['site','year','zone','date']}[kind]
      levels=sorted(set(tuple(row) for row in s.loc[s.status.eq(role),cols].itertuples(index=False,name=None)))
      def incidence(frame):
       return np.array([[float(r.status==role and tuple(getattr(r,c) for c in cols)==lev) for lev in levels] for r in frame.itertuples()])
      B=np.vstack([incidence(s),-d['W']@incidence(v)])
      np.testing.assert_allclose(d['K'][j],B@B.T,atol=1e-12)
      self.assertGreaterEqual(np.linalg.eigvalsh(d['K'][j]).min(),-1e-10)
    np.testing.assert_allclose(d['K'][-1][len(s):,len(s):],.3**abs(years[:,None]-years[None,:]))
 def test_jacobian_unequal_aggregation_zero_offset(self):
    _,v,years,off,_=self.fix;v=v.drop(index=[2,7,14]).reset_index(drop=True)
    extra=v.iloc[[0]].copy();extra['date']='2013-07-10';v=pd.concat([v,extra],ignore_index=True);v.loc[0,'density']=0
    z=np.log(v.mass.to_numpy());h,W,A=m.aggregation(v,years,np.exp(z),off)
    def exact(z):
      b=v.assign(b=v.density.to_numpy()*np.exp(z));rolemeans=b.groupby(['year','status','site','zone']).b.mean().groupby(['year','status','site']).mean().groupby(['year','status']).mean().unstack()
      return np.log((rolemeans.loc[years,'mpa']+off)/(rolemeans.loc[years,'reference']+off)).to_numpy()
    np.testing.assert_allclose(h,exact(z),atol=1e-12)
    eps=1e-5;fd=np.column_stack([(exact(z+np.eye(len(z))[i]*eps)-exact(z-np.eye(len(z))[i]*eps))/(2*eps) for i in range(len(z))])
    np.testing.assert_allclose(W,fd,atol=1e-9);np.testing.assert_array_equal(W[:,0],0)
    np.testing.assert_allclose((h-W@z)+W@z,h,atol=1e-12)
    _,w0,_=m.aggregation(v,years,np.exp(z),0);np.testing.assert_allclose(w0.sum(axis=1),0,atol=1e-12)
 def test_gls_inverse_objective_gradient(self):
    d,f=self.d,self.f;b,bc,obj=independent(d['y'],d['X'],f['C'])
    np.testing.assert_allclose(f['beta'],b,rtol=1e-8,atol=1e-9);np.testing.assert_allclose(f['bc'],bc,rtol=1e-8,atol=1e-9)
    np.testing.assert_allclose(f['inv']@f['C'],np.eye(len(d['y'])),atol=1e-9);self.assertAlmostEqual(f['obj'],obj,places=7)
    def loss(q):return independent(d['y'],d['X'],d['S']+sum(x*k for x,k in zip(q,d['K']))+np.eye(len(d['y']))*1e-9)[2]
    eps=1e-6;fd=np.array([(loss(f['q']+np.eye(7)[j]*eps)-loss(f['q']-np.eye(7)[j]*eps))/(2*eps) for j in range(7)])
    np.testing.assert_allclose(f['gradient'],fd,rtol=1e-5,atol=1e-5)
 def test_reused_target_not_independent_observation(self):
    # Re-request an existing latent coordinate without adding an observation.
    # build correctly forbids duplicated visit observations; prediction reuse is separate.
    d=dict(self.d);n=len(d['Xt']);take=np.r_[np.arange(n),0]
    d['Xt']=d['Xt'][take];d['Ktt']=[k[np.ix_(take,take)] for k in d['Ktt']];d['Kcross']=[k[take] for k in d['Kcross']];d['W']=np.c_[d['W'],np.zeros(len(d['h']))]
    self.assertEqual(d['ns'],len(self.fix[0]));self.assertEqual(len(d['y']),len(self.fix[0])+len(self.fix[2]))
    mu,C=m.latent(d,self.f)
    self.assertAlmostEqual(mu[0],mu[-1],places=10);self.assertAlmostEqual(C[0,0]+C[-1,-1]-2*C[0,-1],0,places=10)
 def test_count_reuse_exact_and_fixed_expansion(self):
    # Duplicate annual aggregate shares all physical resampling weights exactly.
    s,v,yr,off,c=self.fix;h,W,A=m.aggregation(v,yr,v.mass.to_numpy(),off);A=np.concatenate([A,A[:1]],axis=0)
    S,rep=m.count_cov(v,A,v.mass.to_numpy(),off,c,np.random.default_rng(5),draws=100)
    np.testing.assert_array_equal(rep[0],rep[-1]);np.testing.assert_allclose(S[0],S[-1],atol=1e-14)
    np.testing.assert_allclose(S,np.cov(rep-np.arange(len(rep))[:,None]),atol=1e-12)
 def test_conditional_latent_against_augmented_normal_equations(self):
    d,f=self.d,self.f;G=sum(q*k for q,k in zip(f['q'][:6],d['Ktt']));H=np.c_[sum(q*k for q,k in zip(f['q'][:6],d['Kcross'])),-G@d['W'].T];T=np.c_[d['Xt'],np.zeros((len(d['Xt']),2))]
    # Universal kriging saddle system is independent of production conditional formula.
    M=np.block([[f['C'],d['X']],[d['X'].T,np.zeros((8,8))]])
    rhs=np.vstack([H.T,T.T]);sol=np.linalg.solve(M,rhs)
    mu=sol[:len(d['y'])].T@d['y'];C=G-rhs.T@sol
    gotmu,gotC=m.latent(d,f);np.testing.assert_allclose(gotmu,mu,atol=1e-8);np.testing.assert_allclose(gotC,C,atol=1e-8)
 def test_rank_df_and_exhausted_optimizer_fail_closed(self):
    d=self.d
    with self.assertRaises(ValueError):m.fit(d['y'],np.c_[d['X'],d['X'][:,0]],d['S'],d['K'])
    with self.assertRaises(ValueError):m.fit(np.ones(2),np.eye(2),np.eye(2),[np.eye(2)])
    fail=OptimizeResult(success=False,fun=1.,message='forced failure',x=np.ones(7),nit=1)
    with patch.object(m,'minimize',return_value=fail) as mocked:
      with self.assertRaises(RuntimeError):m.fit(d['y'],d['X'],d['S'],d['K'])
      self.assertEqual(mocked.call_count,3)
 def test_actual_slsqp_retry_and_attempt_logging(self):
    original=m.minimize;methods=[];d=self.d
    def wrapped(fun,x,*args,**kwargs):
      methods.append(kwargs['method'])
      if len(methods)<3:return OptimizeResult(success=False,fun=1.,message='forced LBFGS failure',x=self.f['q'].copy(),nit=1)
      return original(fun,x,*args,**kwargs)
    with patch.object(m,'minimize',side_effect=wrapped):f=m.fit(d['y'],d['X'],d['S'],d['K'])
    self.assertEqual(methods,['L-BFGS-B','L-BFGS-B','SLSQP']);self.assertEqual(f['optimizer_attempts'],3);self.assertTrue(f['success'])
    b,bc,_=independent(d['y'],d['X'],f['C']);np.testing.assert_allclose(f['beta'],b,atol=1e-8)
 def test_new_support_guards(self):
    s,v,yr,off,c=self.fix
    with self.assertRaises(ValueError):build((s,pd.concat([v,v.iloc[[0]]],ignore_index=True),yr,off,c))
    bad=v.copy();bad.loc[0,'method']='unsupported'
    with self.assertRaises(ValueError):build((s,bad,yr,off,c))
    for rho in [-.1,1.,float('nan')]:
      with self.assertRaises(ValueError):build(rho=rho)
 def test_smoke_gaussian_solve(self):
    files=list(Path('/tmp/joint_biomass_smoke/fits').glob('*.npz'))
    if not files:self.skipTest('Optional smoke fixture unavailable')
    for path in files:
      with np.load(path) as z:
        C=z['S']+np.einsum('k,kij->ij',z['q'],z['K'])+np.eye(len(z['y']))*1e-9;b,bc,_=independent(z['y'],z['X'],C)
        np.testing.assert_allclose(b,z['beta'],atol=1e-8);np.testing.assert_allclose(bc,z['beta_cov'],atol=1e-8)
 def test_single_site_component_disabled(self):
    s,v,yr,off,c=self.fix;s=s[~s.site.eq('mpaB')].copy();v=v[~v.site.eq('mpaB')].copy();d=build((s,v,yr,off,c))
    np.testing.assert_array_equal(d['K'][0],0);self.assertEqual(np.linalg.matrix_rank(d['X']),8)
 def test_invalid_positive_mass_and_missing_counts_rejected(self):
    s,v,yr,off,c=self.fix;v=v.copy();v.loc[0,'mass']=np.nan
    with self.assertRaises(AssertionError):build((s,v,yr,off,c))
    with self.assertRaises(AssertionError):build((s,self.fix[1],yr,off,c.iloc[3:]))
 def test_actual_support_and_saved_gaussian_fits(self):
    main,policy=self.completed_run()
    base=ROOT/'outputs/biomass_next_steps_20261102';s=pd.read_csv(base/'uncertainty/size_observation_registry.csv');v=pd.read_csv(base/'borrowing/visit_support.csv');v=v[v.tolerance_days.eq(7)]
    self.assertFalse(s.sample_id_physical.duplicated().any());self.assertFalse(v.duplicated(['series_id','site','status','year','zone','date']).any())
    ann=pd.read_csv(base/'borrowing/annual.csv');ann=ann[ann.tolerance_days.eq(7)&ann.t.gt(0)&ann.supported]
    admitted=v.merge(ann[['series_id','year']],on=['series_id','year'],validate='many_to_one')
    self.assertFalse(admitted.method.eq('unsupported').any());self.assertTrue(np.isfinite(admitted.mass).all())
    main=ROOT/'outputs/joint_borrowed_biomass_20261102';files=sorted((main/'fits').glob('*.npz'));self.assertEqual(len(files),18)
    for path in files:
      with self.subTest(file=path.name),np.load(path) as z:
        C=z['S']+np.einsum('k,kij->ij',z['q'],z['K'])+np.eye(len(z['y']))*1e-9;b,bc,obj=independent(z['y'],z['X'],C)
        np.testing.assert_allclose(z['beta'],b,rtol=1e-7,atol=1e-8);np.testing.assert_allclose(z['beta_cov'],bc,rtol=1e-7,atol=1e-8)
        self.assertEqual(np.linalg.matrix_rank(z['X']),8);self.assertTrue(np.isfinite(C).all())
    bs=pd.read_csv(main/'bootstrap_fits.csv');tr=pd.read_csv(main/'joint_trends.csv');self.assertEqual(len(bs),18*policy['bootstraps']);self.assertTrue(bs.converged.all())
    for row in tr.itertuples():
      g=bs[bs.series_id.eq(row.series_id)&bs.rho.eq(row.rho)];self.assertEqual(len(g),policy['bootstraps'])
      piv=(g.slope-row.estimate)/g.se;np.testing.assert_allclose([row.low,row.high],[row.estimate-np.quantile(piv,.975)*row.se,row.estimate-np.quantile(piv,.025)*row.se],atol=1e-10)
 def test_actual_rebuild_support_jacobians_and_smoke(self):
    main,policy=self.completed_run()
    base=ROOT/'outputs/biomass_next_steps_20261102';s=pd.read_csv(base/'uncertainty/size_observation_registry.csv');v=pd.read_csv(base/'borrowing/visit_support.csv');v=v[v.tolerance_days.eq(7)]
    counts=pd.read_csv(base/'density/biomass_roster_count_units.csv');off=pd.read_csv(ROOT/'outputs/all_program_station_refits_20261102/offsets.csv').set_index('series_id').offset
    self.assertFalse(counts.duplicated(['biomass_series_id','physical_unit_id']).any());self.assertTrue(counts.density.ge(0).all())
    ix=pd.read_csv(ROOT/'outputs/joint_borrowed_biomass_20261102/fit_index.csv')
    for row in ix[ix.rho.eq(0)].itertuples():
      with self.subTest(series=row.series_id),np.load(ROOT/'outputs/joint_borrowed_biomass_20261102/fits'/f'{row.key}.npz') as saved:
        yr=saved['years'];vv=v[v.series_id.eq(row.series_id)&v.year.isin(yr)].reset_index(drop=True);ss=s[s.series_id.eq(row.series_id)].reset_index(drop=True);cc=counts[counts.biomass_series_id.eq(row.series_id)]
        si=int(row.key.split('_')[0][1:]);d=m.build(ss,vv,yr,off[row.series_id],cc,np.random.default_rng(policy['seed']+si*1000),0,pool_df=policy['pool_df'])
        for k in ['y','X','S','K']:np.testing.assert_allclose(d[k],saved[k],atol=1e-11)
        # Vectorized finite differences of each logmass coordinate using aggregate coefficients.
        bm=vv.density.to_numpy()*d['mass'];tot=np.einsum('jri,i->jr',d['A'],bm);eps=1e-5
        def moved(e):
          new=tot[:,:,None]+d['A']*bm[None,None,:]*np.expm1(e)
          return np.log((new[:,0,:]+off[row.series_id])/(new[:,1,:]+off[row.series_id]))
        np.testing.assert_allclose(d['W'],(moved(eps)-moved(-eps))/(2*eps),atol=1e-9)
    smoke=Path('/tmp/joint_biomass_smoke/fits')
    for path in smoke.glob('*.npz'):
      with np.load(path) as z:
        C=z['S']+np.einsum('k,kij->ij',z['q'],z['K'])+np.eye(len(z['y']))*1e-9;b,bc,_=independent(z['y'],z['X'],C)
        np.testing.assert_allclose(b,z['beta'],atol=1e-8);np.testing.assert_allclose(bc,z['beta_cov'],atol=1e-8)
 def test_count_resampling_independent_reconstruction(self):
    _,v,yr,off,c=self.fix;h,W,A=m.aggregation(v,yr,v.mass.to_numpy(),off)
    S,rep=m.count_cov(v,A,v.mass.to_numpy(),off,c,np.random.default_rng(9),draws=70)
    units=sorted(c.physical_unit_id.unique());ww=np.random.default_rng(9).exponential(size=(len(units),70));lookup={u:w for u,w in zip(units,ww)};dd=[]
    for r in v.itertuples():
      g=c[c.site.eq(r.site)&c.year.eq(r.year)&c.zone.eq(r.zone)&c.date.eq(r.date)];weight=np.stack([lookup[u] for u in g.physical_unit_id]);dd.append((g.density.to_numpy()[:,None]*weight).sum(axis=0)/weight.sum(axis=0))
    expected=[]
    for i in range(70):
      vv=v.assign(density=np.array(dd)[:,i]);hh,_,_=m.aggregation(vv,yr,v.mass.to_numpy(),off);expected.append(hh)
    np.testing.assert_allclose(rep,np.array(expected).T,atol=1e-12)
 def test_guard_duplicate_raw_size(self):
    s,v,yr,off,c=self.fix;s=pd.concat([s,s.iloc[[0]]],ignore_index=True)
    with self.assertRaises((ValueError,AssertionError)):build((s,v,yr,off,c))
 def test_guard_duplicate_year(self):
    s,v,yr,off,c=self.fix
    with self.assertRaises((ValueError,AssertionError)):build((s,v,np.r_[yr,yr[0]],off,c))
 def test_guard_invalid_rho(self):
    with self.assertRaises((ValueError,AssertionError)):build(rho=1.2)
 def test_same_count_covariance_across_rho_scenarios(self):
    self.completed_run()
    main=ROOT/'outputs/joint_borrowed_biomass_20261102/fits'
    with np.load(main/'s0_rho0.npz') as z,np.load(main/'s0_rho0.3.npz') as q:np.testing.assert_array_equal(z['S'],q['S'])

if __name__=='__main__':
 print('SOURCE',m.__file__,hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest(),flush=True)
 suite=unittest.defaultTestLoader.loadTestsFromTestCase(JointTests)
 result=unittest.TextTestRunner(verbosity=2).run(suite)
 import json
 record=dict(tests=result.testsRun,skipped=len(result.skipped),failures=len(result.failures),errors=len(result.errors),success=result.wasSuccessful() and not result.skipped,source_sha256=hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest())
 (ROOT/'outputs/joint_borrowed_biomass_20261102/independent_checks.json').write_text(json.dumps(record,indent=2))
 sys.exit(0 if record['success'] else 1)
