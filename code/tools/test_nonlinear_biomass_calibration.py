"""Independent algebra, generator-moment and scope checks for nonlinear review."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json,hashlib,unittest
import numpy as np,pandas as pd
import nonlinear_biomass_calibration as m
from joint_borrowed_biomass import build,fit


def fixture():
    rows=[];cc=[]
    for role in ['mpa','reference']:
        for site in [role+'A',role+'B']:
            for year in [2013,2015,2018,2021]:
                for zone in ['INNER','MID']:
                    i=len(rows);date=f'{year}-07-01';density=.2+(i%7)/5;mass=8+i%11
                    rows.append(dict(status=role,site=site,year=year,date=date,zone=zone,density=density,mass=mass,mean_mass=mass,relative_ss=1.2,df=3.,n=4.,sampling_cv2=.4,method='measured_same_visit'))
                    for tr,k in enumerate([.5,1.,1.5]):cc.append(dict(site=site,year=year,date=date,zone=zone,method='swath',density=density*k,physical_unit_id=f'{i}_{tr}'))
    s=pd.DataFrame(rows);v=s.copy();years=np.array([2013,2015,2018,2021]);counts=pd.DataFrame(cc);d=build(s,v,years,.5,counts,np.random.default_rng(1));s['sampling_cv2']=np.expm1(d['known'])*s.n
    f=fit(d['y'],d['X'],d['S'],d['K']);return s,v,years,counts,d,f

class Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s,cls.v,cls.years,cls.counts,cls.d,cls.f=fixture()
    def test_exact_arithmetic_by_group_means(self):
        d=self.d;v=self.v;rng=np.random.default_rng(3);z=rng.normal(2,.3,(len(v),20));dd=m.density_replicates(v,self.counts,rng,20);r=m.exact_ratio(d['A'],dd,z,.5)
        for j in range(20):
            vv=v.assign(biomass=dd[:,j]*np.exp(z[:,j]));a=vv.groupby(['year','status','site','zone']).biomass.mean().groupby(['year','status','site']).mean().groupby(['year','status']).mean().unstack()
            np.testing.assert_allclose(r[:,j],np.log((a.mpa+.5)/(a.reference+.5)))
    def test_exact_equals_linear_in_single_cell_zero_offset_case(self):
        A=np.zeros((1,2,2));A[0,0,0]=1;A[0,1,1]=1;z=np.random.default_rng(5).normal(size=(2,100));density=np.array([[2.],[3.]])
        np.testing.assert_allclose(m.exact_ratio(A,density,z,0)[0],np.log(2/3)+z[0]-z[1],atol=1e-14)
    def test_zero_counts_remain_zero(self):
        A=np.zeros((1,2,2));A[0,0,0]=1;A[0,1,1]=1;z=np.ones((2,3))*20;assert np.isfinite(m.exact_ratio(A,np.zeros((2,3)),z,.5)).all();np.testing.assert_array_equal(m.exact_ratio(A,np.zeros((2,3)),z,.5),np.zeros((1,3)))
    def test_reexpression_preserves_nonlinear_residual(self):
        d=self.d;z=d['logm']+.2*np.sin(np.arange(len(self.v)));dd=m.reexpress(d,d['y'],self.v,self.years,.5,z)
        c0=d['h']-d['W']@d['logm'];c1=dd['h']-dd['W']@dd['logm']
        np.testing.assert_allclose(dd['y'][d['ns']:]-c1,d['y'][d['ns']:]-c0,atol=1e-12)
        for k in dd['K']:assert np.linalg.eigvalsh(k).min()>-1e-8
    def test_paired_gaussian_moments_and_curvature_identity(self):
        n=16000;d=self.d;f=self.f
        Y,detail=m.simulated_vectors(d,f,self.v,self.s,self.counts,.5,np.random.default_rng(9),n)
        expected=d['X']@f['beta'];sd=np.sqrt(np.diag(f['C']));z=(Y['gaussian'].mean(axis=1)-expected)/(sd/np.sqrt(n));assert max(abs(z))<5
        C=np.cov(Y['gaussian']);assert np.linalg.norm(C-f['C'])/np.linalg.norm(f['C'])<.07
        np.testing.assert_allclose(Y['exact_mass'][:d['ns']],Y['gaussian'][:d['ns']])
        np.testing.assert_allclose(Y['exact_mass'][d['ns']:]-Y['gaussian'][d['ns']:],-detail['curvature'],atol=1e-12)
        relative=np.exp(detail['size_error']-.5*d['known'][:,None]);assert max(abs(relative.mean(axis=1)-1))<.025
    def test_saved_support_and_hashes(self):
        O=Path('outputs/nonlinear_biomass_calibration_20261102')
        if not (O/'completed.json').exists():self.skipTest('Review run incomplete')
        support=pd.read_csv(O/'support.csv');old=pd.read_csv('outputs/joint_borrowed_biomass_20261102/support.csv').set_index('series_id')
        for r in support.itertuples():assert r.years==old.loc[r.series_id,'years']
        for p,h in json.loads((O/'input_hashes.json').read_text()).items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
        fits=pd.read_csv(O/'paired_fits.csv');assert np.isfinite(fits[['estimate','se','pivot']]).all().all();assert (fits.se>0).all();assert fits.groupby(['series_id','rho','replicate']).size().eq(3).all()
if __name__=='__main__':
    result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Tests));O=Path('outputs/nonlinear_biomass_calibration_20261102');O.mkdir(exist_ok=True,parents=True)
    record=dict(success=result.wasSuccessful() and not result.skipped,tests=result.testsRun,skipped=len(result.skipped),failures=len(result.failures),errors=len(result.errors));(O/'checks.json').write_text(json.dumps(record,indent=2))
    raise SystemExit(0 if record['success'] else 1)
