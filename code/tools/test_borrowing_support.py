"""Independent regression tests for portable rebuild and generated CSVs."""
import json, unittest, os
from pathlib import Path
import numpy as np
import pandas as pd
from rebuild_borrowing_support import fit, aggregate, Engine
O=Path(os.environ.get('BORROWING_TEST_OUTPUT','outputs/biomass_next_steps_20261102/borrowing'))

class BorrowingTests(unittest.TestCase):
    def test_calendar_hac_pairwise(self):
        z=pd.DataFrame(dict(year=[2001,2002,2006,2007,2010,2013],t=[1,2,6,7,10,13],biomass_lnratio=[.2,.7,.4,1.2,.8,1.1]))
        ans=fit(z);X=np.column_stack([np.ones(len(z)),z.t]);y=z.biomass_lnratio.to_numpy()
        beta=np.linalg.lstsq(X,y,rcond=None)[0];u=y-X@beta
        # Pairwise calendar distances, independent of zero-filled score implementation.
        W=np.maximum(0,1-abs(z.year.to_numpy()[:,None]-z.year.to_numpy()[None,:])/4)
        inv=np.linalg.inv(X.T@X)
        V=inv@X.T@((u[:,None]*u[None,:])*W)@X@inv*len(z)/(len(z)-2)
        self.assertAlmostEqual(ans['estimate'],beta[1],places=12)
        self.assertAlmostEqual(ans['se'],np.sqrt(V[1,1]),places=12)
        sequential=z.copy();sequential.year=np.arange(6)
        self.assertGreater(abs(fit(sequential)['se']-ans['se']),1e-6)

    def test_no_skip_unknown_and_equal_hierarchy(self):
        base=pd.DataFrame([dict(series_id='s',year=2010,t=1,revised=0.,original=0.)])
        v=pd.DataFrame([
            ('s',2010,'mpa','A','INNER',10.),('s',2010,'mpa','A','INNER',20.),
            ('s',2010,'mpa','A','MID',30.),('s',2010,'reference','B','INNER',2.),
            ('s',2010,'reference','C','INNER',8.)],
            columns=['series_id','year','status','site','zone','biomass_point'])
        a=aggregate(v,base,{'s':1.}).iloc[0]
        self.assertEqual(a.protected_biomass,22.5);self.assertEqual(a.reference_biomass,5.)
        v.loc[4,'biomass_point']=np.nan;a=aggregate(v,base,{'s':1.}).iloc[0]
        self.assertFalse(a.supported);self.assertTrue(np.isnan(a.reference_biomass))

    def test_actual_outputs_lossless_and_uniform(self):
        v=pd.read_csv(O/'visit_support.csv');a=pd.read_csv(O/'donor_attempts.csv')
        self.assertEqual(len(v),4*1086)
        need=v[v.method.isin(['unsupported','borrowed_other_site','borrowed_other_visit'])]
        self.assertEqual(len(a),2*len(need))
        self.assertTrue(a.groupby(['scenario','visit_key']).kind.nunique().eq(2).all())
        for r in a.itertuples():
            expected=(r.donor_sample_groups>0 and r.cal_clusters>=2 and r.max_gap<=r.cal_max_gap+r.tolerance_days)
            self.assertEqual(bool(r.eligible),bool(expected))
            for d in json.loads(r.donor_evidence):
                dk=json.loads(d['sample_key']);tk=json.loads(r.visit_key)
                self.assertEqual(dk[0],tk[0]);self.assertEqual(dk[2:5],tk[2:5])
                self.assertEqual(d['site']==r.site,r.kind=='other_visit')
                self.assertEqual(d['gap_days'],abs((pd.Timestamp(d['date'])-pd.Timestamp(r.date)).days))
                self.assertNotEqual(d['sample_key'],r.visit_key)
        for (sc,k),g in a.groupby(['scenario','visit_key']):
            eligible=g[g.eligible].sort_values('priority');selected=g[g.selected]
            self.assertEqual(len(selected),int(len(eligible)>0))
            if len(selected):self.assertEqual(selected.kind.iloc[0],eligible.kind.iloc[0])

    def test_actual_outer_exclusion(self):
        c=pd.read_csv(O/'outer_site_calibration.csv')
        for r in c.itertuples():
            self.assertNotEqual(r.outer_site,r.cal_site)
            for donor in json.loads(r.donor_keys):
                d=json.loads(donor)
                self.assertNotEqual(d[1],r.outer_site)
                self.assertEqual(d[0],r.series_id);self.assertEqual(d[2],r.status)
                self.assertEqual(int(d[3]),r.cal_year);self.assertEqual(d[4],r.cal_zone)

    def test_variance_and_decomposition(self):
        v=pd.read_csv(O/'visit_support.csv');m=v[(v.tolerance_days==0)&(v.method=='measured_same_visit')]
        self.assertEqual(int(m.singleton.sum()),84)
        self.assertTrue(m.loc[m.singleton,'empirical_sample_variance'].isna().all())
        self.assertTrue(m.loc[m.singleton,'zero_empirical_variance'].all())
        e=pd.read_csv(O/'effects.csv')
        np.testing.assert_allclose(e.total_change,e.value_change+e.coverage_change,atol=1e-12)
        a=pd.read_csv(O/'annual.csv')
        self.assertTrue((a.supported==(a.unknown_visits==0)).all())
        self.assertTrue(a.loc[a.unknown_reference_visits>0,'reference_biomass'].isna().all())

if __name__=='__main__': unittest.main(verbosity=2)
