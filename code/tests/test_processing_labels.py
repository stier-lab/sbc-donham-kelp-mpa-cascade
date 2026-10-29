import importlib.util
import unittest
from pathlib import Path

path=Path(__file__).resolve().parents[1]/'tools/check_processing_labels.py'
spec=importlib.util.spec_from_file_location('processing_labels',path)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

class ProcessingLabelsTest(unittest.TestCase):
    def test_ci_label_does_not_imply_no_before_data(self):
        row=dict(MPA='Example',Taxa='Example',Source='PISCO',Resp='Den',AnalysisType='CI',Model='Mean')
        coverage={**row,'primary_design':'CI','n_before':'4','n_after':'8'}
        out=module.build_routes([row],[coverage])[0]
        self.assertTrue(out['three_years_each_side'])
        self.assertEqual(out['construction_route'],'PISCO_direct_CI')
        self.assertEqual(out['pretrend_failure_status'],'not_established_by_processing_label')
    def test_duplicate_identity_fails(self):
        row=dict(MPA='Example',Taxa='Example',Source='PISCO',Resp='Den',AnalysisType='CI',Model='Mean')
        coverage={**row,'primary_design':'CI','n_before':'4','n_after':'8'}
        with self.assertRaises(ValueError):module.build_routes([row,row],[coverage])

if __name__=='__main__':unittest.main()
