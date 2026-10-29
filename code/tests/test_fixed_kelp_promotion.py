import importlib.util,unittest
from pathlib import Path
p=Path(__file__).resolve().parents[1]/'promote_fixed_kelp.py'
s=importlib.util.spec_from_file_location('promote',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
class PromotionTests(unittest.TestCase):
 def setUp(self):
  self.base=b'CA_MPA_Name_Short,year,taxon_name,source,resp,status,value\nA,2014,Macrocystis pyrifera,KFM,Bio,Inside,100\nA,2014,Macrocystis pyrifera,KFM,Den,Inside,2\nA,2014,Panulirus interruptus,KFM,Den,Inside,3\n'
 def test_only_biomass_changes(self):
  out,changes=m.promote(self.base,self.base.replace(b',100\n',b',95\n'));self.assertEqual(len(changes),1);self.assertEqual(out.splitlines()[2:],self.base.splitlines()[2:])
 def test_other_change_rejected(self):
  with self.assertRaises(AssertionError):m.promote(self.base,self.base.replace(b',Inside,3',b',Inside,4'))
 def test_identity_change_rejected(self):
  with self.assertRaises(AssertionError):m.promote(self.base,self.base.replace(b'A,2014',b'A,2015'))
if __name__=='__main__':unittest.main()
