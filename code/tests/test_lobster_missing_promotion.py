import importlib.util, unittest
from pathlib import Path
p=Path(__file__).resolve().parents[1]/'promote_lobster_missing_lengths.py'
s=importlib.util.spec_from_file_location('promotion',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
class Promotion(unittest.TestCase):
 def setUp(self):
  self.before=b'CA_MPA_Name_Short,year,taxon_name,source,resp,status,value\nA,2014,Panulirus interruptus,LTER,Bio,mpa,1\nA,2014,Panulirus interruptus,LTER,Den,mpa,2\n'
 def test_scoped_value(self):
  after=self.before.replace(b'mpa,1',b'mpa,1.1');out,changes=m.promote(self.before,after)
  self.assertEqual(len(changes),1);self.assertEqual(out.splitlines()[-1],self.before.splitlines()[-1])
 def test_reject_unrelated(self):
  with self.assertRaises(AssertionError):m.promote(self.before,self.before.replace(b'mpa,2',b'mpa,3'))
 def test_reject_membership(self):
  with self.assertRaises(AssertionError):m.promote(self.before,self.before.replace(b'2014',b'2015'))
if __name__=='__main__':unittest.main()
