import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'tools'))
from biomass_source_guard import require_allowed
class Guard(unittest.TestCase):
 def test_aliases_blocked(self):
  for s in ['a-k-01','a-k-14','a-k-21','a-k-37']:
   with self.assertRaises(ValueError):require_allowed(s,'M. franciscanus')
 def test_known_valid_and_distinct_table(self):
  require_allowed('a-k-33','S. purpuratus')
  require_allowed('a-k-01','M. pyrifera','KFM_Macrocystis_RawData_1984-2023.csv')
 def test_explicit_exclusion_applies_all_tables(self):
  with self.assertRaises(ValueError):require_allowed('a-k-21','M. pyrifera','original')
if __name__=='__main__':unittest.main()
