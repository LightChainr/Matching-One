import copy,json,sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"scripts"))
from kagome_exact_relation_control import DEFAULT_OUTPUT,build_result,validate_result


class KagomeControlTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.result=build_result()
 def test_committed_result(self):
  checked=json.loads(DEFAULT_OUTPUT.read_text());self.assertEqual(checked,self.result);validate_result(checked)
 def test_relation_and_unique_root(self):
  self.assertEqual(self.result["polynomial_coefficients_ascending"],[1,0,-3,1]);self.assertEqual(self.result["sturm_open_root_count_in_physical_window"],1);self.assertEqual(self.result["endpoint_signs"],[1,-1])


if __name__=="__main__":unittest.main()
