import copy,json,sys,unittest
from fractions import Fraction
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"scripts"))
from pslq_synthetic_false_positive_calibration import DEFAULT_OUTPUT,build_result,synthetic_values,validate_result
class SyntheticCalibrationTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.result=build_result()
 def test_generator_is_reproducible_and_bounded(self):
  values=synthetic_values(100,20260830);self.assertEqual(len(values),100);self.assertTrue(all(Fraction(55,100)<=x<=Fraction(65,100) for x in values))
 def test_committed_result(self):
  checked=json.loads(DEFAULT_OUTPUT.read_text());self.assertEqual(checked,self.result);validate_result(checked)
 def test_all_rows_have_exact_diagnostics(self):self.assertEqual(len(self.result["rows"]),100);self.assertTrue(all("/" in row["root_distance"] for row in self.result["rows"]))
 def test_tampering_fails(self):
  changed=copy.deepcopy(self.result);changed["generator"]["seed"]+=1
  with self.assertRaisesRegex(ValueError,"does not exactly reproduce"):validate_result(changed)
if __name__=="__main__":unittest.main()
