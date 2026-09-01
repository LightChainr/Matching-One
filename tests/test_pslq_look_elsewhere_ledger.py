import copy,json,sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"scripts"))
from pslq_look_elsewhere_ledger import DEFAULT_OUTPUT,build_result,mobius,primitive_polynomial_count,validate_result
class LookElsewhereLedgerTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.result=build_result()
 def test_mobius_and_small_counts(self):self.assertEqual([mobius(i) for i in range(1,7)],[1,-1,-1,0,-1,1]);self.assertEqual(primitive_polynomial_count(1,1),3)
 def test_committed_counts_recover_prior_searches(self):self.assertEqual(self.result["primitive_polynomial_counts_by_degree"]["1"],12175);self.assertEqual(self.result["primitive_polynomial_counts_by_degree"]["2"],3355121)
 def test_committed_result(self):
  checked=json.loads(DEFAULT_OUTPUT.read_text());self.assertEqual(checked,self.result);validate_result(checked)
 def test_tampering_fails(self):
  changed=copy.deepcopy(self.result);changed["method_interval_count"]=1
  with self.assertRaisesRegex(ValueError,"does not exactly reproduce"):validate_result(changed)
if __name__=="__main__":unittest.main()
