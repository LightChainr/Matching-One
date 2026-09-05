import json,sys,unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"scripts"))
from standard_constant_stability_audit import DEFAULT_OUTPUT,build_result


class StandardConstantStabilityTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.result=build_result()
 def test_committed_result(self):
  # One rebuild only: validate_result() would rebuild a second time for the same check.
  checked=json.loads(DEFAULT_OUTPUT.read_text());self.assertEqual(checked,self.result)
  self.assertEqual(sum(len(r["endpoint_checks"]) for r in self.result["rows"]),144)
 def test_all_witnesses_and_endpoints_stable(self):self.assertTrue(self.result["conclusion"]["all_closest_witnesses_stable"]);self.assertTrue(self.result["conclusion"]["all_144_endpoint_checks_exclude_zero"])


if __name__=="__main__":unittest.main()
