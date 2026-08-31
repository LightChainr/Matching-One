import importlib.util,json
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[1];PATH=ROOT/"scripts"/"terminal_partition_serial_powers.py";ART=ROOT/"analysis"/"terminal_partition_serial_powers_certificate.json"
SPEC=importlib.util.spec_from_file_location("serial_powers",PATH);M=importlib.util.module_from_spec(SPEC);assert SPEC.loader;SPEC.loader.exec_module(M)
class SerialPowerTests(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.value=M.build_artifact();c.product=M.table()
 def test_artifact(self):self.assertEqual(json.loads(ART.read_text()),self.value)
 def test_histogram(self):self.assertEqual(self.value["index_period_histogram"],{"1,1":12,"1,2":1,"2,1":2})
 def test_period_two(self):self.assertEqual(self.value["profiles"][8]["distinct_powers"],[8,6])
 def test_collapses(self):self.assertEqual(self.value["profiles"][10]["repeat_target"],14);self.assertEqual(self.value["profiles"][11]["repeat_target"],14)
 def test_invalid_index(self):
  with self.assertRaises(ValueError):M.power_profile(15,self.product)
 def test_tamper(self):
  x=json.loads(json.dumps(self.value));x["global_eventual_period"]=1
  with self.assertRaises(ValueError):M.validate_artifact(x)
if __name__=="__main__":unittest.main()
