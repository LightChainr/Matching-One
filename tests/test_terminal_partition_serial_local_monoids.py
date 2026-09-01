import importlib.util,json
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[1];PATH=ROOT/"scripts"/"terminal_partition_serial_local_monoids.py";ART=ROOT/"analysis"/"terminal_partition_serial_local_monoids_certificate.json"
SPEC=importlib.util.spec_from_file_location("serial_local_monoids",PATH);M=importlib.util.module_from_spec(SPEC);assert SPEC.loader;SPEC.loader.exec_module(M)
class SerialLocalMonoidTests(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.value=M.build_artifact();c.product=M.table()
 def test_artifact(self):self.assertEqual(json.loads(ART.read_text()),self.value)
 def test_size_profile(self):self.assertEqual(self.value["local_size_histogram"],{"1":4,"2":7,"15":1})
 def test_global_local(self):
  record=next(x for x in self.value["local_monoids"] if x["identity"]==6);self.assertEqual(record["carrier"],list(range(15)));self.assertEqual(record["units"],[6,8])
 def test_proper_locals_commute(self):self.assertTrue(all(x["commutative"] for x in self.value["local_monoids"] if x["size"]<15))
 def test_nonidempotent_identity_fails(self):
  with self.assertRaises(ValueError):M.local_record(8,self.product)
 def test_tamper(self):
  x=json.loads(json.dumps(self.value));x["local_monoids"][0]["size"]=99
  with self.assertRaises(ValueError):M.validate_artifact(x)
if __name__=="__main__":unittest.main()
