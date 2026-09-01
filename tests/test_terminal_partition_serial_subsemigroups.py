import importlib.util,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parents[1];P=R/"scripts"/"terminal_partition_serial_subsemigroups.py";A=R/"analysis"/"terminal_partition_serial_subsemigroups_certificate.json"
S=importlib.util.spec_from_file_location("ss",P);M=importlib.util.module_from_spec(S);assert S.loader;S.loader.exec_module(M)
class T(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.v=M.build_artifact()
 def test_artifact(self):self.assertEqual(json.loads(A.read_text()),self.v)
 def test_counts(self):self.assertEqual(self.v["subsemigroup_count"],416);self.assertEqual(self.v["with_wire_identity_count"],228);self.assertEqual(self.v["without_wire_identity_count"],188)
 def test_histogram_mass(self):self.assertEqual(sum(self.v["size_histogram"].values()),416)
 def test_maximal(self):self.assertEqual(sorted(map(len,self.v["maximal_proper_subsemigroups"])),[9,9,12,12,14])
 def test_tamper(self):
  x=json.loads(json.dumps(self.v));x["subsemigroup_count"]=415
  with self.assertRaises(ValueError):M.validate_artifact(x)
if __name__=="__main__":unittest.main()
