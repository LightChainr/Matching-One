import importlib.util,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parents[1];P=R/"scripts"/"terminal_partition_serial_unit_orbits.py";A=R/"analysis"/"terminal_partition_serial_unit_orbits_certificate.json";S=importlib.util.spec_from_file_location("m",P);M=importlib.util.module_from_spec(S);assert S.loader;S.loader.exec_module(M)
class T(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.v=M.build_artifact()
 def test_artifact(self):self.assertEqual(json.loads(A.read_text()),self.v)
 def test_counts(self):self.assertEqual((len(self.v["left_orbits"]),len(self.v["right_orbits"]),len(self.v["double_orbits"]),len(self.v["conjugation_orbits"])),(11,11,9,11))
 def test_units(self):self.assertEqual(self.v["unit_indices"],[6,8])
 def test_tamper(self):
  x=dict(self.v);x["unit_indices"]=[6]
  with self.assertRaises(ValueError):M.validate_artifact(x)
if __name__=="__main__":unittest.main()
