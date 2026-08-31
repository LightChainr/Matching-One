import importlib.util,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parents[1];P=R/"scripts"/"terminal_partition_serial_subsemigroup_lattice.py";A=R/"analysis"/"terminal_partition_serial_subsemigroup_lattice_certificate.json";S=importlib.util.spec_from_file_location("m",P);M=importlib.util.module_from_spec(S);assert S.loader;S.loader.exec_module(M)
class T(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.v=M.build_artifact()
 def test_artifact(self):self.assertEqual(json.loads(A.read_text()),self.v)
 def test_invariants(self):self.assertEqual((self.v["elements"],self.v["cover_count"],self.v["height"],self.v["width"]),(416,1400,12,82))
 def test_dilworth(self):self.assertEqual(self.v["elements"]-self.v["dilworth_matching_size"],self.v["width"])
 def test_tamper(self):
  x=dict(self.v);x["width"]=81
  with self.assertRaises(ValueError):M.validate_artifact(x)
if __name__=="__main__":unittest.main()
