import importlib.util,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parents[1];P=R/"scripts"/"terminal_partition_serial_operator_inverses.py";A=R/"analysis"/"terminal_partition_serial_operator_inverses_certificate.json";S=importlib.util.spec_from_file_location("m",P);M=importlib.util.module_from_spec(S);assert S.loader;S.loader.exec_module(M)
class T(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.v=M.build_artifact()
 def test_artifact(self):self.assertEqual(json.loads(A.read_text()),self.v)
 def test_regular(self):self.assertEqual(self.v["regular_operator_count"],133)
 def test_profile(self):self.assertEqual(sum(self.v["inverse_count_histogram"].values()),133);self.assertEqual(self.v["inverse_count_histogram"]["81"],1)
 def test_tamper(self):
  x=dict(self.v);x["regular_operator_count"]=132
  with self.assertRaises(ValueError):M.validate_artifact(x)
if __name__=="__main__":unittest.main()
