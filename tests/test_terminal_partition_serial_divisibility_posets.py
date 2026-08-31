import importlib.util,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parents[1];P=R/"scripts"/"terminal_partition_serial_divisibility_posets.py";A=R/"analysis"/"terminal_partition_serial_divisibility_posets_certificate.json";S=importlib.util.spec_from_file_location("m",P);M=importlib.util.module_from_spec(S);assert S.loader;S.loader.exec_module(M)
class T(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.v=M.build_artifact()
 def test_artifact(self):self.assertEqual(json.loads(A.read_text()),self.v)
 def test_counts(self):self.assertEqual((self.v["left"]["class_count"],self.v["right"]["class_count"],self.v["two_sided"]["class_count"]),(6,6,3))
 def test_chain(self):self.assertEqual(self.v["two_sided"]["cover_relations"],[[0,1],[1,2]])
 def test_tamper(self):
  x=dict(self.v);x["status"]="bad"
  with self.assertRaises(ValueError):M.validate_artifact(x)
if __name__=="__main__":unittest.main()
