import importlib.util,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parents[1];P=R/"scripts"/"terminal_partition_serial_subgroups.py";A=R/"analysis"/"terminal_partition_serial_subgroups_certificate.json"
S=importlib.util.spec_from_file_location("sg",P);M=importlib.util.module_from_spec(S);assert S.loader;S.loader.exec_module(M)
class T(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.v=M.build_artifact()
 def test_artifact(self):self.assertEqual(json.loads(A.read_text()),self.v)
 def test_count(self):self.assertEqual(self.v["subgroup_count"],13);self.assertEqual(self.v["size_histogram"],{"1":12,"2":1})
 def test_unique_nontrivial(self):self.assertEqual([x for x in self.v["subgroups"] if len(x["carrier"])>1],[{"carrier":[6,8],"identity":6}])
 def test_tamper(self):
  x=json.loads(json.dumps(self.v));x["subgroup_count"]=12
  with self.assertRaises(ValueError):M.validate_artifact(x)
if __name__=="__main__":unittest.main()
