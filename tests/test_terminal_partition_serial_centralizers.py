import importlib.util,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parents[1];P=R/"scripts"/"terminal_partition_serial_centralizers.py";A=R/"analysis"/"terminal_partition_serial_centralizers_certificate.json";S=importlib.util.spec_from_file_location("m",P);M=importlib.util.module_from_spec(S);assert S.loader;S.loader.exec_module(M)
class T(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.v=M.build_artifact()
 def test_artifact(self):self.assertEqual(json.loads(A.read_text()),self.v)
 def test_profile(self):self.assertEqual(self.v["centralizer_size_histogram"],{"3":6,"4":4,"5":2,"7":2,"15":1})
 def test_center(self):self.assertEqual(self.v["center"],[6]);self.assertEqual(self.v["ordered_commuting_pair_count"],73)
 def test_tamper(self):
  x=dict(self.v);x["center"]=[]
  with self.assertRaises(ValueError):M.validate_artifact(x)
if __name__=="__main__":unittest.main()
