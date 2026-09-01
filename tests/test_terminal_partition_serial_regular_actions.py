import importlib.util,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parents[1];P=R/"scripts"/"terminal_partition_serial_regular_actions.py";A=R/"analysis"/"terminal_partition_serial_regular_actions_certificate.json"
S=importlib.util.spec_from_file_location("ra",P);M=importlib.util.module_from_spec(S);assert S.loader;S.loader.exec_module(M)
class T(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.v=M.build_artifact()
 def test_artifact(self):self.assertEqual(json.loads(A.read_text()),self.v)
 def test_rank_histograms(self):self.assertEqual(self.v["left_rank_histogram"],{"2":4,"5":9,"15":2});self.assertEqual(self.v["left_rank_histogram"],self.v["right_rank_histogram"])
 def test_intersection(self):self.assertEqual(self.v["left_right_intersection"],[[6,6]])
 def test_width_drift(self):
  with self.assertRaises(ValueError):M.compose([0],[0,1])
 def test_tamper(self):
  x=json.loads(json.dumps(self.v));x["left_right_intersection"]=[]
  with self.assertRaises(ValueError):M.validate_artifact(x)
if __name__=="__main__":unittest.main()
