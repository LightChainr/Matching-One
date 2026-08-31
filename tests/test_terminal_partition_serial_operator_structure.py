import importlib.util,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parents[1];P=R/"scripts"/"terminal_partition_serial_operator_structure.py";A=R/"analysis"/"terminal_partition_serial_operator_structure_certificate.json";S=importlib.util.spec_from_file_location("m",P);M=importlib.util.module_from_spec(S);assert S.loader;S.loader.exec_module(M)
class T(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.v=M.build_artifact()
 def test_artifact(self):self.assertEqual(json.loads(A.read_text()),self.v)
 def test_counts(self):self.assertEqual((self.v["operator_count"],self.v["identity_index"],len(self.v["unit_indices"]),self.v["idempotent_count"]),(133,13,4,76))
 def test_ranks(self):self.assertEqual(self.v["transformation_rank_histogram"],{"1":4,"2":89,"5":36,"15":4})
 def test_tamper(self):
  x=dict(self.v);x["idempotent_count"]=75
  with self.assertRaises(ValueError):M.validate_artifact(x)
if __name__=="__main__":unittest.main()
