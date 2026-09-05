import importlib.util,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parents[1];P=R/"scripts"/"terminal_partition_serial_two_sided_actions.py";A=R/"analysis"/"terminal_partition_serial_two_sided_actions_certificate.json"
S=importlib.util.spec_from_file_location("ta",P);M=importlib.util.module_from_spec(S);assert S.loader;S.loader.exec_module(M)


class T(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.v=M.build_artifact();c.p=M.table()
 def test_artifact(self):self.assertEqual(json.loads(A.read_text()),self.v)
 def test_counts(self):self.assertEqual(self.v["distinct_operator_count"],133);self.assertEqual(self.v["rank_histogram"],{"1":4,"2":89,"5":36,"15":4})
 def test_fibers(self):self.assertEqual(self.v["pair_fiber_size_histogram"],{"1":121,"2":8,"16":1,"22":2,"28":1})
 def test_bad_multiplier(self):
  with self.assertRaises(ValueError):M.operator(15,0,self.p)


if __name__=="__main__":unittest.main()
