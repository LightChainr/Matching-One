import importlib.util,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parents[1];P=R/"scripts"/"terminal_partition_serial_operator_commuting.py";A=R/"analysis"/"terminal_partition_serial_operator_commuting_certificate.json";S=importlib.util.spec_from_file_location("m",P);M=importlib.util.module_from_spec(S);assert S.loader;S.loader.exec_module(M)


class T(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.v=M.build_artifact()
 def test_artifact(self):self.assertEqual(json.loads(A.read_text()),self.v)
 def test_center(self):self.assertEqual(self.v["center"],[13])
 def test_count(self):self.assertEqual(self.v["ordered_commuting_pair_count"],3239);self.assertEqual(sum(self.v["centralizer_size_histogram"].values()),133)


if __name__=="__main__":unittest.main()
