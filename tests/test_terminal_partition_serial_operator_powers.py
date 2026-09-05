import importlib.util,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parents[1];P=R/"scripts"/"terminal_partition_serial_operator_powers.py";A=R/"analysis"/"terminal_partition_serial_operator_powers_certificate.json";S=importlib.util.spec_from_file_location("m",P);M=importlib.util.module_from_spec(S);assert S.loader;S.loader.exec_module(M)


class T(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.v=M.build_artifact()
 def test_artifact(self):self.assertEqual(json.loads(A.read_text()),self.v)
 def test_profile(self):self.assertEqual(self.v["index_period_histogram"],{"1,1":76,"1,2":17,"2,1":40})
 def test_bounds(self):self.assertEqual((self.v["maximum_index"],self.v["maximum_period"]),(2,2))


if __name__=="__main__":unittest.main()
