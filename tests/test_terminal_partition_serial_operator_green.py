import importlib.util,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parents[1];P=R/"scripts"/"terminal_partition_serial_operator_green.py";A=R/"analysis"/"terminal_partition_serial_operator_green_certificate.json";S=importlib.util.spec_from_file_location("m",P);M=importlib.util.module_from_spec(S);assert S.loader;S.loader.exec_module(M)


class T(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.v=M.build_artifact()
 def test_artifact(self):self.assertEqual(json.loads(A.read_text()),self.v)
 def test_counts(self):self.assertEqual(self.v["class_counts"],{"L":21,"R":24,"H":112,"J":7,"D":7})
 def test_d_j(self):self.assertEqual({tuple(x) for x in self.v["D_classes"]},{tuple(x) for x in self.v["J_classes"]})


if __name__=="__main__":unittest.main()
