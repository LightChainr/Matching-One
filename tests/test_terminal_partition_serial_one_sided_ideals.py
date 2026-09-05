import importlib.util,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parents[1];P=R/"scripts"/"terminal_partition_serial_one_sided_ideals.py";A=R/"analysis"/"terminal_partition_serial_one_sided_ideals_certificate.json"
S=importlib.util.spec_from_file_location("oi",P);M=importlib.util.module_from_spec(S);assert S.loader;S.loader.exec_module(M)


class T(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.v=M.build_artifact()
 def test_artifact(self):self.assertEqual(json.loads(A.read_text()),self.v)
 def test_counts(self):self.assertEqual(len(self.v["left_ideals"]),16);self.assertEqual(len(self.v["right_ideals"]),16)
 def test_covers_and_height(self):self.assertEqual(len(self.v["left_cover_relations"]),26);self.assertEqual(self.v["left_height"],7);self.assertEqual(self.v["right_height"],7)
 def test_histograms(self):self.assertEqual(self.v["left_size_histogram"],self.v["right_size_histogram"])


if __name__=="__main__":unittest.main()
