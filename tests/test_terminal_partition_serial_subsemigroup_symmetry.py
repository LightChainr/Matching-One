import importlib.util,json,unittest
from pathlib import Path
R=Path(__file__).resolve().parents[1];P=R/"scripts"/"terminal_partition_serial_subsemigroup_symmetry.py";A=R/"analysis"/"terminal_partition_serial_subsemigroup_symmetry_certificate.json";S=importlib.util.spec_from_file_location("m",P);M=importlib.util.module_from_spec(S);assert S.loader;S.loader.exec_module(M)


class T(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.v=M.build_artifact()
 def test_artifact(self):self.assertEqual(json.loads(A.read_text()),self.v)
 def test_profile(self):self.assertEqual(self.v["orbit_count"],165);self.assertEqual(self.v["orbit_size_histogram"],{"1":34,"2":71,"4":60})
 def test_stable(self):self.assertEqual((self.v["reversal_stable_count"],self.v["lane_swap_stable_count"]),(54,128))


if __name__=="__main__":unittest.main()
