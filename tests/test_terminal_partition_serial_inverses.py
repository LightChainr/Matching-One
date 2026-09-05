import importlib.util,json
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[1]; PATH=ROOT/"scripts"/"terminal_partition_serial_inverses.py"; ART=ROOT/"analysis"/"terminal_partition_serial_inverses_certificate.json"
SPEC=importlib.util.spec_from_file_location("serial_inverses",PATH); M=importlib.util.module_from_spec(SPEC); assert SPEC.loader; SPEC.loader.exec_module(M)


class SerialInverseTests(unittest.TestCase):
 @classmethod
 def setUpClass(c): c.value=M.build_artifact()
 def test_artifact(self): self.assertEqual(json.loads(ART.read_text()),self.value)
 def test_all_regular(self): self.assertEqual(self.value["regular_element_indices"],list(range(15)))
 def test_profile(self): self.assertEqual(self.value["inverse_count_profile"],[9,6,6,4,4,6,1,4,1,6,4,4,4,4,4])
 def test_units(self): self.assertEqual(self.value["unit_indices"],[6,8])
 def test_inverse_relation_symmetric(self):
  inv=[set(x) for x in self.value["inverse_sets"]]
  self.assertTrue(all((b in inv[a])==(a in inv[b]) for a in range(15) for b in range(15)))


if __name__=="__main__": unittest.main()
