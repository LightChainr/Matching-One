import importlib.util,json
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[1];PATH=ROOT/"scripts"/"terminal_partition_serial_idempotents.py";ART=ROOT/"analysis"/"terminal_partition_serial_idempotents_certificate.json"
SPEC=importlib.util.spec_from_file_location("serial_idempotents",PATH);M=importlib.util.module_from_spec(SPEC);assert SPEC.loader;SPEC.loader.exec_module(M)


class SerialIdempotentTests(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.value=M.build_artifact();c.product=M.table()
 def test_artifact(self):self.assertEqual(json.loads(ART.read_text()),self.value)
 def test_sector(self):self.assertEqual(self.value["excluded_indices"],[8]);self.assertEqual(len(self.value["idempotent_generated_sector"]),14)
 def test_rank(self):self.assertEqual(self.value["idempotent_rank"],3);self.assertEqual(len(self.value["minimum_idempotent_generating_sets"]),3)
 def test_all_minimum_sets_close(self):
  target=frozenset(self.value["idempotent_generated_sector"])
  for seed in self.value["minimum_idempotent_generating_sets"]:self.assertEqual(M.closure(seed,self.product),target)
 def test_invalid_index(self):
  with self.assertRaises(ValueError):M.closure([15],self.product)


if __name__=="__main__":unittest.main()
