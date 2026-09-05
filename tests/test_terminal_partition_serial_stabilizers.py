import importlib.util,json
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[1];PATH=ROOT/"scripts"/"terminal_partition_serial_stabilizers.py";ART=ROOT/"analysis"/"terminal_partition_serial_stabilizers_certificate.json"
SPEC=importlib.util.spec_from_file_location("serial_stabilizers",PATH);M=importlib.util.module_from_spec(SPEC);assert SPEC.loader;SPEC.loader.exec_module(M)


class SerialStabilizerTests(unittest.TestCase):
 @classmethod
 def setUpClass(c):c.value=M.build_artifact()
 def test_artifact(self):self.assertEqual(json.loads(ART.read_text()),self.value)
 def test_cancellation(self):self.assertEqual(self.value["left_cancellative_indices"],[6,8]);self.assertEqual(self.value["right_cancellative_indices"],[6,8])
 def test_translation_ranks(self):self.assertEqual(self.value["left_translation_ranks"],[5,5,5,2,2,5,15,5,15,5,5,5,5,2,2]);self.assertEqual(self.value["left_translation_ranks"],self.value["right_translation_ranks"])
 def test_identity_stabilizer(self):self.assertEqual(self.value["left_stabilizers"][6],[6]);self.assertEqual(self.value["right_stabilizers"][6],[6])


if __name__=="__main__":unittest.main()
