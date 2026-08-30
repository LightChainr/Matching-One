import json
from pathlib import Path
import unittest

from p334_tm_translation_orbit_hall import build_result


ROOT = Path(__file__).resolve().parents[1]


class P334TMTranslationOrbitHallTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build_result()

    def test_semantic_correction_is_explicit(self):
        correction = self.result["semantic_correction"]
        self.assertIn("fixed projective line", correction["missing_gate"])
        self.assertIn("588", correction["effect"])

    def test_minimal_rows_need_both_resources(self):
        rows = self.result["minimal_HNF_obstruction"]["rows"]
        self.assertEqual(len(rows), 4)
        for row in rows:
            self.assertEqual(
                row["two_carrier_base_transport"]["Hall_deficiency"], 564
            )
            self.assertEqual(
                row["transverse_mark_only"]["Hall_deficiency"], 384
            )
            self.assertEqual(
                row["one_carrier_plus_transverse_mark"]["Hall_deficiency"], 0
            )

    def test_orbit_compression_closes_first_smith_gate(self):
        row = self.result["next_Smith_gate"]["row"]
        self.assertEqual(row["Smith_invariants"], [2, 4])
        old = row["two_carrier_base_transport_orbits"]
        new = row["one_carrier_plus_transverse_mark_orbits"]
        self.assertEqual(old["compression_factor"], 8)
        self.assertEqual(old["Hall_deficiency"], 3264)
        self.assertEqual(new["maximum_matching"], 5760)
        self.assertEqual(new["Hall_deficiency"], 0)

    def test_checked_artifact_reproduces(self):
        checked = json.loads(
            (
                ROOT
                / "results/p334-tm-translation-orbit-hall/latest.json"
            ).read_text()
        )
        self.assertEqual(checked, self.result)


if __name__ == "__main__":
    unittest.main()
