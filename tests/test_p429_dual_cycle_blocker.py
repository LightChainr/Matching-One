from __future__ import annotations
import copy
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from verify_p429_dual_cycle_blocker import Geometry, verify, verify_case

class DualCycleBlockerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.payload = json.loads((ROOT / "results/p429-dual-cycle-blocker/certificate.json").read_text())

    def test_both_real_cycle_certificates(self):
        result = verify(self.payload)
        self.assertEqual([x["cycle_lengths"] for x in result], [[20,25],[19,43]])
        self.assertTrue(all(x["two_cycle_certificate"] for x in result))

    def test_full_pair_regeneration(self):
        expected = json.loads((ROOT / "results/p429-dual-cycle-blocker/verification.json").read_text())
        self.assertEqual(verify(self.payload, full_pairs=True), expected)
        self.assertEqual(sum(x["all_safe_pair_checks"] for x in expected), 29756)

    def test_wrong_winding_rejected(self):
        x = copy.deepcopy(self.payload)
        x["cases"][0]["cycles"][0]["winding"] = [0,0]
        with self.assertRaises(ValueError): verify(x)

    def test_duplicate_cycle_vertex_rejected(self):
        x = copy.deepcopy(self.payload)
        x["cases"][0]["cycles"][0]["vertices"].append(x["cases"][0]["cycles"][0]["vertices"][0])
        with self.assertRaises(ValueError): verify(x)

    def test_nonwhite_vertex_rejected(self):
        x = copy.deepcopy(self.payload)
        mask = int(x["cases"][0]["occupied_mask_hex"], 16)
        x["cases"][0]["cycles"][0]["vertices"][0] = next(v for v in range(425) if mask >> v & 1)
        with self.assertRaises(ValueError): verify(x)

    def test_unproven_safe_overlap_rejected(self):
        x = copy.deepcopy(self.payload)
        x["cases"][0]["cycles"][1] = copy.deepcopy(x["cases"][0]["cycles"][0])
        with self.assertRaises(ValueError): verify(x)

    def test_wrong_occupation_count_rejected(self):
        x = copy.deepcopy(self.payload)
        x["cases"][0]["occupied_count"] += 1
        with self.assertRaises(ValueError): verify(x)

    def test_wrong_singleton_set_rejected(self):
        x = copy.deepcopy(self.payload)
        x["cases"][0]["singleton_triggers"] = [27]
        with self.assertRaises(ValueError): verify(x)

    def test_duplicate_case_id_rejected(self):
        x = copy.deepcopy(self.payload)
        x["cases"][1]["id"] = x["cases"][0]["id"]
        with self.assertRaises(ValueError): verify(x)

    def test_unknown_schema_rejected(self):
        x = copy.deepcopy(self.payload); x["schema"] = "unknown"
        with self.assertRaises(ValueError): verify(x)

    def test_non_hnf_period_rejected(self):
        with self.assertRaises(ValueError): Geometry([[4,0],[1,4]])

    def test_overlap_only_at_singleton_is_allowed(self):
        x = {
            "id":"N16-A", "periods":[[4,0],[0,4]], "occupied_mask_hex":hex(12463),
            "occupied_count":8, "singleton_triggers":[9],
            "cycles":[{"vertices":[4,9,10,11],"winding":[1,0]},
                      {"vertices":[8,9,14,15],"winding":[1,0]}],
        }
        self.assertTrue(verify_case(x)["two_cycle_certificate"])
        self.assertEqual(verify_case(x)["singleton_count"],1)

if __name__ == "__main__": unittest.main()
