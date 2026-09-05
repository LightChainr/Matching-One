
from __future__ import annotations
from pathlib import Path
import shutil
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import validate_p50_fullcurve_quantity_contract as contract  # noqa: E402


class P50FullcurveQuantityContractTests(unittest.TestCase):
    def copied_root(self) -> tuple[tempfile.TemporaryDirectory, Path]:
        directory = tempfile.TemporaryDirectory()
        root = Path(directory.name)
        destination = root / contract.CONTRACT
        destination.parent.mkdir(parents=True)
        shutil.copy(ROOT / contract.CONTRACT, destination)
        return directory, root

    def test_repository_blobs_are_frozen(self) -> None:
        payload, _ = contract.load_contract(ROOT)
        contract.validate_repository_files(ROOT, payload)

    def test_four_topology_anchors_are_exact_identities(self) -> None:
        _, validated = contract.load_contract(ROOT)
        self.assertEqual(
            list(validated), ["thermal_even", "matching_function", "P4_S", "P4_D"]
        )
        self.assertTrue(all(
            (row["transform"].scale, row["transform"].offset) == (1.0, 0.0)
            for row in validated.values()
        ))

    def test_response_coordinates_are_not_descriptor_fields(self) -> None:
        payload, _ = contract.load_contract(ROOT)
        descriptor_fields = set(next(iter(payload["topology_anchors"].values())))
        self.assertTrue(set(payload["response_coordinates"]).isdisjoint(descriptor_fields))
        self.assertIn("root_gap_lineage", payload["response_coordinates"])


if __name__ == "__main__":
    unittest.main()
