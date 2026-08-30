from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_n65_charged_activity_net import score, symmetric_eigen  # noqa: E402


class N65ChargedActivityNetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = score(
            ROOT / "results/p337-n65-charged-source-reveal/latest.json"
        )

    def test_jacobi_eigenmodes_reconstruct_trace(self) -> None:
        matrix = [[1.0, 0.4, 0.1], [0.4, 2.0, -0.2], [0.1, -0.2, 3.0]]
        values, vectors = symmetric_eigen(matrix)
        self.assertAlmostEqual(sum(values), 6.0)
        for value, vector in zip(values, vectors):
            image = [sum(matrix[i][j] * vector[j] for j in range(3))
                     for i in range(3)]
            self.assertLess(max(abs(image[i] - value * vector[i]) for i in range(3)), 1e-12)

    def test_net_activity_transform_preserves_full_quadratic(self) -> None:
        a = self.payload["orientation_contrast"]["A"]["quadratic"]
        d = self.payload["orientation_contrast"]["D"]["quadratic"]
        self.assertAlmostEqual(a["full"], 12.153148777336366)
        self.assertAlmostEqual(d["full"], 1.5091860895963012)
        self.assertGreater(a["net_only"], 11.0)
        self.assertLess(a["activity_only"], 0.3)
        self.assertLess(d["net_only"], 1.4)

    def test_eigenmode_contributions_sum_to_quadratic(self) -> None:
        for channel in ("A", "D"):
            block = self.payload["orientation_contrast"][channel]
            contribution = sum(mode["quadratic_contribution"]
                               for mode in block["correlation_eigenmodes"])
            self.assertAlmostEqual(contribution, block["quadratic"]["full"])

    def test_dependency_groups_are_not_pooled(self) -> None:
        crosswalk = self.payload["p334_crosswalk"]
        self.assertIn("not independent", crosswalk["same_archive_coordinates"])
        self.assertIn("no quadratic", crosswalk["different_dependency_group"])


if __name__ == "__main__":
    unittest.main()
