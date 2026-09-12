"""Finite mathematical controls; none checks document wording or file hashes."""
import importlib.util
from pathlib import Path
import unittest

_spec = importlib.util.spec_from_file_location(
    "axial_ring_gluing", Path(__file__).resolve().parents[1] / "scripts" / "axial_ring_gluing.py")
ring = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ring)


class AxialRingTests(unittest.TestCase):
    def test_four_by_two_ring_and_winding_are_different_events(self):
        """Prevents replacing the sufficient ring event (17) by all winding sets (35)."""
        result = ring.exact_census(4, 1)
        self.assertEqual(sum(result["good_counts_by_occupation"]), 17)
        self.assertEqual(sum(result["winding_counts_by_occupation"]), 35)
        self.assertEqual(result["implication_failures"], 0)

    def test_lift_detects_the_periodic_edge_not_just_connectivity(self):
        """An occupied row winds; deleting one vertex breaks its cycle."""
        width = 7
        self.assertTrue(ring.horizontal_winding((1 << width)-1, width, 1))
        self.assertFalse(ring.horizontal_winding((1 << width)-2, width, 1))

    def test_nondivisible_circumference_does_not_leave_a_seam_gap(self):
        """Prevents dropping the remainder cell when circumference is not a scale multiple."""
        self.assertEqual(ring.cell_boundaries(11, 2), [0, 3, 5, 7, 9, 11])
        result = ring.uneven_cell_controls()
        self.assertEqual(result["implication_failures"], 0)
        self.assertGreater(result["ring_positive"], 0)


if __name__ == "__main__":
    unittest.main()
