"""Tests for the #595 channel-leverage design scan.

Each test names the wrong number it would stop us believing.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p205_channel_leverage import (  # noqa: E402
    BASELINE_CHANNEL,
    BASIS_CHANNELS,
    COMBINATION_FLOOR,
    P_GRID_POINTS,
    P_HALF_WIDTH,
    P_REF,
    SPLIT_PARITY,
    baseline_key,
    best_combination,
    decide,
    grid_boundary_binds,
    probability_grid,
    readouts,
    selection_control,
    symmetric_eigen,
)


class Grid(unittest.TestCase):
    def test_p_ref_is_on_the_grid_exactly(self) -> None:
        """Stops us believing an amplification measured against the wrong baseline.

        Every number this scan exports is a ratio to the historical readout at
        the frozen ``p_ref``.  If the grid missed ``p_ref`` the baseline would be
        silently taken at a neighbouring probability, and the whole cost table
        would be scaled by an unknown factor.
        """

        mp.mp.dps = 40
        for half, points in (("0.035", 21), ("0.10", 21), ("0.001", 3)):
            grid = probability_grid(half, points)
            self.assertEqual(len(grid), points)
            self.assertEqual(baseline_key(grid), mp.nstr(mp.mpf(P_REF), 12))
            self.assertLess(abs(grid[0] + grid[-1] - 2 * mp.mpf(P_REF)), mp.mpf("1e-30"))
        with self.assertRaises(ValueError):
            probability_grid("0.035", 20)

    def test_boundary_binding_is_detected(self) -> None:
        """Stops us believing an optimum that is really the edge of what we looked at.

        A maximum at a grid endpoint means the scan never saw the optimum, so
        the reported amplification is a lower bound.  Missing that would let a
        design recommendation be quoted as if the range had been exhausted.
        """

        mp.mp.dps = 40
        grid = probability_grid("0.035", 3)
        keys = [mp.nstr(p, 12) for p in grid]

        def scan(profile):
            return {"subsets": {"full": {
                key: {
                    "channels": {name: {"leverage": mp.mpf(profile[name][i])} for name in BASIS_CHANNELS},
                    "combination": {"leverage": mp.mpf(profile["combination"][i])},
                }
                for i, key in enumerate(keys)
            }}}

        interior = {name: [1, 5, 1] for name in list(BASIS_CHANNELS) + ["combination"]}
        self.assertEqual(
            set(grid_boundary_binds(scan(interior), grid).values()), {False}
        )
        edge = dict(interior)
        edge["Sp"] = [9, 5, 1]
        edge["Dp"] = [1, 5, 9]
        binds = grid_boundary_binds(scan(edge), grid)
        self.assertTrue(binds["Sp"])
        self.assertTrue(binds["Dp"])
        self.assertFalse(binds["M"])


class Eigen(unittest.TestCase):
    def test_jacobi_reconstructs_the_matrix(self) -> None:
        """Stops us believing a combination leverage built on a wrong spectrum.

        ``best_combination`` divides by eigenvalues.  A wrong decomposition would
        not fail loudly -- it would return a plausible number, and that number is
        the headline of this scan.  Checked by independent means: reassemble
        ``sum lambda_k v_k v_k^T`` and compare to the input.
        """

        mp.mp.dps = 40
        matrix = [
            [mp.mpf("4"), mp.mpf("1"), mp.mpf("-2"), mp.mpf("0.5")],
            [mp.mpf("1"), mp.mpf("3"), mp.mpf("0.25"), mp.mpf("-1")],
            [mp.mpf("-2"), mp.mpf("0.25"), mp.mpf("6"), mp.mpf("0.75")],
            [mp.mpf("0.5"), mp.mpf("-1"), mp.mpf("0.75"), mp.mpf("2")],
        ]
        values, vectors = symmetric_eigen(matrix)
        for i in range(4):
            for j in range(4):
                rebuilt = mp.fsum(
                    values[k] * vectors[k][i] * vectors[k][j] for k in range(4)
                )
                self.assertLess(abs(rebuilt - matrix[i][j]), mp.mpf("1e-25"))
        for k in range(4):
            norm = mp.fsum(component**2 for component in vectors[k])
            self.assertLess(abs(norm - 1), mp.mpf("1e-25"))


class Combination(unittest.TestCase):
    def test_combination_matches_the_closed_form_on_a_diagonal_covariance(self) -> None:
        """Stops us believing a Hotelling leverage that is not v^T Sigma^-1 v.

        On a diagonal covariance the optimum is the root-sum-square of the
        per-coordinate leverages.  A wrong normalization here would misprice the
        sample saving by its square, which is the number a production decision
        would be made on.
        """

        mp.mp.dps = 40
        vector = [mp.mpf("3"), mp.mpf("-4"), mp.mpf("1"), mp.mpf("0")]
        covariance = [
            [mp.mpf("1") if i == j else mp.mpf("0") for j in range(4)] for i in range(4)
        ]
        covariance[1][1] = mp.mpf("4")
        got = best_combination(vector, covariance)
        expected = mp.sqrt(mp.mpf(9) + mp.mpf(16) / 4 + 1)
        self.assertLess(abs(got["leverage"] - expected), mp.mpf("1e-25"))
        self.assertEqual(got["retained"], 4)

    def test_combination_never_beats_its_best_single_coordinate_downward(self) -> None:
        """Stops us believing a combination that is worse than a channel it contains.

        The optimum over a space includes every coordinate axis, so it can never
        resolve worse than the best single channel.  If it did, the sign or the
        inverse would be wrong somewhere and the direction reported to #583
        would be meaningless.
        """

        mp.mp.dps = 40
        vector = [mp.mpf("0.4"), mp.mpf("-1.1"), mp.mpf("0.05"), mp.mpf("2.0")]
        covariance = [
            [mp.mpf("2.0"), mp.mpf("0.6"), mp.mpf("-0.3"), mp.mpf("0.1")],
            [mp.mpf("0.6"), mp.mpf("1.5"), mp.mpf("0.2"), mp.mpf("-0.4")],
            [mp.mpf("-0.3"), mp.mpf("0.2"), mp.mpf("0.9"), mp.mpf("0.05")],
            [mp.mpf("0.1"), mp.mpf("-0.4"), mp.mpf("0.05"), mp.mpf("3.0")],
        ]
        combined = best_combination(vector, covariance)["leverage"]
        singles = [
            abs(vector[i]) / mp.sqrt(covariance[i][i]) for i in range(4)
        ]
        self.assertGreaterEqual(combined, max(singles) - mp.mpf("1e-25"))

    def test_unresolvable_directions_are_dropped_not_inverted(self) -> None:
        """Stops us believing leverage manufactured out of a singular direction.

        A covariance direction the sample cannot resolve has near-zero
        eigenvalue; inverting it turns numerical noise into apparent signal and
        would report an unbuyable sample saving.
        """

        mp.mp.dps = 40
        vector = [mp.mpf("1"), mp.mpf("1"), mp.mpf("1"), mp.mpf("1")]
        covariance = [
            [mp.mpf("1") if i == j else mp.mpf("0") for j in range(4)] for i in range(4)
        ]
        covariance[3][3] = COMBINATION_FLOOR / 100
        got = best_combination(vector, covariance)
        self.assertEqual(got["retained"], 3)
        self.assertLess(abs(got["leverage"] - mp.sqrt(mp.mpf(3))), mp.mpf("1e-25"))


class SelectionDiscipline(unittest.TestCase):
    def _scan(self, profile, keys):
        return {"subsets": {
            subset: {
                key: {
                    "channels": {
                        name: {"leverage": mp.mpf(profile[subset][name][i])}
                        for name in BASIS_CHANNELS
                    },
                    "combination": {"leverage": mp.mpf(profile[subset]["combination"][i])},
                }
                for i, key in enumerate(keys)
            }
            for subset in profile
        }}

    def test_a_readout_that_only_wins_in_sample_does_not_transport(self) -> None:
        """Stops us believing an amplification that is selection optimism.

        This scan maximizes over 105 readouts against a finite sample, so it
        *will* find a winner even from pure noise.  A control that cannot report
        failure is decoration; this plants a readout that wins on one half and
        loses on the other, and requires the control to say so.
        """

        mp.mp.dps = 40
        grid = probability_grid("0.035", 3)
        keys = [mp.nstr(p, 12) for p in grid]
        flat = {name: [1, 1, 1] for name in list(BASIS_CHANNELS) + ["combination"]}
        lucky = {**{name: [1, 1, 1] for name in list(BASIS_CHANNELS) + ["combination"]},
                 "Sp": [1, 9, 1]}
        unlucky = {**{name: [1, 1, 1] for name in list(BASIS_CHANNELS) + ["combination"]},
                   "Sp": [1, 0.5, 1]}
        scan = self._scan({"full": flat, "half_0": lucky, "half_1": unlucky}, keys)
        control = selection_control(scan, grid)
        self.assertFalse(control["both_directions_transport"])
        forward = control["directions"][0]
        self.assertEqual(forward["readout"], "Sp")
        self.assertLess(mp.mpf(forward["held_out_amplification"]), 1)
        self.assertGreater(mp.mpf(forward["optimism_factor"]), 1)

    def test_a_post_hoc_grid_is_refused_a_verdict_by_code(self) -> None:
        """Stops us believing a design rule justified by a grid we widened after looking.

        The predeclared grid's optimum sits on its boundary, so the temptation to
        widen and re-read is real and immediate.  The refusal has to be in the
        code, not in the discipline of whoever writes the note.
        """

        controls = {
            n: {"both_directions_transport": True,
                "directions": [{"held_out_amplification": "5.0"}]}
            for n in (325, 425)
        }
        self.assertIn("NO_VERDICT", decide(controls, predeclared=False)["verdict"])
        self.assertNotIn("NO_VERDICT", decide(controls, predeclared=True)["verdict"])
        self.assertIn(
            "worst_held_out_amplification", decide(controls, predeclared=True)
        )
        self.assertNotIn(
            "worst_held_out_amplification", decide(controls, predeclared=False)
        )

    def test_a_spread_that_does_not_reproduce_is_called_noise(self) -> None:
        """Stops us believing the 9000x channel spread when it does not survive a split.

        If the held-out amplification never beats the baseline, the honest answer
        is that the observed spread is sampling noise -- including Sp's advantage
        over M, which #589's cost table already leans on.
        """

        controls = {
            n: {"both_directions_transport": False,
                "directions": [{"held_out_amplification": "0.8"}]}
            for n in (325, 425)
        }
        self.assertIn("NOT_REPRODUCIBLE", decide(controls, predeclared=True)["verdict"])


class DeclaredManifest(unittest.TestCase):
    def test_the_dependent_channel_is_excluded_and_the_baseline_is_included(self) -> None:
        """Stops us believing a combination covariance that is singular by construction.

        ``obs`` defines ``D = M/2``, so keeping both would make the 4x4 amplitude
        covariance exactly rank-deficient and hand ``best_combination`` a
        direction with zero eigenvalue in every sample.
        """

        self.assertNotIn("D", BASIS_CHANNELS)
        self.assertIn(BASELINE_CHANNEL, BASIS_CHANNELS)
        self.assertEqual(len(set(BASIS_CHANNELS)), len(BASIS_CHANNELS))
        self.assertEqual(P_GRID_POINTS % 2, 1)
        self.assertEqual(SPLIT_PARITY, 2)
        self.assertGreater(mp.mpf(P_HALF_WIDTH), 0)

    def test_readouts_enumerates_every_candidate_the_selection_sees(self) -> None:
        """Stops us believing an optimism estimate computed over too few candidates.

        The held-out control is only honest if it selects from the same pool the
        headline number was chosen from.  A readout visible to the cost table but
        invisible to the selection would understate the optimism.
        """

        mp.mp.dps = 40
        grid = probability_grid("0.035", 3)
        keys = [mp.nstr(p, 12) for p in grid]
        profile = {name: [1, 2, 3] for name in list(BASIS_CHANNELS) + ["combination"]}
        scan = {"subsets": {"full": {
            key: {
                "channels": {name: {"leverage": mp.mpf(profile[name][i])} for name in BASIS_CHANNELS},
                "combination": {"leverage": mp.mpf(profile["combination"][i])},
            }
            for i, key in enumerate(keys)
        }}}
        flat = readouts(scan, "full")
        self.assertEqual(len(flat), (len(BASIS_CHANNELS) + 1) * len(keys))
        self.assertIn(("combination", keys[0]), flat)


if __name__ == "__main__":
    unittest.main()
