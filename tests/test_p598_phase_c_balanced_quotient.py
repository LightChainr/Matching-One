"""Tests for the #600 Phase C balanced-quotient comparison.

The claim under test is a *contingent* theorem, not an identity: the protected
dictionary must factor through the R-orbit quotient at every width, the R-odd
sector must carry no Hankel energy either way, and the exposed dictionary must
*break* at the two odd widths (where ``halves_linked`` is not R-even).  Each
test names the wrong number it would stop us believing.
"""

from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p598_phase_c_balanced_quotient import (  # noqa: E402
    CONTINUITY_ARTIFACT,
    decide,
    width_report,
)

WIDTHS = (4, 5, 6, 7, 8)
# Catalan(w) microscopic states and r_reflect(w) orbits, declared in #598.
STATES = {4: 14, 5: 42, 6: 132, 7: 429, 8: 1430}
ORBITS = {4: 10, 5: 26, 6: 76, 7: 232, 8: 750}
ODD_DIM = {w: STATES[w] - ORBITS[w] for w in WIDTHS}


def _report(width: int):
    return width_report(width)


class QuotientIsWellDefined(unittest.TestCase):
    def test_orbit_counts_match_the_reflection_partition(self) -> None:
        """Stops us believing a quotient built on the wrong partition of states."""
        for width in WIDTHS:
            report = _report(width)
            self.assertEqual(report["states"], STATES[width])
            self.assertEqual(report["quotient_states"], ORBITS[width])

    def test_odd_sector_dimension_is_the_complement(self) -> None:
        """Stops us believing the odd sector has a dimension it cannot have."""
        for width in WIDTHS:
            report = _report(width)
            self.assertEqual(report["odd_sector"]["dimension"], ODD_DIM[width])

    def test_quotient_rows_do_not_depend_on_representative(self) -> None:
        """Stops us believing a quotient that changes with the chosen representative."""
        for width in WIDTHS:
            report = _report(width)
            self.assertEqual(
                report["quotient_is_well_defined"]["max_representative_dependence"],
                0.0,
            )

    def test_all_sources_are_reflection_even(self) -> None:
        """Stops us believing the source lifting is nontrivial when it is not."""
        for width in WIDTHS:
            report = _report(width)
            self.assertTrue(report["parity"]["all_sources_are_R_even"])


class ProtectedFactorsThroughQuotient(unittest.TestCase):
    def test_hankel_spectra_agree_on_the_protected_dictionary(self) -> None:
        """Stops us believing the protected task has a microscopic-only direction."""
        for width in WIDTHS:
            report = _report(width)
            diff = report["dictionaries"]["protected"]["readings"]["restriction"][
                "differences"
            ]
            self.assertLessEqual(
                diff["hankel_spectrum"]["max_absolute"], 1e-12, msg=f"w={width}"
            )

    def test_numerical_rank_is_preserved(self) -> None:
        """Stops us believing reduction drops or creates a resolvable direction."""
        for width in WIDTHS:
            report = _report(width)
            rank = report["dictionaries"]["protected"]["readings"]["restriction"][
                "differences"
            ]["numerical_rank"]
            self.assertEqual(rank["A"], rank["B"], msg=f"w={width}")


class OddSectorIsInert(unittest.TestCase):
    def test_reach_energy_has_no_odd_component(self) -> None:
        """Stops us believing the R-odd states are reachable from R-even sources."""
        for width in WIDTHS:
            report = _report(width)
            self.assertLessEqual(
                report["odd_sector"]["reach_odd_energy_relative"], 1e-14, msg=f"w={width}"
            )

    def test_observe_energy_has_no_odd_component(self) -> None:
        """Stops us believing the R-odd states are observable in R-even readouts."""
        for width in WIDTHS:
            report = _report(width)
            self.assertLessEqual(
                report["odd_sector"]["observe_odd_energy_relative"],
                1e-14,
                msg=f"w={width}",
            )


class ExposedDictionaryHasTeeth(unittest.TestCase):
    def test_restriction_disagrees_at_odd_widths(self) -> None:
        """Stops us believing a comparison with no control -- odd-width
        ``halves_linked`` is not R-even, so restriction must disagree there."""
        for width in (5, 7):
            report = _report(width)
            diff = report["dictionaries"]["exposed"]["readings"]["restriction"][
                "differences"
            ]
            self.assertGreater(
                diff["hankel_spectrum"]["max_absolute"], 1e-6, msg=f"w={width}"
            )

    def test_restriction_agrees_at_even_widths(self) -> None:
        """Stops us believing the odd-width break is a generic numerical failure
        rather than the R-parity of ``halves_linked`` specifically."""
        for width in (4, 6, 8):
            report = _report(width)
            diff = report["dictionaries"]["exposed"]["readings"]["restriction"][
                "differences"
            ]
            self.assertLessEqual(
                diff["hankel_spectrum"]["max_absolute"], 1e-12, msg=f"w={width}"
            )

    def test_symmetrization_restores_agreement_at_odd_widths(self) -> None:
        """Stops us believing the odd-width disagreement is anything but the odd
        part: symmetrizing the readout removes exactly that part and the chains
        agree again."""
        for width in (5, 7):
            report = _report(width)
            sym = report["dictionaries"]["exposed"]["readings"]["symmetrization"]
            self.assertFalse(sym["identical_to_restriction"], msg=f"w={width}")
            self.assertLessEqual(
                sym["differences"]["hankel_spectrum"]["max_absolute"],
                1e-12,
                msg=f"w={width}",
            )


class Decision(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls._blocks = [width_report(width) for width in WIDTHS]
        cls._decision = decide(cls._blocks)

    def test_three_claims_are_all_true(self) -> None:
        """Stops us believing any one leg of the theorem failed."""
        self.assertTrue(self._decision["protected_factors_through_the_quotient"])
        self.assertTrue(self._decision["odd_sector_is_inert"])
        self.assertTrue(self._decision["exposed_disagrees_at_odd_widths"])
        self.assertTrue(self._decision["exposed_agrees_at_even_widths"])


@unittest.skipUnless(CONTINUITY_ARTIFACT.exists(), "committed #588 artifact absent")
class ContinuityGate(unittest.TestCase):
    def test_port_reproduces_committed_balanced_scores(self) -> None:
        """Stops us trusting a numpy port that silently changed the numbers.

        The committed artifact stores ten significant digits, so the gate must
        hold to that resolution rather than to machine precision.
        """
        from p598_phase_c_balanced_quotient import continuity_gate

        for width in WIDTHS:
            gate = continuity_gate(width)
            if gate.get("status") != "run":
                self.fail(f"w={width}: gate status {gate.get('status')}")
            self.assertLessEqual(
                gate["hankel_spectrum_above_tolerance"]["max_absolute_difference"],
                1e-7,
                msg=f"w={width} spectrum above tolerance",
            )
            self.assertLessEqual(
                gate["max_score_difference"], 1e-7, msg=f"w={width} score"
            )


class ArtifactShape(unittest.TestCase):
    def test_produced_artifact_has_the_declared_schema(self) -> None:
        """Stops us shipping an artifact without its provenance shape."""
        from p598_phase_c_balanced_quotient import DEFAULT_OUTPUT

        if not DEFAULT_OUTPUT.exists():
            self.skipTest("artifact not yet generated")
        payload = json.loads(DEFAULT_OUTPUT.read_text())
        self.assertEqual(payload["schema"], "matching-one.p598-balanced-quotient.v1")
        self.assertEqual(payload["issue"], 600)
        self.assertEqual([b["width"] for b in payload["blocks"]], list(WIDTHS))


if __name__ == "__main__":
    unittest.main()
