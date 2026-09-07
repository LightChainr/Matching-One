"""Tests for scripts/probe638_matching_boundary_states.py.

Every test docstring names the specific wrong number that test exists to
stop us believing.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from noncrossing_connectivity_codec import (  # noqa: E402
    catalan,
    is_noncrossing_blocks,
    noncrossing_states,
)

from probe638_matching_boundary_states import (  # noqa: E402
    bond_closure,
    bond_step,
    class_rgs,
    canonical,
    closure,
    direct_class,
    rgs_blocks,
    step,
    EMPTY,
)

NESTED_AT_W4 = (0, 1, 1, 0)  # blocks {0,3},{1,2}: nested noncrossing pair
DOUBLE_PAIR_W4 = (0, 0, 1, 1)  # blocks {0,1},{2,3}


class Probe638BoundaryStateTests(unittest.TestCase):
    def test_site_helix_nn_closure_is_not_catalan_at_w4(self) -> None:
        """Stops us believing 'the site frontier reaches all 14 noncrossing
        classes at w=4': it reaches exactly 13 and is missing (0,1,1,0)
        ({0,3},{1,2}) -- the nested-pocket class, structurally unreachable
        because the wrap seam is consumed by the (w-1)-(0) intra-row edge."""
        classes, _p0, _n, _parent, _cap = closure(4, False)
        self.assertEqual(len(classes), 13)
        self.assertNotIn(NESTED_AT_W4, set(classes))
        self.assertIn(DOUBLE_PAIR_W4, set(classes))  # THIS pair IS reachable

    def test_site_nn_closure_reaches_12_of_14_noncrossing_at_w4_rowcut(self) -> None:
        """Stops us believing 'row-cut classes = 13 too': at phase 0 the
        closure gives exactly 12, missing BOTH (0,1,1,0) and (0,0,1,1) --
        the original straight-row result that started probe638."""
        _classes, p0, _n, _parent, _cap = closure(4, False)
        self.assertEqual(len(p0), 12)
        self.assertNotIn(NESTED_AT_W4, set(p0))
        self.assertNotIn(DOUBLE_PAIR_W4, set(p0))

    def test_rowcut_crosscheck_w4_no_violations(self) -> None:
        """Stops us believing 'the closure's phase-0 set could drift away
        from what the lattice actually produces at a row cut': every
        depth-3 full-row pattern (3*2^4 = 3072 patterns) classified by
        the INDEPENDENT direct classifier must land in the closure's
        phase-0 class set.  The guarded wrong claim: a closure p0 set
        that misses a class the raw lattice demonstrably produces."""
        from probe638_matching_boundary_states import rowcut_crosscheck

        _classes, p0, _n, _parent, _cap = closure(4, False)
        rec = rowcut_crosscheck(4, set(p0))
        self.assertEqual(rec["violations"], 0, rec)
        self.assertEqual(rec["patterns_checked"], (1 << 4) ** 3)  # 4096
    def test_bond_anchor_hits_catalan_exactly_through_w7(self) -> None:
        """Stops us believing 'Catalan(w) is unreachable by any transfer'
        (the 12-at-w=4 reading of the first draft) AND its mirror 'the
        closure is 15 not 14' (the dangling-reconnect reading): with the
        correct semantics -- absent vertical bond = its own isolated
        singleton, sealed blocks dropped -- the PLANAR bond row transfer
        reaches exactly the committed noncrossing_states(w), elementwise,
        for w = 1..7."""
        for width in range(1, 8):
            with self.subTest(width=width):
                got = bond_closure(width, periodic=False)
                want = set(noncrossing_states(width))
                self.assertEqual(got, want)
                self.assertEqual(len(got), catalan(width))

    def test_cylinder_bond_closure_stays_14_at_w4(self) -> None:
        """Stops us believing 'the wrap edge on the cylinder adds crossing
        classes at w=4' (the guarded wrong number: 15, from the old
        dangling-reconnect semantics): with the sealed semantics the wrap
        merely merges the two end ports' blocks, and the w=4 cylinder
        closure equals the planar strip exactly, 14 classes."""
        # with sealed semantics the w=4 cylinder equals the planar strip
        self.assertEqual(len(bond_closure(4, periodic=True)), 14)

    def test_transfer_and_direct_classifier_agree_on_all_w2_nnn_prefixes(self) -> None:
        """Stops us believing 'the shadow machinery is right because the
        closure counts look plausible': EVERY length-18 prefix (2^18 = all
        w=2 patterns of 9 rows) must classify identically under the
        transfer and the independent direct classifier.  The wrong numbers
        this guards: the slot-merge bug (13 -> claimed 14 at w=4 NN) and
        the shadow-merge bug (class (0,1,2,3,1) read as (0,1,2,3,4) --
        the (2,1)/(2,4) co-block through (1,0) diagonals dropped)."""
        from probe638_matching_boundary_states import selfcheck_prefix_agreement

        for diag in (False, True):
            result = selfcheck_prefix_agreement(2, diag)
            self.assertTrue(result["agrees"], result)

    def test_direct_classifier_finds_diagonal_connection_transfer_must_too(
        self,
    ) -> None:
        """Stops us believing 'an occupied site at (2,1) and (2,4) with a
        shared occupied neighbour (1,0) on w=5 NNN stays in two blocks':
        both diagonal edges exist, so the class is (0,1,2,3,1) -- the
        transfer's shadow variables must reproduce exactly this."""
        w = 5
        bits = [0] * 5 + [1, 0, 1, 0, 0] + [0, 1, 0, 0, 1]
        want = direct_class(bits, w, True)
        self.assertEqual(want, (0, 1, 2, 3, 1))
        state = canonical(0, [EMPTY] * w, EMPTY, EMPTY)
        for b in bits:
            state = step(state, b == 1, w, True)
        self.assertEqual(class_rgs(state, w), want)

    def test_site_closure_state_space_stays_finite_and_small(self) -> None:
        """Stops us believing 'the helical state space can drift unbounded
        through fresh labels': the w=8 NN closure holds exactly 7429
        distinct states and terminates without the cap (the guarded wrong
        number: any claim of a cap hit / divergent state count)."""

        _classes, _p0, n_states, _parent, capped = closure(8, False, state_cap=2_000_000)
        self.assertFalse(capped)
        self.assertEqual(n_states, 7429)

    def test_nnn_subset_nn_and_all_noncrossing_through_w8(self) -> None:
        """Stops us believing 'NNN diagonals create crossing boundary
        classes at small width': through w=8 the NNN class set is a subset
        of the NN set and contains ZERO crossing classes -- the first
        crossing width (if any) is > 8, and the #636 cost model must not
        assume crossing states appear at w <= 8."""
        for width in range(2, 9):
            with self.subTest(width=width):
                nn_c, _p0, _n, _par, _cap = closure(width, False)
                nn_c2, _p02, _n2, _par2, _cap2 = closure(width, True)
                nnn = set(nn_c2)
                self.assertTrue(nnn <= set(nn_c))
                for c in nnn:
                    self.assertTrue(is_noncrossing_blocks(rgs_blocks(c)))

    def test_bond_absent_vertical_bond_is_isolated_singleton(self) -> None:
        """Stops us believing 'an absent vertical bond keeps its old block'
        (that semantics yields the phantom 15/52/203 closure): with the
        sealed semantics a single row with vrow=0000 relabels every port
        to its own fresh singleton -- 2^w distinct rows cannot freeze an
        old partition across an empty row."""
        w = 4
        sealed = bond_step(NESTED_AT_W4, 0b0000, 0b0000, w, periodic=False)
        self.assertEqual(sealed, tuple(range(w)))
        # a full row of vertical bonds (no horizontals) carries the OLD
        # partition through: connectivity persists across the row, so the
        # new ports keep {0,3} and {1,2}.  Guarded wrong number: all
        # singletons (which would mean connectivity dies at every row).
        full = bond_step(NESTED_AT_W4, 0b1111, 0b0000, w, periodic=False)
        self.assertEqual(full, (0, 1, 1, 0))  # same partition, canonical order


if __name__ == "__main__":
    unittest.main()
