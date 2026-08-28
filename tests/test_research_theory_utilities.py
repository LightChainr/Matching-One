from __future__ import annotations

import sys
import unittest
from fractions import Fraction
from pathlib import Path

SCRIPTS=Path(__file__).resolve().parents[1]/"scripts"
sys.path.insert(0,str(SCRIPTS))

from c4_self_matching_cyclic_geometry import build as build_c4, validate as validate_c4  # noqa:E402
from gaussian_harmonic_arithmetic import harmonic, gmul, norm  # noqa:E402
from thermal_spin4_tower import kac_singular_levels, quasiprimary_count  # noqa:E402
from generic_potts_singlet_spin4_gap import row as spin4_primary_row  # noqa:E402
from v14_scalar_post_l7 import diagonal_x, interchiral_parity  # noqa:E402


class C4SelfMatchingGeometryTests(unittest.TestCase):
    def test_small_and_production_geometries(self) -> None:
        for a,b in ((3,1),(11,3),(9,7),(13,1),(11,7),(17,1),(13,11)):
            g=build_c4(a,b)
            info=validate_c4(g)
            self.assertEqual(info["V"],a*a+b*b)
            self.assertEqual(info["E"],3*info["V"])
            self.assertEqual(info["F"],2*info["V"])
            self.assertEqual(info["even_degree"],8)
            self.assertEqual(info["odd_degree"],4)
            self.assertEqual((info["rotation_multiplier"]**2+1)%info["V"],0)

    def test_rejects_nonprimitive_or_even_period(self) -> None:
        with self.assertRaises(ValueError): build_c4(5,5)
        with self.assertRaises(ValueError): build_c4(4,1)


class GaussianHarmonicTests(unittest.TestCase):
    def test_exact_harmonic_identities(self) -> None:
        for z in ((8,1),(7,4),(17,4),(16,7),(18,1),(17,6)):
            c4,s4=harmonic(z,1); c8,s8=harmonic(z,2); c12,s12=harmonic(z,3)
            self.assertEqual(c4*c4+s4*s4,1)
            self.assertEqual(c8,2*c4*c4-1)
            self.assertEqual(s8,2*s4*c4)
            self.assertEqual(c12,4*c4*c4*c4-3*c4)
            self.assertEqual(s12,s4*(4*c4*c4-1))

    def test_one_plus_i_spin4_sign(self) -> None:
        h=(1,1)
        for z in ((8,1),(7,4),(9,2),(7,6),(12,1),(9,8)):
            child=gmul(z,h)
            self.assertEqual(norm(child),2*norm(z))
            c4,_=harmonic(z,1); cc4,_=harmonic(child,1)
            self.assertEqual(cc4,-c4)
            c8,_=harmonic(z,2); cc8,_=harmonic(child,2)
            self.assertEqual(cc8,c8)
            c12,_=harmonic(z,3); cc12,_=harmonic(child,3)
            self.assertEqual(cc12,-c12)

    def test_frozen_h4_h12_design_rationals(self) -> None:
        dc4_305=harmonic((17,4),1)[0]-harmonic((16,7),1)[0]
        dc12_305=harmonic((17,4),3)[0]-harmonic((16,7),3)[0]
        self.assertEqual(dc4_305,Fraction(12672,18605))
        self.assertEqual(dc12_305/dc4_305,Fraction(-14829638967,8653650625))
        dc4_325=harmonic((18,1),1)[0]-harmonic((17,6),1)[0]
        dc12_325=harmonic((18,1),3)[0]-harmonic((17,6),3)[0]
        self.assertEqual(dc4_325,Fraction(16128,21125))
        self.assertEqual(dc12_325/dc4_325,Fraction(1555994781,858203125))


class ThermalTowerTests(unittest.TestCase):
    def test_singular_levels_and_spin4_gap(self) -> None:
        self.assertEqual(kac_singular_levels()[:3],[2,10,16])
        counts=[quasiprimary_count(i) for i in range(10)]
        self.assertEqual(counts,[1,0,0,1,1,1,2,2,3,4])
        allowed=[]
        for m in range(6):
            n=m+4
            if counts[m] and counts[n]: allowed.append((n,m,n+m))
        self.assertEqual(allowed[:2],[(4,0,4),(7,3,10)])


class CorrectedCriticalPottsBranchTests(unittest.TestCase):
    def test_non_diagonal_spin4_singlet_gap(self) -> None:
        self.assertEqual(spin4_primary_row(2)["x"], Fraction(17, 4))
        self.assertEqual(spin4_primary_row(3)["x"], Fraction(17, 4))
        self.assertEqual(spin4_primary_row(4)["x"], Fraction(6))

    def test_v14_scalar_post_l7_arithmetic(self) -> None:
        self.assertEqual(diagonal_x(4), Fraction(33, 4))
        self.assertEqual(interchiral_parity(4), -1)
        self.assertEqual(diagonal_x(4) - Fraction(21, 4), 3)


if __name__=="__main__":
    unittest.main()
