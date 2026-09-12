"""Checks of the geometric equality classes, not substitutes for the proof."""
from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import rank_two_onset as o

class TwoCycleCoreTests(unittest.TestCase):
    def test_existing_small_diamond_exhaustive_classification(self):
        for L in (2,3):
            result=o.check_small(L,True)
            self.assertEqual(result['rank2_below_onset'],0)
            self.assertEqual(result['minimizers'],4*L*L)

    def test_axis_equality_is_a_cross(self):
        for L in (2,3,4):
            result=o.check_small(L,False)
            self.assertEqual(result['rank2_below_onset'],0)
            self.assertEqual(result['minimizers'],L*L)

    def test_all_constructed_larger_minimizers_are_rank_two(self):
        for L in range(2,9):
            for diamond in (False,True):
                sets=o.predicted_minimizers(L,diamond);adj=o.geometry(L,diamond)[3]
                self.assertEqual(len(sets),4*L*L if diamond else L*L)
                for mask in sets:
                    self.assertEqual(mask.bit_count(),3*L-1 if diamond else 2*L-1)
                    self.assertEqual(o.rank(mask,adj),2)

    def test_diamond_periods_are_the_named_geometry(self):
        for L in (2,3,5):
            _,_,canonical,_=o.geometry(L,True)
            for x,y in ((0,0),(1,2),(-3,7)):
                self.assertEqual(canonical(x,y),canonical(x+L,y+L))
                self.assertEqual(canonical(x,y),canonical(x+L,y-L))

if __name__=='__main__':unittest.main()
