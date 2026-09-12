import sys
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from oblique_winding_corridor import (Torus, corridor, staircase, support,
    squared_distance_to_segment, rectangle_crossing, rank, circle_packing, exhaustive_case)


class ObliqueCorridorTests(unittest.TestCase):
    def test_staircase_endpoint_and_step_bound(self):
        for u in [(17,5),(16,8),(3,0),(-17,6),(0,-9),(-11,-13)]:
            for s in [1,2,3,5]:
                path=staircase(u,s)
                self.assertEqual(path[0],(0,0))
                self.assertEqual(path[-1],u)
                for a,b in zip(path,path[1:]):
                    self.assertLessEqual(max(abs(a[0]-b[0]),abs(a[1]-b[1])),s)
                for z in support(corridor(u,s)):
                    self.assertLessEqual(squared_distance_to_segment(z,u),16*s*s)

    def test_full_oblique_rings_close(self):
        for u,v in [((4,1),(0,4)),((4,2),(0,4)),((-5,2),(-3,-7))]:
            t=Torus(u,v); ev=corridor(u,1)
            occ={t.reduce(z) for z in support(ev)}
            self.assertTrue(all(rectangle_crossing(t,r,occ) for r in ev))
            self.assertGreater(rank(t,occ),0)

    def test_small_twisted_exhaustive(self):
        result=exhaustive_case((4,0),(1,3))
        self.assertEqual(result["implication_failures"],[])
        self.assertGreater(result["ring_configurations"],0)
        self.assertLess(result["ring_configurations"],result["configurations"])

    def test_greedy_translation_packing(self):
        for n in range(5,26):
            for pattern in [{0},{0,1},{0,1,3},{0,2,4}]:
                centers,diff=circle_packing(n,pattern)
                packed=[{(x+c)%n for x in pattern} for c in centers]
                self.assertGreaterEqual(len(centers)*len(diff),n)
                for i,a in enumerate(packed):
                    for b in packed[i+1:]: self.assertFalse(a&b)

    def test_quotient_and_seam_robustness(self):
        t=Torus((4,2),(-1,5))
        self.assertEqual(len(t.vertices()),t.n)
        for z in [(-9,7),(0,0),(13,-8)]:
            for u in [t.u,t.v]:
                self.assertEqual(t.reduce(z),t.reduce((z[0]+u[0],z[1]+u[1])))
        self.assertEqual(rank(t,set()),0)
        self.assertEqual(rank(t,set(t.vertices())),2)


if __name__ == "__main__": unittest.main()
