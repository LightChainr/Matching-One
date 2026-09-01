from __future__ import annotations
import unittest
from fractions import Fraction as F
import landing_minor as lm

class LandingMinorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.l4=lm.l4()
    def test_existing_oracle(self):
        self.assertEqual(self.l4['cross'],self.l4['expected']);self.assertEqual(self.l4['shift_violations'],0)
    def test_l4_matrix(self):
        m=self.l4['matrix'];self.assertEqual(F(m['T_01']),F(-527,16384));self.assertEqual(F(m['T_12']),F(337,32768));self.assertEqual(F(m['A_01']),F(16710343,68719476736));self.assertEqual(F(m['A_12']),F(5521655,137438953472));self.assertEqual(F(m['det']),F(-533831111,140737488355328))
    def test_positive_control(self): self.assertEqual(self.l4['positive_control_a_equals_K_det'],'0')
    def test_transition_counts(self): self.assertEqual(self.l4['state_counts'],{'0_to_1':2048,'1_to_2':1190})
    def test_all_scale_family(self):
        for L,R in ((7,1),(9,2),(11,3),(13,4)):
            z=lm.family(L,R);A,B=z['states'];self.assertEqual(A['transition'],[0,1]);self.assertEqual(B['transition'],[1,2]);self.assertEqual(A['h4'],1);self.assertEqual(B['h4'],1);self.assertEqual(F(B['amid'])-F(A['amid']),F(-2,L**4));self.assertNotEqual(F(z['minor']),0)
    def test_source_formulas(self):
        for L,R in ((7,1),(10,2),(13,4)):
            A,B=lm.family(L,R)['states'];self.assertEqual(F(A['a0']),F(2*(L-2),L**4));self.assertEqual(F(A['a1']),F(2*(L-1),L**4));self.assertEqual(F(B['a0']),F(2*(L-3),L**4));self.assertEqual(F(B['a1']),F(2*(L-2),L**4))

if __name__=='__main__':unittest.main()
