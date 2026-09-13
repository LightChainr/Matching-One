import sys
from pathlib import Path
import unittest
from fractions import Fraction as F
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import p398_tl_fattening as t

class TestFattening(unittest.TestCase):
    def test_bijection_independent_rgs(self):
        for w in range(2,7):
            states=t.state_space(w)
            self.assertEqual(states,t.independent_rgs(w))
            self.assertEqual({t.fatten(s) for s in states},set(t.link_patterns(2*w)))
    def test_each_move(self):
        for s in t.state_space(5):
            for i in range(5):
                self.assertEqual(t.fatten(t.join(s,i)),t.tl_reconnect(t.fatten(s),2*i+1))
                self.assertEqual(t.fatten(t.detach(s,i)),t.tl_reconnect(t.fatten(s),2*i))
    def test_periodic_seam(self):
        s=(0,1,2,3)
        self.assertEqual(t.fatten(t.join(s,3)),t.tl_reconnect(t.fatten(s),7))
    def test_half_step_not_involution(self):
        s=(0,0,1,2,3)
        self.assertEqual(t.half_step(t.half_step(s)),t.canonical((s[-1],)+s[:-1]))
        self.assertNotEqual(t.half_step(t.half_step(s)),s)
    def test_readout_transport(self):
        for s in t.state_space(6):
            k=t.half_step(s)
            self.assertEqual(max(s)+max(k)+2,7)
            self.assertEqual(s[0]==s[-1],k.count(k[0])==1)
    def test_staggered_generator_conjugacy(self):
        ss=t.state_space(4);ix={s:i for i,s in enumerate(ss)};k=[ix[t.half_step(s)] for s in ss]
        gp=t.generator(ss,F(1,4));gm=t.generator(ss,F(-1,4))
        self.assertTrue(all(gp[i][j]==gm[k[i]][k[j]] for i in range(len(ss)) for j in range(len(ss))))
    def test_stationary_normalization(self):
        ss=t.state_space(4);pi=t.stationary(t.generator(ss))
        self.assertEqual(sum(p*(max(s)+1) for p,s in zip(pi,ss)),F(5,2))
        self.assertEqual(t.asm_number(4),42)
        self.assertEqual(sorted(set(p*42 for p in pi)),[1,3,7])
    def test_input_guards(self):
        with self.assertRaises(ValueError):t.link_patterns(5)
        with self.assertRaises(ValueError):t.unfatten((0,1))
        with self.assertRaises(ValueError):t.generator(t.state_space(3),F(2))

if __name__=='__main__':unittest.main()
