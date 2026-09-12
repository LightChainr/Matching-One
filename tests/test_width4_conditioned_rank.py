import os
from pathlib import Path
import sys
import unittest
from fractions import Fraction
from random import Random
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
import width4_site_sources as ws
import width4_conditioned_rank as br
from site_source_graph_checks import lifted_graph
CERTIFICATE=Path(os.environ.get('MATCHING_ONE_RANK_CERTIFICATE',str(ws.DEFAULT_CERTIFICATE)))

class RankBridgeTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.source=ws.load_certificate(CERTIFICATE)

    def test_all_two_row_conditional_path_probabilities(self):
        ps=[[Fraction(1,3),Fraction(2,5),Fraction(1,2),Fraction(3,4)],
            [Fraction(2,3),Fraction(1,5),Fraction(3,7),Fraction(4,5)]]
        for rank in range(3):
            bridge=br.RankBridge(self.source,ps,rank)
            total=Fraction(0)
            for a in range(16):
                for b in range(16):
                    value=bridge.sequential_word_probability([a,b])
                    self.assertEqual(value,bridge.direct_word_probability([a,b]))
                    if lifted_graph([a,b])[0]!=rank:self.assertEqual(value,0)
                    total+=value
            self.assertEqual(total,1)

    def test_conditional_moments_and_thermal_score(self):
        p=Fraction(2,5);ps=[[p]*4]*3
        moments=br.sector_moments(self.source,ps)
        direct=[Fraction(0)]*3;derivative=[Fraction(0)]*3
        for mask in range(1<<12):
            word=[(mask>>(4*y))&15 for y in range(3)]
            r=lifted_graph(word)[0];k=mask.bit_count()
            weight=p**k*(1-p)**(12-k)
            direct[r]+=weight
            derivative[r]+=weight*(Fraction(k)-12*p)/(p*(1-p))
        self.assertEqual(direct,[m['probability'] for m in moments])
        lhs=derivative[2]/direct[2]-derivative[0]/direct[0]
        rhs=(moments[2]['mean']-moments[0]['mean'])/(p*(1-p))
        self.assertEqual(lhs,rhs)
        self.assertGreater(rhs,0)

    def test_samples_end_in_specified_physical_sector(self):
        for rank in (0,1,2):
            bridge=br.RankBridge(self.source,[[Fraction(3,5)]*4]*6,rank)
            for seed in (4,5,6):
                word=bridge.sample(Random(seed))
                self.assertEqual(lifted_graph(word)[0],rank)
                self.assertEqual(bridge.sequential_word_probability(word),bridge.direct_word_probability(word))

    def test_impossible_condition_and_invalid_word(self):
        with self.assertRaises(ValueError):br.RankBridge(self.source,[[Fraction(0)]*4]*2,2)
        bridge=br.RankBridge(self.source,[[Fraction(1,2)]*4]*2,0)
        with self.assertRaises(ValueError):bridge.direct_word_probability([0])
        with self.assertRaises(ValueError):br.integer_choice([0,0],Random(1))

if __name__ == '__main__':unittest.main()
