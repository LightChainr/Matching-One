import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from p333_source_landing_doublet_width4 import landing_gram_jet
from p398_rooted_gr1_completion import selected_completion_families
from p398_width8_source_spectrum import analyze, model


class WidthEightSourceSpectrum(unittest.TestCase):
    def test_width_four_same_readouts(self):
        _,_,functions,*_ = model(4)
        old_a = selected_completion_families()["rooted_charge1"]["columns"]
        old_l,_ = landing_gram_jet()
        old = np.array([[complex(row[0],row[1]),complex(old_l[j][-2],old_l[j][-1])]
                        for j,row in enumerate(old_a)])
        np.testing.assert_array_equal(functions,old)

    def test_width_eight_exact_dimension_and_duality(self):
        value = analyze(8)
        self.assertEqual(value["states"],1430)
        self.assertEqual(value["exact_character_sector_dimension"],186)
        self.assertEqual(value["source_krylov"]["reachable_rank_lower_bound"],186)
        for sign in ("minus","plus"):
            self.assertEqual(value["kreweras"]["source_"+sign+"_krylov"]["reachable_rank_lower_bound"],93)
        self.assertLess(value["max_cross_ray_spectral_residue"],1e-9)
        self.assertLess(value["max_eigenpair_residual"],1e-10)


if __name__=="__main__":
    unittest.main()
