"""The missing off-diagonal is recoverable from aligned deletions."""
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from recover_n580_covariance import jackknife_covariance

class RecoveryTests(unittest.TestCase):
    def test_covariance_from_deletions(self):
        # n=3, centered columns (-1,0,1),(-2,0,2),(1,0,-1).
        got=jackknife_covariance([[1,2,3],[2,4,6],[3,2,1]])
        self.assertAlmostEqual(got[0][0],4/3)
        self.assertAlmostEqual(got[1][2],-8/3)
        self.assertEqual(got[1][2],got[2][1])

    def test_joint_permutation_preserves_covariance(self):
        a=[[1,2,5],[3,-1,7],[0,2,-3]]
        b=[[v[i] for i in (2,0,1)] for v in a]
        self.assertEqual(jackknife_covariance(a),jackknife_covariance(b))

if __name__=='__main__':
    unittest.main()
