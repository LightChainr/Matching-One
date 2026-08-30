# P250 adaptive x spatial-spectrum fresh production

The fresh N505 production used 200,000 replicas in 400 covariance batches on
HZ, with a counter stream independent of the 20k pilot.  It finished in 111.83
seconds; all typed primary/complement controls had zero failures.

The adaptive response has a spatial shape.  A frozen cluster-preserving
randomization of the conditional response over all nonzero CRT residues rejects
spatial exchangeability (`p=0.0007998`, 5,000 permutations).  The largest mode
is `k=1` on the plus child and `k=10` on the minus child.  This is a joint,
multiple-frequency-controlled result; it is not a post hoc single-mode z score.

The shape nevertheless remains inside the ordinary positive Fourier cone.  The
full 100-coordinate covariance has rank 100, the minimum cone distance is zero,
and the frozen bootstrap p value is 1.0.  Because residue zero was excluded, the
positive spectral completion is nonunique: adding a common offset to all
Fourier weights changes only the missing `C(0)`.  The fit therefore excludes no
ordinary cone state and identifies neither a unique spectrum nor a state count.

Separately, the N325/N425 source coefficient `alpha=1.78731(1007)` does not
transfer to N505.  The two N505 conditional ratios are `1.41538(1576)` and
`1.41586(1604)`; the stable two-hand transfer statistic is
`chi-square=514.99/2`, `p=1.48e-112`.  That is a scale/geometry amplitude result,
not spatial rejection.

The production question is closed without more samples: adaptive spatial
heterogeneity is present, but no cone-external mechanism is required.
