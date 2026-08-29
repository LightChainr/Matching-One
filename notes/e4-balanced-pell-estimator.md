# E4-balanced Pell estimator

This freezes the deterministic part of Issue #160 after its two stated gates: the norm-5 score and the first Pell modular filter are now complete.

For each Pell pair straddling the hexagonal elliptic point, the script evaluates

`a_4 = L^-4 E4(1/2 + i x/(2m))`

from a fixed 96-term q-series and constructs positive normalized weights that annihilate `a_4`. The calculation reads no root estimate. It therefore supplies a prospective estimator rather than a post-target fit.

The generated table records the exact Pell residuals, E4 signs, weights, numerical cancellation residual, the remaining `L^-7` coefficient, and its generation-to-generation ratio. The weights rapidly approach `w_- = 0.071453117982...` and `w_+ = 0.928546882018...`; the scalar ratio approaches `(2+sqrt(3))^-7`.

Execution remains gated by covariance discipline: reconstruct both roots and combine them inside every delete-one replicate. First score the frozen H4 null. Only after that passes may the predeclared scalar law be scored. Failure of the H4 null diagnoses an extra spin-4/topological contribution and invalidates the scalar interpretation, not the arithmetic.
