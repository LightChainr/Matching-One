# P28 production threshold-profile tail elimination freeze

This is the empirical production step that the exact synthetic profile tools deliberately left open.
It consumes the committed P43/P49/P50/P57 batch histograms and launches no Monte Carlo.

The four source sizes (`N=130,145,170,185`) freeze the following choices before any held-out
tail-model score is generated:

- the equal-weight `K_minus/K_plus` Beta order-statistic mixture;
- mixture mean and standard deviation as the location/scale convention;
- the two-sided grid `|z|=2.5,2.75,3,3.25,3.5`;
- at least 10,000 expected mixture thresholds beyond the outer endpoint per curve and at least
  100 per batch;
- separate left/right nuisance intercepts and slopes;
- fixed tail exponents `4/3`, `2` (Gaussian), and `1` (exponential).

The held-out sizes are `N=265,290,325,425`.  Each orientation and each tail contributes the
three parameter-free curvature contrasts left after fitting an intercept and decay constant.
Within a size, all 100 batches are deleted synchronously across both orientations and both tails,
so the score retains the same-curve and shared-permutation covariance.  Held-out archives have
disjoint RNG domains and enter as independent covariance blocks.

A model is `eliminated` at frozen `p<0.01`, `survives` otherwise, or `underpowered` when a count
gate fails or fewer than 75% of the nominal covariance directions are resolved.  `Survives` means
only not rejected.  This is a Level-S production statistical elimination score, not an exact tail
theorem or a cross-microscopic universality proof.
