## #437 / #419: topology survives algebraically; the literal six-level readout is noise-limited

Branch `experiment/p437-high-pass-mc-pilot-20260831`, freeze `7fd8aa0`, local
host amendment/execution `7a9ce54`. The checked results and full 38-coordinate
covariance are in
`results/local-20260831/P437-N112-six-noise-20k/`.

- **Mechanism space changed:** the N112 complex C3 Etop field has an exact
  fifth mixed difference `-1/3` on a declared 32-point configuration subcube.
  It is not degree<=4; its positive high-pass self-energy is strictly nonzero.
- **Measured result:** 20k fresh common-noise replicas give energy
  `-0.01224 +/- 0.50950`, versus `Var(F)=.080647 +/- .000447`.
  Sampling variance inflates by **1.2983 million**. Euler zero control is
  compatible with zero (z=-.093); even the known degree-five positive response
  `.298004` is unresolved (`-8.052 +/- 7.381`).
- **Cost/decision:** 7.88 local wall seconds / 69.20 CPU seconds. With the
  exact energy ceiling `1/9`, the pilot-variance 5-sigma projection is at least
  10.51M samples even at maximal signal. Stop at 20k; change estimator, not N.
- **Observer/sector/source/geometry:** square-bond p=.5; three N112 rho
  children; complex C3 r1 Etop; self source `conj(F(X))`; six nested product
  noise levels. Unconditional HP mean is identically zero and is not a target.
- **Dependency:** one new block `p437-N112-six-noise-fresh-20k-20260831`;
  prior `2402a33` results are not pooled or replayed. Host amendment occurred
  before sampling; HZ authentication repair did not become the workstream.
- **Not proved:** no independent external-source bridge, field identity,
  physical rank, or universal impossibility of efficient high-pass acquisition.
- **Next升权观测:** a positive/conditional-difference representation with
  pointwise low-degree annihilation and measured variance, preserving the
  current spectral-energy estimand rather than silently changing it.
