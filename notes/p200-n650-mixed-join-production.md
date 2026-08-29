# N650 mixed-join production interface

This is the dedicated Phase D runner for the exact `C2 x C5` mixed-join
observable.  It does not add another scale-law coordinate after the N580
q2/Jordan tie.

For both N650 lineages, one counter produces a frozen
`Binomial(650,p_ref)` occupancy count and a conditional Fisher--Yates
permutation.  The same prefix labels feed both HNF geometries.  Each colour
and orientation is evaluated at the four join corners in memory; only the
batch sums of the frozen state are written.

Primary state:

```text
(ES,ED,OS,OD)
```

from the connected residual `R_c=J_full,c-J_local,c`.  Stored integer sums
are twice these averages, avoiding rounding.  The scorer recomputes a single
delete-one `4x4` covariance and zero-vector GLS.  Ambient-H1 mixed rows are
accumulated simultaneously as a secondary diagnostic.

The exact tiny self-test exhausts the real `N1->N2/N5->N10` cover, locks the
partition-residual histogram signs, `499/1024` local mean,
`681/512` odd variance, typed-layer-swap oddness, and both join orders.  One
new caveat emerged: ambient-H1 after an artificial identification depends on
the lift chosen for that identification.  Production freezes the raw
column-HNF representative displacement.  This does not affect the primary
partition-rank residual; ambient-H1 remains convention-labelled secondary
evidence.

The local 2,000-sample smoke already gives an overwhelming primary rejection
of the factor-additive zero bridge, dominated by negative `ES` and `OS`.
Therefore the 100M command is frozen for reproducibility but should normally
be stopped after the mandatory 20k calibration if the scorer returns
`stop_recommended=true`.  A rejection means nonlocal mixed factor
connectivity, not chronological path memory or Jordan structure by itself.

Build and gate:

```bash
g++ -O3 -DNDEBUG -std=c++17 -fopenmp \
  src/p200_n650_mixed_join_mc.cpp \
  -o build/p200/p200_n650_mixed_join_mc
build/p200/p200_n650_mixed_join_mc --self-test
```

The full commands, isolated RNG domain, resource estimate, and stop rules are
frozen in `experiments/p200_n650_mixed_join_100m_20260829.yaml`.
