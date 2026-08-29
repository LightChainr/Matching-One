# Exact FK critical-manifold Q-score

On a finite square graph the random-cluster weight is `Q^k v^b`.  The exact
coordinates

```text
eta = log(v/sqrt(Q)),   t = log Q
```

turn it into the two-parameter exponential family

```text
exp[eta b + t T],       T=k+b/2=J/2,   J=2k+b.
```

Thus the self-dual critical manifold is `eta=0`, and one ordinary
`Q=1,p=1/2` ensemble contains every derivative along its `Q` tangent.  The
integer histogram in `J` is sufficient; no nearby-Q simulation is required
for the derivative.

`scripts/exact_fk_q_score_oracle.py` exhausts the 256 configurations of the
`L=2` square-bond torus, records the full `sqrt(Q)^J` histogram by wrapping
observable, and differentiates each expectation in two independent exact
ways.  Direct differentiation of the numerator/partition-function ratio agrees
through third order with

```text
H1 = X,
H2 = X^2-kappa2,
H3 = X^3-3 kappa2 X-kappa3,
X  = T-<T>.
```

The same enumeration evaluates the mixed thermal/Q derivative with the exact
centered score

```text
(b-<b>)(T-<T>)-Cov(b,T).
```

The open-primal and closed-dual wrapping probabilities agree at `Q=1`, while
their difference has a nonzero critical-manifold Q tangent.  Finite-torus
duality therefore does not erase the parameter-space direction: topology is
already a sensitive covector for the `Q` score.

In fact the complete `L=2` wrapping-difference numerator factorizes in
`x=sqrt(Q)` as

```text
x^5 (x-1)(x+1)(x+2)(x^2+6x+16).
```

The percolation zero `x=1` is simple.  Differentiating this factorization and
dividing by the `Q=1` partition sum gives the exact tangent `69/256`.  This is
an algebraic finite-volume topological response, not a fitted derivative.

This establishes the estimator, not a logarithmic field.  A fixed lattice
observable has only the measure derivative `Cov(O,T)`.  The derivative of its
generic-Q representation projector, normalization, and explicit insertion
must be added before taking a collision limit.  The next exact gate is the
Vasseur--Jacobsen--Saleur energy/two-cluster positive control; only after that
should the common spin-4 differential be applied.
