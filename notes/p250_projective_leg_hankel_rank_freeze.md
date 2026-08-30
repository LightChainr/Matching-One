# P250 model-free multivariate Hankel rank freeze

## Why this test

The common-state rank-1/2/3 exponential fits and the canonical rank-five Weyl
fit are model classes.  This reanalysis removes their spectral assumptions.
For the six monomials of total degree at most two, it constructs

```text
H[(channel,u),v] = G_channel(u+v).
```

Every path-independent finite-dimensional commuting realization has
`rank(H)<=state dimension`, whether its transfer generators are diagonal,
defective, or Jordan.  Rejecting `rank(H)<=5` is therefore a dimension lower
bound, not another root fit.

## Frozen rank chart

For each rank `r=1..5`, choose the maximum-volume `r x r` submatrix `P` of the
full mean.  Hold that coordinate chart fixed in every delete-one sample and
score the complete Schur complement

```text
S - R P^-1 Q = 0.
```

These are a locally complete set of vanishing `(r+1)`-minor constraints.  The
test is run for each charge row separately, for the two-charge plus and minus
blocks, and for the four-channel shared block.  The block formulation allows
channel-specific left functionals/amplitudes but requires one common transfer
state.

Every real and imaginary Schur coordinate enters one joint 400-batch
delete-one covariance.  The primary probability is the finite-batch Hotelling
correction after correlation eigenmodes below `1e-10` of the largest are
discarded; the asymptotic chi-square probability is retained as a secondary
diagnostic.  At `alpha=0.01`, rejection through rank five gives lower bound
six.

## Flat-extension boundary

The available order-one plateau is `rank(H_degree2)=rank(H_degree1)<=3` and is
equivalent here to survival of the rank-three null.  A failure only says that
more than three states are needed.

The stronger order-two plateau requires the degree-three monomial matrix and
therefore moments through total degree six; the current diamond stops at total
degree four.  Likewise, the stream stores one row per endpoint displacement,
so `G_xy` and `G_yx` are not separately observed.  Endpoint Hankel consistency
is an exact construction gate, not evidence against path/context memory.  A
path-dependence claim would require ordered-path rows, not reinterpretation of
this dataset.
