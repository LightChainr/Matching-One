# Q-velocity spectroscopy for the two live spin-4 fields

## Exact target

Let `u=beta^2` and use the dense/FK convention

```text
sqrt(Q) = -2 cos(pi u),       Q=1 at u=2/3.
```

Then `dQ/du=2*pi*sqrt(3)` at percolation. The four-leg primary has

```text
x_22(u) = 1 + (3/2)(u+1/u),
```

while the thermal dimension is `x_epsilon(u)=3/(2u)-1`. A fixed level-4
descendant changes `x` by 4 but does not change its Q velocity. Therefore

```text
V_(2,2):       x=17/4,  dx/dQ=-5 sqrt(3)/(16 pi) = -0.1722902798...
thermal Q4:    x=21/4,  dx/dQ=-9 sqrt(3)/(16 pi) = -0.3101225037...
separation:                  sqrt(3)/(4 pi)       =  0.1378322239...
```

These are exact continuum-family fingerprints. PR #260 already established the
Q=1 dimension gap; this note adds its transverse velocity rather than counting
the same spectrum fact again.

## Amplitude-free two-size scorer

For one consistently defined generic-Q field,

```text
O_Q(L)=A(Q) L^[-y(Q)],
dQ log O(L)=dQ log A - x'(Q) log L.
```

At two sizes the unknown normalization derivative cancels:

```text
V_Q(L2,L1) = -[(dQ O2/O2)-(dQ O1/O1)] / log(L2/L1).
```

The implementation consumes the full covariance of
`[O1,dQO1,O2,dQO2]` and propagates it with the exact nonlinear gradient. Its
synthetic oracle changes both `A(1)` and `dQ log A` and still recovers each
target to roundoff.

This cancellation is a one-field statement. In a two-field mixture there are
two unrelated amplitude derivatives unless a representation calculation ties
them together. Two-size velocity must not be sold as mixture tomography by
quietly dropping those terms.

## Why the present archives cannot be scored

`results/fk-q-score/latest.json` contains an exact L=2 critical-manifold
measure score for fixed topology functions. It has neither a second size nor
the explicit Q derivative of a spin-4 projector/field normalization.

The N130/N170 site-percolation files containing `score_t` and `score_lambda`
belong to a Bernoulli anisotropy parameterization. They are not the FK Potts-Q
score `T=k+b/2`.

Hence the current lattice output is deliberately `NOT_SCOREABLE`. In general,

```text
dQ O_field
 = Cov(O_field, T)                 # measure score
 + E[dQ O_field_definition].       # projector/normalization/insertion
```

The second row is exactly where the unresolved physical Potts multiplicity and
selection rule live. Replacing it by zero would prejudge the question.

## Minimal acquisition

Use square-lattice critical FK at Q=1 and two sizes `L=8,16`. For each of two
controlled field families archive

```text
observable_sum,
measure_score_T_sum,
observable_times_score_T_sum,
explicit_field_definition_derivative_sum,
```

with aligned batches and joint covariance. The two fields are:

1. the character-weighted four-leg `V_(2,2)` insertion;
2. the energy field with the same source-frozen Q4 differential applied across Q.

Samples remain variance-selected. This is a semantic acquisition first: a
billion replicas of the wrong fixed observable would not estimate a field
velocity.

## Reproduction

```bash
python3 scripts/q1_spin4_velocity_oracle.py \
  --output results/q-velocity-spin4/latest.json
```

Once an input satisfying `matching-one.q-velocity-two-size-input.v1` exists,
pass it with `--input`. The scorer refuses measure-only inputs.
