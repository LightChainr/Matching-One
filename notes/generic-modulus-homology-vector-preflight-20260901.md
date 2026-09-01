# Generic-modulus full homology-vector preflight

This closes the no-sampling design step for the stabilizer-free route in the
2026-09-01 compute portfolio.  It does not launch the PR #546 queue.

## Concrete geometry

Use the area-100 base period matrix

```text
P = [[10,3],[0,10]],  tau=3/10+i.
```

Exact Gram-automorphism enumeration gives the orientation-preserving lattice
stabilizer `{+I,-I}`.  The same exact result holds for four equal-area signed
stencil cells:

```text
shear -/+          [[10,2],[0,10]], [[10,4],[0,10]]
oblique-aspect -/+ [[-10,0],[-3,-10]], [[-9,-7],[4,-8]]
```

The first span is `0.2`; the second is
`0.0443582710678 + 0.113496642391 i`.  They are linearly independent and all
five cells remain well away from the square/hexagonal elliptic stabilizers.

## Observable and normalizer contract

At every aligned batch and `p` save

```text
P0,
P1[(a,b)] for every observed primitive unoriented winding line,
P1_tail,
P2,
d_p of every preceding column.
```

Then compute `q=P2-P0`, `E=P0+P2`, `d_p q=d_p P2-d_p P0`, and
`d_p E=d_p P0+d_p P2`.  Refind the pooled root `q(p*)=0` inside the same
covariance replicate, set `D=d_p q(p*)`, subtract the saved continuum
homology-vector tangent, and form

```text
U_h = N^(13/8) directional_contrast[d_p P_h] / (2 D)
```

for every homology channel `h`.  This keeps P0/P1/P2, the rank-1 winding
vector, q/E, the normalizer, pooled root and U in one typed pipeline.

## Theory-vector result

The deterministic Q=1 continuum calculation keeps 16 explicit primitive
lines plus a declared tail channel and the rank-0/rank-2 completion.  It
freezes three stacked two-direction templates:

- `mu_KdV`: directional transport of the reflection-even `D2 D0` response;
- `mu_Q4_Jordan`: directional transport of the reflection-even weight-8
  covariant template `E4 * D2 D0`;
- `mu_embedding`: zero after exact continuum-tangent subtraction.

Using the continuum multinomial covariance as a per-sample design proxy and
profiling one continuum-tangent leakage column per direction, the nonzero
templates have

```text
|cos(mu_KdV,mu_Q4_Jordan)| = 0.6681699700
unit-Fisher shape D^2      = 0.6636600599
```

Thus the two operator templates are not rank-degenerate in this vector
design.  The historical N100 one-thread timing gives a four-cell lower-bound
cost of `0.0112135142` CPU-hours per one-million-sample cell-equivalent and a
conditional shape-value proxy `59.1839497` per CPU-hour.  This is a timing
lower bound because sparse full-vector output overhead has not been measured.

The strict three-model maximin portfolio `V` is still zero: the embedding
template is exactly zero after continuum subtraction and no nonzero physical
operator-amplitude floor is supplied.  That zero is distinct from the positive
KdV-versus-Q4/Jordan shape separation.

## Outcome

The route is ready for one bounded N100 sparse full-vector covariance pilot on
the four frozen cells.  It is not yet a large-production authorization.

The requested Huawei machine TV2N0X was started, but the existing key was
rejected with `Permission denied (publickey)`.  No key reset was attempted.
The owned tunnel was stopped, TV2N0X was returned to `Ready`, and the exact
preflight was completed locally in 12.68 seconds with zero random samples.

Reproduce once with:

```text
python3 scripts/generic_modulus_homology_vector_preflight.py \
  analysis/generic_modulus_homology_vector_preflight.json \
  --json-output results/generic-modulus-homology-vector-preflight/latest.json \
  --report-output results/generic-modulus-homology-vector-preflight/REPORT.md
python3 -m unittest discover -s tests \
  -p 'test_generic_modulus_homology_vector_preflight.py' -v
```
