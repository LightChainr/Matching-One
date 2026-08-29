# One N580/N650 acquisition for q2, Jordan, and morphism memory

Status: post-reveal freeze for Issue #200.  Existing N145/N290 and P57 data
supply source means and covariance; only N580 and one path-marked N650 block
are new.

## 1. N580 is the radial minimal-polynomial shot

Use the two N145 orientations and their exact dyadic children:

```text
N145: (12,1), (9,8)
N580: (24,2), (18,16)
```

The same endpoint is reached directly by `2i` or by two `1+i` steps through
N290.  The affine-clock theorem gives two fixed full-state predictions:

```text
q2:     x580=-(1/2)x145+(3/2)x290
Jordan: x580=-x145+2x290.
```

The committed P180 result already contains both four-vector means and full
source covariance.  Score all four typed coordinates jointly.  Selecting a
model from `T_Su` alone is not allowed.

## 2. N650 is one endpoint with two marked factorizations

Starting from N65, use

```text
path A: 65 --(2-i)--> 325 --(1+i)--> 650
path B: 65 --(1+i)--> 130 --(2-i)--> 650.
```

Both products equal `3+i`.  The final unmarked designs are exactly

```text
(23,11), (17,19),
```

and `chi4(3+i)=(7+24i)/25`.  There is only one unmarked final graph and one
unmarked histogram.  Treating the two path labels as independent endpoint
replicates would be pseudoreplication.

For each final configuration, additionally retain the allowed intermediate
deck/character flag for both filtrations.  Their difference is the marked
commutator

```text
C_mark=H_[(2-i) then (1+i)]-H_[(1+i) then (2-i)].
```

Ordinary q2 and a continuum Jordan clock both predict `C_mark=0`.  The risky
morphism-memory prediction is a nonzero rank-one vector aligned with the
already frozen P57 conjugation-odd r=2..6 template.  Fit one amplitude to that
fixed direction, then score the four covariance-orthogonal components.  A
nonzero but misaligned defect means higher rank, not the proposed one-state
cover memory.

## 3. Existing norm-5 covariance already fixes the two N650 forecasts

Let `x1,x2,x5` be the width-corrected thermal jets at N65,N130,N325.  For each
order r=2..6 the two source-only continuations are

```text
q2 via (1,5):  x10=(9x5-x1)/8
q2 via (1,2):  x10=(9x2-4x1)/5

Jordan via (1,5): x10=x5+[log2/log5](x5-x1)
Jordan via (1,2): x10=x2+[log5/log2](x2-x1).
```

Their disagreement is not a new fitted statistic.  It is exactly the already
archived P57 cocycle residual:

```text
Delta_path(q2)    =(9/8) R_q2,
Delta_path(Jordan)=[1+log2/log5] R_J.
```

The script propagates the existing five-by-five residual covariance by these
fixed factors.  Descriptively, the archived q2 source block is already strained
(`p=0.0139`), while Jordan remains live (`p=0.0732`).  The future N650 target
must agree with both source paths under the selected model; it does not repair
their existing disagreement by refitting.

## 4. Three-way decision in one acquisition

```text
q2:
  N580 passes (-1/2,3/2),
  N650 passes both q2 source paths,
  C_mark=0.

Jordan:
  N580 passes (-1,2),
  N650 passes both Jordan source paths,
  C_mark=0.

morphism memory:
  one unmarked radial law survives,
  C_mark is nonzero and rank-one aligned with the P57 odd template.

higher rank:
  neither affine law survives, or the marked defect is not template-aligned.
```

This separates a radial minimal polynomial from categorical memory.  A marked
commutator cannot be absorbed by changing the scalar correction exponent.

## Minimal production contract

- N580: the two listed orientations, full threshold-rank histograms, at least
  100 delete-one blocks, counters independent of P50.
- N650: the two listed endpoint orientations on one common stream; archive one
  unmarked histogram plus both intermediate path flags configurationwise.
- Recompute centers, widths, and finite-N r=2..6 jets inside every delete-one
  replicate.
- Score N580 first, then the two N650 endpoint forecasts, then `C_mark=0`, then
  the one-dimensional P57-template alignment.

## Boundary

Gaussian multiplication, affine coefficients, and the transformation of the
existing cocycle covariance are exact.  The two-state lattice closure and the
alignment of a marked defect with the P57 odd template are high-risk,
falsifiable hypotheses.  Existing data are reused only to freeze source
predictions, not counted again as new evidence.
