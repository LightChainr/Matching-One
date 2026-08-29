# Issue #159: score the primitive-sector Pell/hex bridge

## Why no new pilot was run

PR #222 already contains the exact N=4 primitive-sector oracle and the
maximum requested local pilot: 200,000 square-bond configurations in 100
batches at both `N=30` and `N=56`. Its batch counts are sufficient to
reconstruct all three sector probabilities and their covariance. Repeating
that stream would duplicate evidence, so this scorer consumes those committed
sufficient statistics directly.

This is a **post-reveal operational bridge audit**, not a newly preregistered
score and not independent evidence from PR #222.

## Typed observable and basis transport

The observable is fixed-p square-bond percolation at `p=1/2`. Each complete
bond configuration is classified by its primitive rank-one homology subgroup
in the period-matrix column basis. The ordered positive-rho registry is

```text
l0=(1,0), l1=(0,1), l2=(1,-1).
```

Both Pell geometries use normalized periods
`omega1=1, omega2=1/2+i*y`, with positive imaginary part. Hence transport from
the `D=-2`, N30 basis to the `D=+1`, N56 basis is exactly the identity; no
modular relabeling or sign choice is available after seeing the means. The
positive-rho action

```text
A = [[0,-1],[1,1]]
```

cycles the unoriented registry `l0 -> l1 -> l2 -> l0` on both sides.

After subtracting the Pinson--Arguin continuum sector values at each actual
modulus, the fixed coordinates are

```text
C = r0-(r1+r2)/2       primary nontrivial real character
Q = sqrt(3)/2*(r2-r1) reflection/convention null
S = r0+r1+r2           scalar control.
```

The score JSON contains the full block-diagonal 6x6 covariance in order
`(C30,Q30,S30,C56,Q56,S56)` and the two 3x3 within-design covariances.

## What passes

The exact N=4 oracle passes with counts
`75/57/24/24/1/75` and zero invariant failures. The transported nontrivial
character is decisively nonzero:

```text
C30 = 0.00754883 +/- 0.00099631, z=7.58
C56 = 0.00378652 +/- 0.00097123, z=3.90
joint zero chi2=72.61/2, p=1.71e-16.
```

The reflection null also passes:

```text
Q30 z=-1.21, Q56 z=-0.59
joint zero chi2=1.814/2, p=0.404.
```

The scalar coordinate is separately nonzero
(`chi2=43.45/2`, `p=3.68e-10`). This does not algebraically leak into C,
because C and S are distinct exact representation coordinates.

## High-information phase failure

The important new result is negative. Along `Re(tau)=1/2`, the ordinary
H4/Q4 simple-zero phase reverses between these Pell sides. The task-provided
phase coordinates are approximately

```text
N30: -0.150698
N56: +0.040410.
```

The repository's direct E4 q-series independently gives
`E4(tau)/E4(i)=-0.1490664,+0.0362565`, confirming the same sign reversal.
Because the sector basis transports by the identity, the observed positive
`C30,C56` cannot be repaired by a basis swap. Therefore C fails the
conditional ordinary-H4 simple-zero sign gate. Its same-sign pattern is
compatible with an even-in-E4, H8-like response; the real three-line
character alone already aliases H4 and H8.

The historically frozen conditional square-bond rule from commit `46f3a6f`
would use `A_D=N^2 C_D` and predict `A_Dminus2/A_Dplus1 -> -2`. Applied
post-reveal to this newly proposed C bridge, it gives

```text
A_Dminus2 =  6.79394 +/- 0.89668
A_Dplus1  = 11.87452 +/- 3.04577
observed ratio = +0.57214
A_Dminus2 + 2*A_Dplus1 = 30.54299 +/- 6.15719, z=4.96.
```

This is recorded as an exploratory historical-rule audit, not disguised as a
preregistered C score.

## Decision

The primitive-character bridge passes: it is exact, nontrivial, covariance
aware, and convention controlled. The H4-specific Pell filter fails. The
present C statistic must not promote square-site H4 production under #159.
Scientifically, this distinguishes two claims that PR #222 alone left close
together:

1. a non-scalar C3 finite-size response exists — supported;
2. that response carries the ordinary H4/E4 simple-zero phase — unsupported
   by this observable and pilot.

## Reproduction

```bash
PYTHONPATH=scripts python3 scripts/score_p159_pell_hex_filter.py \
  --batches results/local-20260829/P156-square-bond-primitive-pilot/result.batches.csv \
  --source-result results/local-20260829/P156-square-bond-primitive-pilot/result.json \
  --output results/local-20260829/P159-pell-hex-h4-filter/score.json

python3 tests/test_p159_pell_hex_filter.py
```
