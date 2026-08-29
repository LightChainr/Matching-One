# Square-bond primitive-sector pilot for issue #156

Status: exact tiny oracle plus exploratory 200k-per-geometry control.  This is
not an exponent fit and not an H4-versus-H8 discrimination.

## Design

The calculation uses square-bond percolation at the exact critical point
`p=1/2`.  Open bonds are classified by the exact homology image of the full
configuration in the declared period basis.  The three target engine lines

```text
l0=(1,0), l1=(0,1), l2=(1,-1)
```

map to Pinson--Arguin paper types

```text
{1,0}, {0,1}, {1,1}
```

because the papers use `a*omega1-b*omega2` while the engine uses
`u*omega1+v*omega2`.  Rank zero, rank two, and all other primitive rank-one
directions are retained as separate mutually exclusive categories.

The pilot geometries are the first useful Pell members on opposite sides of
the hexagonal modulus:

| design | period matrix rows | N | tau | samples |
|---|---|---:|---|---:|
| D=-2 | `[[6,3],[0,5]]` | 30 | `1/2+5i/6` | 200,000 |
| D=+1 | `[[8,4],[0,7]]` | 56 | `1/2+7i/8` | 200,000 |

Each design uses 100 independently seeded batches of 2,000 configurations.
The seed is `20260829`.  The full batch counts and covariance-bearing result
are archived under
`results/local-20260829/P156-square-bond-primitive-pilot/`.

## N=4 exact oracle

For the fundamental Pell quotient `P=[[2,1],[0,2]]`, `N=4` and eight bond
variables, all 256 configurations give

| category | exact count |
|---|---:|
| rank zero | 75 |
| `l0` | 57 |
| `l1` | 24 |
| `l2` | 24 |
| rank-one other | 1 |
| rank two | 75 |
| incompatible rank-one directions | 0 |

The single other direction is engine winding `(1,-2)`.  The count proves that
the three shortest sectors are observables, not an exhaustive replacement for
the full rank-one classification.

## Continuum subtraction and real C3 coordinates

Let `r_j=P_j-pi_j(tau_N)`, using the PR #213 general-complex-tau
Pinson--Arguin baseline.  The reported coordinates are

```text
C = r0 - (r1+r2)/2
Q = sqrt(3)/2 * (r2-r1)
S = r0+r1+r2.
```

`C` is the real nontrivial C3 contrast, `Q` is the reflection-odd null, and
`S` is the scalar diagnostic.  The two complex nontrivial DFT modes are
conjugate for real probabilities and are therefore not counted twice.

## Results

| N | empirical `(P0,P1,P2)` | C | Q reflection null | S scalar |
|---:|---|---:|---:|---:|
| 30 | `(0.117800, 0.127500, 0.125975)` | `0.00754883 +/- 0.00099631` (`z=7.58`) | `-0.00132069 +/- 0.00109102` (`z=-1.21`) | `0.00611481 +/- 0.00105245` (`z=5.81`) |
| 56 | `(0.128330, 0.120350, 0.119640)` | `0.00378652 +/- 0.00097123` (`z=3.90`) | `-0.00061488 +/- 0.00104169` (`z=-0.59`) | `0.00331517 +/- 0.00106500` (`z=3.11`) |

The continuum baselines used were

```text
N30: (0.1107291776903850, 0.1272155037499346, 0.1272155037499346)
N56: (0.1247005965715636, 0.1201521157922650, 0.1201521157922650).
```

The reflection-null coordinate is statistically consistent with zero on both
sides, so the winding sign, line order, and covariance transformation pass the
intended control.  The real nontrivial contrast is positive on both sides and
is resolved in this inexpensive pilot.

As a fixed-coordinate diagnostic only, `N*C` is about `0.2265` and `0.2120`.
No exponent is fitted from those two numbers.  More importantly, the scalar
residual is also resolved: `N*S` is about `0.1834` and `0.1856`.  Thus the
continuum approximation has visible finite-size error in both representation
sectors.  This does **not** contaminate `C`: the `(C,Q)` doublet is exactly
orthogonal to the scalar `S` coordinate before any statistical weighting.

The representation statement is consequently stronger than a generic
"positive signal." Under the 60-degree cycle of the three unoriented lines,
spin 12 has phase one and belongs to the scalar coordinate, whereas spin 4 and
spin 8 occupy the two conjugate nontrivial modes.  The resolved `C` therefore
excludes a purely scalar or pure-H12 explanation of this sector residual.  It
does not distinguish H4 from H8, which remain conjugate aliases in the real
three-sector data.

## Boundary and next gate

This result cleanly detects a nontrivial C3 doublet beyond the scalar sum and
shows that the continuum/sign conventions work in production-like data.  It
does not distinguish H4 from H8: the real three-sector probabilities carry one
nontrivial real doublet, while the nominal complex modes are conjugate.  It
also does not identify the separately resolved scalar finite-size term.

A larger N418/N780 run is justified if the next question is the stability and
radial scaling of the nontrivial doublet, independently of the scalar channel.
The present two sizes must not be used for a free-exponent fit.

## Reproduction

```bash
PYTHONPATH=scripts python3 scripts/square_bond_primitive_pilot.py \
  --samples 200000 --batches 100 --workers 10 \
  --seed 20260829 --dps 80 \
  --output-prefix results/local-20260829/P156-square-bond-primitive-pilot/result
```

The CLI appends `.json` and `.batches.csv` to the declared prefix, matching the
committed archive names.
