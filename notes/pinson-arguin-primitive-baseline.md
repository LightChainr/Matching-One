# Pinson--Arguin primitive homology-sector baseline

Status: continuum exact formula / high-precision numerical control for issue
#156.  This is not a finite-lattice result and not an H4 identification.

## Observable and convention

Write the torus as

```text
C / (Z + tau Z),  Im(tau)>0.
```

Arguin and Morin-Duchesne--Saint-Aubin label a primitive type `{a,b}`
by a curve in the physical class

```text
a*omega1 - b*omega2.
```

The minus sign is part of the paper convention.  Repository winding vectors
use `u*omega1+v*omega2`, so the interface map is

```text
engine (u,v) -> paper {u,-v},
```

followed by gcd saturation and the unoriented identification
`{a,b}={-a,-b}`.  In particular, the positive-hexagonal Pell lines

```text
engine: (1,0), (0,1), (1,-1)
paper:  {1,0}, {0,1}, {1,1}
```

are the three equal shortest sectors at `tau=1/2+i*sqrt(3)/2`.

The event is the homology image of the full FK configuration.  A type
`{a,b}` configuration may contain multiple parallel nontrivial clusters, but
its image has rational rank exactly one.  Rank-zero and rank-two cross
configurations are separate events.  Coarse `either`, `both`, or directional
flags are not this observable.

## Q=1 direct formula

The percolation specialization has `e0=2/3` and Coulomb-gas coupling
`g=8/3`.  Put

```text
y = Im(tau)
u = a - b*tau
A = 2*pi*|u|^2/(3*y)
C = sqrt(2/(3*y))/|eta(tau)|^2.
```

Then the normalized critical continuum probability is

```text
pi_tau({a,b}) = C sum_{k in Z} exp(-A*k^2)
                  [cos(2*pi*k/3)-cos(pi*k)].
```

The `k=0` term vanishes.  No additional partition-function divisor is
needed at `Q=1`.  If the direct sum is truncated at `|k|<=K`, the implementation
uses the rigorous absolute bound

```text
4*C*exp[-A*(K+1)^2] / (1-exp[-A*(2*K+3)]).
```

The Dedekind product is evaluated after exact `S`/`T` reduction of its
absolute square.  Its omitted log-product tail is bounded explicitly.

## Independent theta form

For `X=y/|a-b*tau|^2`, with

```text
theta3(z) = sum_n exp(i*pi*z*n^2)
theta2(z) = sum_n exp(i*pi*z*(n+1/2)^2),
```

the same probability is

```text
[theta3(iX/6)-theta3(3iX/2)-2*theta2(3iX/2)]
----------------------------------------------------------------.
                    2*|a-b*tau|*|eta(tau)|^2
```

`scripts/pinson_arguin_primitive.py` evaluates this through independent
imaginary-modulus theta sums and modular transformations.  Tests require the
direct and theta paths to agree at arbitrary complex `tau`.

The exact modular laws in the paper convention are

```text
pi_tau({a,b}) = pi_(tau+1)({a+b,b})
pi_tau({a,b}) = pi_(-1/tau)({-b,a}).
```

They are tested both as paper labels and through the repository period-basis
law `P'=P U`, `w'=U^-1 w`.

## Frozen numerical controls

The machine-readable values are in
`predictions/p156_pinson_arguin_baselines_20260829.json`.  Selected values are:

| modulus | paper sector(s) | probability |
|---|---|---:|
| `i` | `{1,0}`, `{0,1}` | `0.16941543532134688938260796919875445000145337645375` |
| `i` | `{1,1}`, `{1,-1}` | `0.020979928575590629661470611008187000992928553334531` |
| `1/2+i` | `{1,0}` | `0.16815464971788045003554385835401755689282702879997` |
| `1/2+i` | `{0,1}`, `{1,1}` | `0.10005678718797952632971693447464511165342474778502` |
| `1/2+i` | `{1,-1}` | `0.0015190922810096831710181718634974911858075569375806` |
| `1/2+i*sqrt(3)/2` | `{1,0}`, `{0,1}`, `{1,1}` | `0.12166379946598032273800506616100306406127073342455` |
| `1/2+5i/6` (N30) | `{1,0}` | `0.11072917769038501312751860913865363519163224173796` |
| `1/2+5i/6` (N30) | `{0,1}`, `{1,1}` | `0.12721550374993464294805538691418227631716778412260` |
| `1/2+7i/8` (N56) | `{1,0}` | `0.12470059657156361356306944379819518283006971050975` |
| `1/2+7i/8` (N56) | `{0,1}`, `{1,1}` | `0.12015211579226502308743616256051247672799257647673` |

At `tau=i`, summing all unoriented primitive sectors gives

```text
0.38094744914033735446061273329420244691096663128102.
```

For `Q=1`, trivial and cross weights are equal, so each remaining sector has
probability

```text
0.30952627542983132276969363335289877654451668435949.
```

This supplies a normalization and cross-topology regression without treating
cross configurations as primitive rank one.

## Interpretation boundary

The formulas are critical continuum FK probabilities.  They apply to
isotropic site and bond percolation through `Q=1` universality; square-bond
percolation at `p=1/2` is therefore the clean first control.  They are not
exact finite-lattice probabilities, do not apply off criticality, and require
the physical isotropic modulus rather than a bare combinatorial aspect ratio
when the lattice metric is anisotropic.

Nothing here derives a joint primal/matching `S/D` identity.  The first #156
use is sector-by-sector subtraction from square-bond data at the actual Pell
modulus.

## Primary sources

- L.-P. Arguin, *Homology of Fortuin--Kasteleyn clusters of Potts models on
  the torus*, arXiv:hep-th/0111193, DOI 10.1023/A:1019979326380.
- A. Morin-Duchesne and Y. Saint-Aubin, *Critical exponents for the homology
  of Fortuin--Kasteleyn clusters on a torus*, arXiv:0812.2925,
  DOI 10.1103/PhysRevE.80.021130.
- H. T. Pinson, *Critical percolation on the torus*, Journal of Statistical
  Physics 75 (1994), DOI 10.1007/BF02186762.

