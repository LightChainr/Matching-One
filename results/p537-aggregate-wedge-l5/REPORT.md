# P537: the complete mixed landing response contracts like a local `L^-4` defect from L4 to L5

## Result

The complete axis-L5 (`N=25`) population has now been scored for the same
radius-one axis-minus-diagonal landing character and the same canonical pair
source used by PR #545.  At the 80-digit solution of the exact finite matching
polynomial,

```text
p_L5  = 0.5919882565183338446109686802119288790...
Psi4  = -4.0685187141747587287287153549704744e-7
```

where

```text
Psi4 = T_01 A_12 - T_12 A_01
```

is invariant under the common thermal/Schur shear `A -> A-beta*T`.  The full
finite-population matrix is

```text
             0 -> 1                         1 -> 2
T   -1.8419547511865175e-2       -6.1909017310334680e-3
A   +7.4687524437144010e-5       +4.7190898429077420e-5
```

The L4 control reproduces PR #545's four matrix entries and
`Psi4=-2.5014630411224358e-6`.

The latest PR #545 interface identifies the actual root-conditioned finite
mixed response as

```text
S  = T_01 + T_12
C4 = 2*Psi4/S
G4 = (1/S) * partial_logit(C4)
   = 2*p*(1-p)*(Psi4_prime*S-Psi4*S_prime)/S^3.
```

The same complete L5 population gives

```text
C4 = +3.3063343736796936e-5
G4 = +5.2774014587435467e-3.
```

Thus the finite root-Schur residual and its root-coordinate derivative are
both nonzero at L5; the L4 positive sign is retained.

## New finite-size fingerprint

Define the thermal-orthogonal source coefficient

```text
chi_perp = Psi4 / (T_01^2 + T_12^2)
```

and `A_perp=Psi4/sqrt(T_01^2+T_12^2)`.  These quantities are immune to a
common thermal counterterm and avoid the poles of ratios such as `A_12/T_12`.
The two complete populations give

| quantity | L4 | L5 |
|---|---:|---:|
| `Psi4` | `-2.5014630411e-6` | `-4.0685187142e-7` |
| `chi_perp` | `-2.5874324444e-3` | `-1.0774479208e-3` |
| `L^4 chi_perp` | `-0.6623827058` | `-0.6734049505` |
| `L^6 A_perp` | `-0.3295273782` | `-0.3271420183` |
| `L^8 Psi4` | `-0.1639358819` | `-0.1589265123` |
| `L^2 ||T||` | `+0.4974878954` | `+0.4858028116` |
| `L^6 C4` | `+0.5597088537` | `+0.5166147459` |
| `L^4 G4` | `+3.5861547177` | `+3.2983759117` |
| `L^10 (Psi4_prime S-Psi4 S_prime)` | `-1.4907566175` | `-1.5902413797` |
| `L^2 S` | `-0.5857898469` | `-0.6152612311` |

The raw L4-to-L5 effective powers are approximately

```text
||T||      : L^-2.1065
A_perp     : L^-6.0326
chi_perp   : L^-3.9260
Psi4       : L^-8.1391
C4         : L^-6.3590
G4         : L^-4.3749
```

With only two sizes these are not exponent estimates.  Their joint near-
integer pattern is nevertheless a concrete mechanism fingerprint: the
thermal row behaves approximately as `L^-2`, while the source direction
orthogonal to it carries an additional `L^-4=N^-2`.  The first Schur residual
then behaves near `L^-6`, and its root-coordinate derivative restores two
powers, leaving the complete `G4` close to `L^-4`.  This is exactly the
normalization of a bounded local contribution to the canonical pair average.

The first scale comparison therefore favors the following interpretation:
finite rank two is real, but the radius-one aggregate defect may be a local
`N^-2` correction rather than a persistent macroscopic transmission channel.
The next production run should test the single predeclared primary quantity
`L^4*G4` at a held-out larger size, with `L^4*chi_perp` as an auxiliary
coordinate.  It should not scan descriptors,
collar radii, or alternative minors.

## Exact population and computation

The producer fixes the thermal site `z=0`, enumerates every one of the
`2^24` off-site backgrounds, and evaluates the complete source in both the
`z=0` and `z=1` states.  It retains exact integer polynomial coefficients for
the matching mean, source mean, two landing masses, and two source-midpoint
masses.  Sixteen local shards completed in 11.4 seconds wall time on the Mac.

The scorer uses 80-digit arithmetic only after the complete integer tables
are summed.  L4 is an exact semantic control; L5 is the new scientific result.
No Monte Carlo sample, free exponent fit, collar-radius scan, GPU job, or
descriptor search was used.

## Reproduction

```sh
c++ -O3 -std=c++17 \
  experiments/p537-aggregate-wedge-l5-20260901/aggregate_wedge_exact.cpp \
  -o /tmp/p537-aggregate-wedge

mkdir /tmp/p537-l5-shards
seq 0 15 | xargs -P 16 -I{} /tmp/p537-aggregate-wedge \
  experiments/p537-landing-matrix-preflight-20260901/kernel.tsv \
  /tmp/p537-l5-shards/shard-{}.tsv 5 {} 16

python3 experiments/p537-aggregate-wedge-l5-20260901/score.py \
  /tmp/p537-l5-shards/shard-*.tsv \
  --output /tmp/p537-l5-score.json
```

## Boundary

This is an exact-coefficient two-size finite-population result, scored at
high precision rather than by a new root-isolation certificate.  It does not prove an
`L^-4` asymptotic law, a vanishing infinite-volume wedge, or a continuum-field
identity.  It changes the next question from “is rank one exact?” to the
prospective one-number test “does `L^4*G4` stabilize beyond L5?”
