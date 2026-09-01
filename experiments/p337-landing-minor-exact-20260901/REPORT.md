# P337: exact ordinary-four-arm landing minor rejects finite pure-thermal rank one

## Decision

`FINITE_PURE_THERMAL_RANK_ONE_REJECTED`.

The canonical pair source and the ordinary Bernoulli thermal score are not
collinear across the first essential birth (`0 -> 1`) and completion
(`1 -> 2`) landing sectors. This remains true after the existing spin-four
landing projection and after subtracting any common thermal counterterm from
the source. The finite rank-one cancellation route in Issue #537 therefore
hits its declared stop condition.

This is a finite algebra/mechanism result. It does **not** prove that the
resulting wedge has a nonzero infinite-volume limit, identify `x=17/4` or
`x=21/4`, or establish the final rate `T_N=o(D/A_N)`.

## Matrix tested

For a forced-vacant root `z`, let

```text
S_mid = K_rest + 1/2 - N p,
a_mid = [a(z=0)+a(z=1)]/2,
a = N^-2 sum_(x!=y) g_xy.
```

Use the repository's fixed-root ordinary four-arm landing character
`ell_4=1_axis-1_diagonal`. For the two rank transitions `tau=01,12`, define

```text
T_tau = sum E_-z[ell_4 S_mid 1(tau)],
A_tau = sum E_-z[ell_4 (a_mid-Ea) 1(tau)].
```

The transfer matrix is

```text
        0->1       1->2
T     [ T_01       T_12 ]
A     [ A_01       A_12 ].
```

A pure thermal source has `A=T` up to a common scalar and therefore has zero
minor. Replacing `a` by `a-beta K+constant` performs the row operation
`A -> A-beta T`, leaving the determinant unchanged. Shifting the landing
registry by `pi/4` negates both columns and also preserves nonzero rank. The
transition map to `(q,E)` has determinant two, so an invertible change to
`(q,E-Rq)` does not restore rank one.

## All-scale physical witness family

Fix `R>=1`, `L>=2R+5`, root `(0,0)`, and put `s=R+2`. On the axis `L x L`
torus define

```text
A_LR = {(x,0):1<=x<L} union {(x,s):1<=x<L},
B_LR = {(0,y):1<=y<L} union {(x,s):0<=x<L}.
```

Both rest configurations contain `2L-2` occupied sites. Inside the radius-R
box, `A_LR` has two separated black east/west arms and two complementary white
north/south arms; `B_LR` has the rotated axis landing. In both cases
`ell_4=+1`, with exactly two retained black and two retained white landing
components.

Adding the root gives

```text
A_LR: rank 0 -> 1,
B_LR: rank 1 -> 2.
```

A direct port-partition count gives only `g=1/4` nonzero unordered pair terms:

| state | number of nonzero unordered pairs | source `a` |
|---|---:|---:|
| `A, z=0` | `4L-8` | `2(L-2)/L^4` |
| `A, z=1` | `4L-4` | `2(L-1)/L^4` |
| `B, z=0` | `4L-12` | `2(L-3)/L^4` |
| `B, z=1` | `4L-8` | `2(L-2)/L^4` |

Thus

```text
a_mid(B)-a_mid(A) = -2/L^4,
S_mid(A)=S_mid(B)=2L-3/2-L^2 p.
```

The two-state minor is exactly

```text
[p^(2L-2) (1-p)^((L-1)^2)]^2
  * (2L-3/2-L^2 p) * (-2/L^4).
```

For every `L>=4` and `p>=1/2`, the thermal factor is strictly negative, so the
minor is strictly positive. The obstruction is therefore not a radius-one
artifact: `R` can increase with `L` while retaining four separated landing
arms to the declared radius.

The script explicitly reconstructs the family at `(L,R)=(7,1),(9,2),(11,3),
(13,4)` and verifies the closed formulas.

## Complete axis-L4 aggregate and matching-root certificate

The exhaustive fixed-root enumeration independently reproduces the already
merged landing oracle:

```text
pivotal=3121, axis=892, diagonal=474, both=88,
landed=1278, h4=418,
registry-shift violations=0.
```

It then aggregates all `ell_4 != 0` states in the `0->1` and `1->2` sectors.
At `p=1/2`,

```text
T_01 = -527/16384,
T_12 =  337/32768,
A_01 =  16710343/68719476736,
A_12 =  5521655/137438953472,

det = -533831111/140737488355328 != 0.
```

The exact axis-L4 matching polynomial has one root in `[0.59,0.60]`. A
standard-library Sturm certificate finds no determinant root in the same
interval and both endpoint determinant values are negative. At the matching
root midpoint the numerical matrix is

```text
T_01 = -0.03048295815741206,
T_12 = -0.006128907272948605,
A_01 =  0.0002808942330364974,
A_12 =  0.000138537661833922,
det  = -2.501463041122436e-6.
```

The `a=K` positive control gives `A=T` and determinant zero exactly.

An independent verifier recomputes all 65,536 L4 source values with a
physical-edge adjacency BFS and all ranks with a potential union-find. It
reproduces the same four matrix entries and determinant exactly.

## Consequence for Issue #537

The proposed exact rank-one cancellation of the ordinary three-packet landing
block is false. The next object should be retained explicitly as the
thermal-gauge-invariant birth/completion wedge

```text
Psi_4 = T_01 A_12 - T_12 A_01.
```

The appropriate remaining problem is to control the signed critical/near-
critical scaling of the **aggregate** `Psi_4`, including kernel reconnection,
rank/readout pivotality and the pooled-root Schur term. It is no longer
appropriate to assume all ordinary landing contributions cancel and jump
directly to a four-packet absolute remainder.

The all-scale witnesses have exponentially small individual Bernoulli weight;
they falsify the exact algebraic cancellation but are not an asymptotic lower
bound. Any claim that rank one emerges only after critical averaging must now
state and prove that quantitative suppression separately.

## Reproduction

```bash
cd experiments/p337-landing-minor-exact-20260901
python3 -m unittest -v test_landing_minor.py test_root_certificate.py
python3 landing_minor.py --out /tmp/p337-landing-minor.json
python3 root_certificate.py --out /tmp/p337-landing-root.json
python3 verify_landing_minor.py --out /tmp/p337-landing-minor-verify.json
```

The implementation is Python standard library only and refuses to overwrite
its output. No random configurations, cloud job, GPU, free exponent, distance
grid or N25 extension is used.
