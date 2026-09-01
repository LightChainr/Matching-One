# P337: ordinary-four-arm landing minor rejects finite pure-thermal rank one

## Decision

`FINITE_PURE_THERMAL_RANK_ONE_REJECTED`.

The canonical pair source and the ordinary Bernoulli thermal score are not
collinear across first essential birth (`0 -> 1`) and completion (`1 -> 2`)
landing sectors. The result survives the existing spin-four landing projection
and subtraction of any common thermal/root counterterm. It reaches the
explicit stop condition in Issue #537 for the exact finite rank-one
cancellation route.

This is a finite algebra/mechanism result. It does **not** prove a nonzero
infinite-volume wedge, identify `x=17/4` or `x=21/4`, or establish the final
rate `T_N=o(D/A_N)`.

## Matrix

For a forced-vacant root `z`, use

```text
S_mid = K_rest + 1/2 - N p,
a_mid = [a(z=0)+a(z=1)]/2,
a = N^-2 sum_(x!=y) g_xy.
```

Let `ell_4=1_axis-1_diagonal` be the already merged fixed-root four-arm
landing character. For `tau in {01,12}` define

```text
T_tau = sum E_-z[ell_4 S_mid 1(tau)],
A_tau = sum E_-z[ell_4 (a_mid-Ea) 1(tau)].
```

The tested matrix is

```text
        0->1       1->2
T     [ T_01       T_12 ]
A     [ A_01       A_12 ].
```

For a pure thermal source, `A` is a common multiple of `T`, so every minor
vanishes. Replacing `a` by `a-beta K+constant` performs the row operation
`A -> A-beta T`; the determinant is invariant. A `pi/4` registry shift
negates the landing character and also preserves nonzero rank. The transition
map to `(q,E)` has determinant two, so an invertible root-conditioned change of
basis does not restore rank one.

## All-scale physical witnesses

Fix `R>=1`, `L>=2R+5`, root `(0,0)`, and `s=R+2`. On the axis `L x L` torus
set

```text
A_LR = {(x,0):1<=x<L} union {(x,s):1<=x<L},
B_LR = {(0,y):1<=y<L} union {(x,s):0<=x<L}.
```

Both rest configurations have `2L-2` occupied sites. Inside the radius-R box,
`A_LR` has two separated black east/west arms and complementary white
north/south arms; `B_LR` has the rotated axis landing. Both have `ell_4=+1`
and exactly two retained black and two retained white landing components.
Adding the root gives

```text
A_LR: rank 0 -> 1,
B_LR: rank 1 -> 2.
```

A direct canonical port-partition count gives only `g=1/4` nonzero unordered
pair terms:

| state | nonzero unordered pairs | source `a` |
|---|---:|---:|
| `A,z=0` | `4L-8` | `2(L-2)/L^4` |
| `A,z=1` | `4L-4` | `2(L-1)/L^4` |
| `B,z=0` | `4L-12` | `2(L-3)/L^4` |
| `B,z=1` | `4L-8` | `2(L-2)/L^4` |

Hence

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
minor is strictly positive. This is not a radius-one artifact: `R` may grow
with `L` while the four landing arms remain separated to radius `R`.

The script reconstructs `(L,R)=(7,1),(9,2),(11,3),(13,4)` and checks the closed
form exactly.

## Complete L4 aggregate

The exhaustive fixed-root enumeration reproduces the merged landing oracle:

```text
pivotal=3121, axis=892, diagonal=474, both=88,
landed=1278, h4=418, registry-shift violations=0.
```

It then aggregates every `ell_4 != 0` state in the two transition sectors. At
`p=1/2`,

```text
T_01 = -527/16384,
T_12 =  337/32768,
A_01 =  16710343/68719476736,
A_12 =  5521655/137438953472,

det = -533831111/140737488355328 != 0.
```

The positive control `a=K` has `A=T` and determinant zero exactly.

A separate implementation used physical-edge adjacency BFS for all source
values and a potential union-find for homology rank. It agreed on all 65,536
L4 states and reproduced the four fractions above exactly.

## Consequence for #537

The proposed exact cancellation of the ordinary three-packet landing block is
false. The surviving thermal-gauge-invariant object should be retained as

```text
Psi_4 = T_01 A_12 - T_12 A_01.
```

The remaining problem is the signed critical/near-critical scaling of the
**aggregate** `Psi_4`, with kernel reconnection, rank/readout pivotality and the
pooled-root Schur term kept together. It is no longer justified to assume all
ordinary landing contributions cancel and jump directly to a four-packet
absolute remainder.

The explicit witnesses have exponentially small individual Bernoulli weight;
they falsify an exact algebraic identity but are not an asymptotic lower bound.
Any claim that rank one emerges only after critical averaging now needs a
separate quantitative suppression theorem.

## Reproduction

```bash
cd experiments/p337-landing-minor-exact-20260901
python3 landing_minor.py --out /tmp/p337-landing-minor.json
python3 -m unittest -v test_landing_minor.py
python3 verify_landing_minor.py --out /tmp/p337-landing-minor-independent.json
```

The scripts use only the Python standard library. No random configurations,
GPU, cloud job, free exponent, distance grid or N25 extension was used.
