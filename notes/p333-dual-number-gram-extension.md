# The Q=1 connectivity module has a universal first-order radical

This note carries out the first-order-thickening proposal in Issue #333 for
the raw set-partition connectivity Gram module.  The result is exact and more
restrictive than the proposal: the dual-number Gram extension exists, but by
itself it is universal across every nonconstant connectivity direction and
therefore cannot identify a particular LCFT logarithmic pair.

## Exact diagonal form

For set partitions `pi,sigma` of `n` marked points, use

```text
G_Q(pi,sigma) = Q^|pi join sigma|.
```

The incidence/Moebius congruence on the partition lattice gives the exact
diagonal factors

```text
(Q)_k = Q(Q-1)...(Q-k+1)
```

once for every partition with `k` blocks.  Hence the multiplicity of `(Q)_k`
is the Stirling number `S(n,k)`.

Set `Q=1+epsilon`, `epsilon^2=0`.  The one-block factor is a unit,

```text
(Q)_1 = 1+epsilon,
```

whereas every `k>=2` factor has exactly one zero at `Q=1` and

```text
(Q)_k = epsilon (-1)^(k-2) (k-2)!.
```

Therefore the Gram module over `Q[epsilon]/(epsilon^2)` has exactly

```text
one unit pivot,
S(n,k) epsilon-pivots of sign (-1)^(k-2),  k=2,...,n,
no zero pivot at first order.
```

Equivalently, the endpoint Gram radical has dimension `Bell(n)-1`, and its
first-jet pairing is nondegenerate.  This statement is invariant under every
regular change of connectivity basis.

## Direct basis-free checks

The accompanying oracle does not assume the diagonal formula when checking
it.  It expands

```text
G(1+epsilon)=G0+epsilon G1,
```

uses the coordinate-sum-zero radical of the rank-one matrix `G0`, restricts
`G1` to that radical, and computes its inertia by exact rational congruence
elimination.  The direct and Moebius predictions agree:

| marked points | Bell dimension | radical inertia `(+, -, 0)` |
|---:|---:|---:|
| 2 | 2 | `(1,0,0)` |
| 3 | 5 | `(3,1,0)` |
| 4 | 15 | `(8,6,0)` |
| 5 | 52 | `(25,26,0)` |

The absence of a zero direction is the important part; the alternating
signature records cluster-number parity in the exactified basis.

## Mechanism consequence

The raw first-order Gram thickening is **not selective**.  At four marked
points it produces fourteen first-order radical directions, not a preferred
rank-two singlet/[2] extension.  Thus an `epsilon` pivot is a necessary piece
of a confluent construction but is far too common to be an LCFT fingerprint.

This separates three structures that should not be conflated:

1. the universal Gram extension proved here;
2. the action of the partition/transfer algebra on that extension;
3. a spectral collision with nonzero scaling-dimension velocity.

Only the latter two can select the Vasseur--Jacobsen--Saleur logarithmic pair.
For Matching One, the next exact object should therefore be the simultaneous
dual-number representation of the Gram form **and one transfer/algebra
generator**, followed by the similarity class of its action on the radical.
Repeating scalar Q-score derivatives cannot supply that information.

## A sharp selector for a future transfer generator

The universal radical still supplies one exact selection rule.  Let a regular
operator family `T(Q)` be compatible with the Gram form.  On the endpoint
radical, write `H` for the nondegenerate first-jet pairing.  Expanding Gram
compatibility to first order gives

```text
H T0 = T0^T H.
```

Thus `T0` is self-adjoint for the indefinite form `H`.  If

```text
(T0-lambda)v=0,     (T0-lambda)w=v
```

is a nontrivial Jordan chain, then

```text
<v,v>_H
 = <v,(T0-lambda)w>_H
 = <(T0-lambda)v,w>_H
 = 0.
```

So the bottom state of every possible Jordan chain must lie on the exact null
cone of the first-radical form.  The gate is sharp: on the two-dimensional
form `H=[[0,1],[1,0]]`, the Jordan matrix `[[2,1],[0,2]]` is H-self-adjoint and
its bottom vector is isotropic.  The machine artifact includes this rational
oracle.

This changes the next search from "find an epsilon pivot" to "find a transfer
eigenvector that is H-isotropic and has a generalized partner".  At four legs
the signature `(8,6)` makes such vectors possible, but not automatic.

There is also an immediate grade constraint.  In the exactified partition
basis, fixed `k` has a definite first-jet sign
`(-1)^(k-2)`.  Consequently an isotropic bottom must mix block-number grades
of opposite sign; no pure fixed-`k` connectivity direction can be the bottom
of a nontrivial Jordan chain.  The two-mark radical is positive definite, so
this mechanism is impossible there.  Indefiniteness first appears at three
marks, with signature `(3,1)`.  At four marks, any candidate must mix the
positive `k=2,4` sector with the negative `k=3` sector.  This is a concrete
connectivity-tensor test for a future VJS/Matching-One transfer eigenvector.

## Claim boundary

This proves a finite algebra statement.  It does not prove a transfer-matrix
Jordan block, a continuum field identity, or a logarithmic coupling.  Its
positive value is to show exactly why the unadorned first-order Gram extension
cannot decide any of those questions.
