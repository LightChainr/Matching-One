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

## Exact obstruction for join-only dynamics

There is a further inexpensive closure.  For every marked-point partition
`tau`, let

```text
J_tau |pi> = |pi join tau>.
```

Then exactly

```text
J_tau J_sigma = J_(tau join sigma) = J_sigma J_tau,
J_tau^2 = J_tau.
```

The connectivity-join algebra is therefore a commutative semilattice algebra
generated by idempotents.  In characteristic zero its generators are
simultaneously diagonalizable, so no element of this join-only algebra can
carry a nontrivial Jordan block.  The direct four-mark oracle verifies all six
pair joins are commuting, idempotent, and self-adjoint for the first-jet form.

This is a useful design obstruction: adding edges, merging blocks, or taking
new linear combinations of join statistics cannot realize the desired
extension.  A successful microscopic generator must contain structure outside
the join semilattice, such as detach/cut operations, a Q-dependent colour
projector, or another operation that changes rather than only coarsens the
connectivity partition.

## Detach is sufficient for an algebraic transient, but not the physical gate

Adding the standard detach move immediately changes the finite algebra.  The
exact deterministic join/detach semigroup has no defective element at two
marks.  At three marks the shortest defective word is

```text
D0, J01, D1
```

(operations applied from left to right).  It maps the chain

```text
{012} -> {02|1} -> {0|1|2},
```

so the differences of consecutive states form a length-two zero-eigenvalue
Jordan chain.  This is the minimal finite positive control that a
non-coarsening history/morphism operation can create nilpotent state.

It does **not** yet pass the physical Gram gate.  Exhausting the entire
deterministic semigroup gives 42 elements at three marks and 1,577 at four;
none of the defective elements is self-adjoint for the first-jet Gram form.
Thus a single bare join/detach word can create algebraic memory, but a
Gram-compatible physical block requires a weighted sum of histories, an
explicit Q-dependent action, or a larger state category.  This sharply
separates "noncommutative morphism exists" from "LCFT Jordan block exists".

## A minimal Gram-compatible weighted-history Jordan block

The weighted-sum option can be realized exactly, still at three marks.  Put

```text
D = detach point 1,
J = join points 0 and 1.
```

On the four-dimensional endpoint radical, define

```text
K = D + J - D J - J D = (D-J)^2.
```

Direct rational calculation gives

```text
K != 0,
rank K = 1,
K^2 = 0,
H K = K^T H.
```

With radical coordinates

```text
v=(1,-1,-1,0),   w=(0,0,0,1),
```

the exact chain is

```text
K w=v,   K v=0,   <v,v>_H=0,   <v,w>_H=1.
```

This is the first finite positive control in this route that simultaneously
has a nontrivial Jordan chain and passes the first-jet Gram gate.  Its form is
suggestive: the nilpotent is a signed, connected-history subtraction of two
incompatible idempotent operations, not a single trajectory.

The signs matter.  `K` is not a positive stochastic transfer matrix, and this
calculation does not identify a continuum LCFT field.  It does establish a
precise microscopic mechanism worth transporting into a physical transfer
construction: subtract disconnected one-step histories from the two possible
ordered detach/join histories.

## Claim boundary

This proves a finite algebra statement.  It does not prove a transfer-matrix
Jordan block, a continuum field identity, or a logarithmic coupling.  Its
positive value is to show exactly why the unadorned first-order Gram extension
cannot decide any of those questions.
