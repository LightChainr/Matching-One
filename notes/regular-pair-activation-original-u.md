# A regular local interaction and its forced Q-activation packet

This is the fixed continuation of the completed
[local pair insertion](local-four-port-original-u-interface.md), not another
fit to its occupation archive. The question is whether that first tangent
extends to a coefficientwise-Q1-regular microscopic tensor family, and
what the canonical singlet completion does to the original observer.
The [contract](../analysis/p337_regular_pair_activation_contract.json)
fixes the completion and its one numerical question before new counts.

## 1. Fix the completion before any observable is scored

On an ordered cut write `Pi_j=i P_j i^dagger` for the same symmetric
unequal-pair embedding used previously. Define

```
Kreg = average_C4 (Pi_2+Pi_0) = average_C4 i(I-P1)i^dagger.
```

At each vertex add `epsilon*Kreg/N` to the vacant summand. All rank
observables still refer to the original occupation, not to colour-diagram
connections. The required algebra and actual four-path lattice witness
are independent gates, and the numerical source is not revised if either
fails the predeclared contraction table.

Regularity in the restricted family `K2+c(Q)K0` fixes only `c(1)=1`.
We select `c(Q)=1` identically. This is the representation-theoretic
complement of the standard block, not a claim that regularity selects a
unique microscopic interaction. In particular a finite counterterm
`c(Q)=1+alpha*(Q-1)+...` can change a Q derivative even when it does not
spoil regularity.

For a vacant x, partition its four incident edge-nodes by their components
in the exterior hypergraph. An edge from x to a vacant neighbour is an
independent singleton. It must not be discarded, and equal colour labels
must not be mistaken for equal exterior components.

The contract fixes the following derivative of the relative colour weight:

| Exterior component pattern | `a_x = d_logQ beta_reg,x(1)` |
|---|---:|
| Four distinct components | 1 |
| One equal pair, opposite NS or EW | 1 |
| One equal pair, adjacent | 1/2 |
| Two pairs, NS\|EW | -1/2 |
| Two pairs, NE\|SW or NW\|ES | -1/4 |
| Three-plus-one or all four equal | 0 |

For occupied x the source is zero. The source is the site average
`a(A)=sum_x a_x(A)/N`, not an extensive sum. The fixed origin supplies
`<a>,<qa>,<Ea>` and their thermal derivatives by exact translation
invariance. This shortcut applies only to first-source joint moments;
it is not a formula for two-insertion moments.

## 2. Why this is a full mixed derivative, not a raw attribution

Let `s_Q(A)=sum_x beta_reg,x(Q,A)/N` be the relative first epsilon
coefficient of the finite occupation weight. The colour algebra will
establish, configuration by configuration,

```
s_1(A)=0,                 d_logQ s_Q(A)|1 = a(A).
```

For any unchanged original observable O, at fixed thermal coordinate h,

```
d_epsilon <O> = Cov_Q(O,s_Q).
d_logQ d_epsilon <O>|1 = Cov_1(O,a).
```

The second identity includes normalization. Derivatives of the base
measure or its normalization multiply the identically zero source s_1;
only the derivative of the source remains. The same argument holds after
any h derivatives. Since `s_1=0` for every h, moving the base Q-dependent
root also multiplies a zero first-source packet. This is why the mixed
derivative simplifies; it is not permission to omit root motion caused
by a itself.

Use the original common-root notation

```
M=(<q>_axis+<q>_tilted)/2,     Y=(<E>_axis-<E>_tilted)/Delta,
D=M_h,                        R=Y_h/D,
Delta=1152/625,                A25=25^(13/8)/2.
```

Define `jM` as the geometry-average covariance with a and `jY` as the
geometry difference of E covariances divided by Delta. Then at the saved
original Q1 root

```
W/A25 = jY_h/D - Y_hh*jM/D^2
        - R*jM_h/D + R*M_hh*jM/D^2,
h_(logQ,epsilon) = -jM/D.
```

Here `W=d_logQ d_epsilon U|Q1,epsilon0` along the original common-root
family. All four terms are mandatory. In contrast,
`d_epsilon U|Q1=0` is exact for the completed interaction and is not the
numerical null being tested. At Q1 the derivatives with respect to Q and
logQ coincide at this order.

## 3. The remaining completion freedom has a named response

For `c(Q)=1+alpha*(Q-1)+...`, the first Q-activated source changes by
`alpha*t`, where t is the old nonnegative local pairing mark. Since the
previously scored source was `S=-t`, linearity gives

```
W_alpha = W_canonical - alpha*V_old,
V_old = +0.0018155512845251097  (site-average N25 units).
```

No alpha is estimated or selected here. This identity exposes the finite
counterterm ambiguity instead of turning it into another fit. A nonzero
canonical W would establish a concrete regular local Q-activation route;
it would not make that number independent of the microscopic completion.

## 4. Frozen single calculation

Only the missing integer `4a, q*4a, E*4a` sums are to be collected on each
complete N25 population, after the algebra, geometry and producer commits.
The original root, slope, U/A and q/E polynomials are imported unchanged.
There is no new root solve, random block, Q grid or finite-counterterm
scan. A strictly signed rational interval rejects the zero mixed response
for this one completion; an interval containing zero ends this question
unresolved at the already fixed arithmetic budget.
