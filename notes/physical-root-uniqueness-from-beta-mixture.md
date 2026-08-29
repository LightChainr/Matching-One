# The physical matching root is unique and simple

Status: exact consequence of the post-#269 two-activation theorem.  This note
closes physical multiple-root mechanisms; it does not analyze the complex
root cloud.

## Theorem

Let an honest finite periodic torus have `N>=1` sites.  For a uniformly random
Newman--Ziff permutation, let

```text
1<=K1<=K2<=N
```

be the first occupation ranks at which the ambient black homology image has
rank at least one and rank two.  Then

```text
M_N(p)=-1+E[H_K1(p)]+E[H_K2(p)],
H_k(p)=P[Binomial(N,p)>=k].
```

For every `k in {1,...,N}`,

```text
H_k'(p)
 = N binom(N-1,k-1) p^(k-1)(1-p)^(N-k)
 = BetaPDF(p;k,N-k+1) > 0                             (1)
```

for every `0<p<1`.  Therefore

```text
M_N'(p)=E[H_K1'(p)+H_K2'(p)]>0                        (2)
```

throughout the open physical interval.

The empty and full configurations have ranks zero and two, hence

```text
M_N(0)=-1,             M_N(1)=+1.                     (3)
```

Continuity, (2), and (3) imply:

> `M_N` has exactly one real root `p*` in `(0,1)`, and it is simple because
> `M_N'(p*)>0`.

This conclusion does not require randomness in either activation rank.  A
deterministic `K1` or `K2` still contributes a strictly positive beta density.
The only nondegeneracy needed is that the two activations occur at genuine
occupation ranks `1,...,N`.

## Exact degeneracy boundary

For a quotient or observable outside the honest-torus hypotheses, extend the
activation convention by

```text
K=0     : direction is already active at p=0, H_0=1;
K=N+1   : direction never activates,          H_(N+1)=0.
```

Write

```text
a_j=P(K_j=0),          b_j=P(K_j=N+1),   j=1,2.
```

Then

```text
M(0)=-1+a_1+a_2,
M(1)=+1-b_1-b_2,

M'(p)=sum_j sum_(k=1)^N P(K_j=k) BetaPDF(p;k,N-k+1).
```

Thus `M'` is strictly positive on `(0,1)` exactly when at least one activation
has positive interior mass.  A root can occur at `p=0` when `a_1+a_2=1`, or
at `p=1` when `b_1+b_2=1`.  If all activation mass is at `0` or `N+1`, the
derivative vanishes identically and constant/endpoint degeneracies are
possible.  Honest tori exclude all of these cases by (3).

## A computable bracket from `K1<=K2`

Define the two marginal activation CDFs

```text
A(p)=E[H_K1(p)],       C(p)=E[H_K2(p)].
```

Since `K1<=K2` configurationwise and binomial tails decrease with their rank,

```text
A(p)>=C(p)             for every p.                   (4)
```

Let `m1,m2` be the unique solutions of

```text
A(m1)=1/2,             C(m2)=1/2.
```

At the physical root, `A(p*)+C(p*)=1`.  Combining this with (4) gives

```text
A(p*)>=1/2>=C(p*),

boxed: m1<=p*<=m2.                                    (5)
```

The bracket is strict if `P(K1<K2)>0`; it collapses to one point when
`K1=K2` almost surely.  It is directly computable from the two marginal rank
histograms and needs neither a polynomial power-basis conversion nor a
complex root solver.

For the exact axis `L=2` histogram reused by the Issue #28 artifact,

```text
m1 = 0.458803...,
p* = 0.541196...,
m2 = 0.614273...,
```

so the onset/completion bracket is visibly nontrivial.

## What this does and does not close

The beta-mixture theorem closes all of the following mechanisms on `(0,1)`:

- two or more physical real roots;
- a physical tangent/double root;
- a higher-multiplicity physical root;
- a sign reversal of the physical slope.

It says nothing about roots outside `[0,1]` or nonreal roots.  A real
polynomial can be strictly increasing on the physical interval while still
having real roots elsewhere and complex-conjugate pairs.  Moreover these are
zeros of the finite matching polynomial, not automatically Fisher or
Lee--Yang zeros of a partition function.  The theorem therefore sharpens the
physical-root statement without reopening the complex-zero program closed by
the earlier local catalogue analysis.
