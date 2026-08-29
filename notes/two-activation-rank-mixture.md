# Exact two-activation semantics of the matching profile

Status: post-#269 exact finite-volume theorem for Issue #28.  No scaling or
single-threshold assumption enters this note.

## 1. A permutation has two ambient-rank activations

Let `v_1,...,v_N` be a uniformly random permutation of the sites and let
`B_n={v_1,...,v_n}`.  Write

```text
V_n = im[H1(K(B_n);Q) -> H1(T^2;Q)],
r_n = dim V_n in {0,1,2}.
```

The occupied cubical complexes are nested.  Every cycle in `K(B_n)` remains
a cycle in `K(B_(n+1))` with the same ambient torus class, so

```text
V_n subset V_(n+1).
```

Attaching new 2-cells cannot remove a nonzero ambient class: if an old cycle
becomes a boundary in the larger subcomplex, it also bounds in the ambient
torus and its ambient class was already zero.  Hence `r_n` is nondecreasing.

The empty and full configurations have ranks zero and two.  Define

```text
K1 = min{n: r_n>=1},
K2 = min{n: r_n=2}.
```

Then `1<=K1<=K2<=N`, with equality allowed when one added site activates two
independent ambient directions at once.  The digital-Alexander theorem gives
`r_white=2-r_n` and the matching sign

```text
q_n=(r_n-r_white)/2=r_n-1.
```

Since an integer-valued nondecreasing rank in `{0,1,2}` is the sum of its two
level indicators,

```text
r_n = 1{n>=K1}+1{n>=K2},
q_n = -1+1{n>=K1}+1{n>=K2}.                            (1)
```

This is configurationwise.  It proves that a permutation has at most two
jumps, including the possible direct jump `-1 -> +1` when `K1=K2`.

## 2. Bernoulli convolution

Generate a Bernoulli-`p` configuration by first drawing a uniform permutation
and then an independent `X~Binomial(N,p)`, occupying its first `X` sites.
For a fixed activation rank `k`, put

```text
H_k(p)=P[X>=k]=sum_(n=k)^N binom(N,n)p^n(1-p)^(N-n).
```

Taking the conditional expectation of (1) and then averaging over the
permutation yields

```text
M_N(p) = -1 + E[H_K1(p)] + E[H_K2(p)],                 (2)

F_N(p)=(1+M_N(p))/2
      = 1/2 E[H_K1(p)+H_K2(p)].                        (3)
```

Therefore `F_N` is exactly the equal mixture of two CDFs:

- onset: the first ambient homology direction activates at `K1`;
- completion: the second direction activates at `K2`.

It is not generically the law of one latent threshold.

## 3. Beta-density mixture

The binomial tail is the CDF of the `k`-th order statistic of `N` independent
uniform labels.  Differentiation gives

```text
dH_k/dp
 = N binom(N-1,k-1) p^(k-1)(1-p)^(N-k)
 = BetaPDF(p;k,N-k+1).                                 (4)
```

Consequently

```text
rho_N(p)=F_N'(p)=M_N'(p)/2
 = 1/2 E[BetaPDF(p;K1,N-K1+1)
        +BetaPDF(p;K2,N-K2+1)].                        (5)
```

The matching root is the mixture-median balance

```text
E[H_K1(p*)]+E[H_K2(p*)]=1.                             (6)
```

Equation (6) does not say that `p*` is the median of either marginal
activation distribution.  It is the median of their equal mixture.

## 4. Midpoint and gap retain the information lost by the mixture

For each permutation define

```text
C=(K1+K2)/2,        G=K2-K1.
```

Then

```text
K1=C-G/2,           K2=C+G/2,
```

so the joint `(C,G)` law is equivalent to the paired `(K1,K2)` law.  The
marginal mixture `F_N` discards their pairing.

The semantics are exact:

- `C` is the activation midpoint in occupation-rank coordinates;
- `G` is the number of microcanonical occupation levels for which `r_n=1`,
  because `r_n=1` precisely for `K1<=n<K2`;
- `G=0` is a simultaneous two-direction activation;
- `E[C]/(N+1)` is the mean of the equal order-statistic mixture;
- integrating the canonical one-cycle probability gives

```text
integral_0^1 P_p(r=1) dp = E[G]/(N+1),                 (7)
```

  since every Bernstein basis polynomial has integral `1/(N+1)`.

Thus the midpoint is a location coordinate and the gap is the exact neutral
one-cycle-window area coordinate.  They should be archived jointly even when
only the one-dimensional mixture is plotted.

## 5. Existing tiny exact histogram certificate

`scripts/two_activation_rank_mixture.py` reuses the existing exact
`threshold_rank_nz.enumerate_exact(axis_integer_torus(2))` histogram.  Its 24
permutations give

```text
(K1,K2)=(2,3): 16,
(K1,K2)=(3,3):  8.
```

The corresponding exact identities are

```text
q_bar[n=0..4] = [-1,-1,-1/3,1,1],
M(p)           = -1+4p^2-2p^4,
F(p)           = 2p^2-p^4,
rho(p)         = 4p-4p^3,

rho = (1/3) BetaPDF(2,3) + (2/3) BetaPDF(3,2).
```

The two-CDF balance root is

```text
p* = sqrt(1-1/sqrt(2)) = 0.541196100146...,
```

which illustrates that the mixture root need not be either marginal median.

The executable oracle independently enumerates all occupied masks, verifies
the microcanonical Bernstein coefficients, checks every rank trace against
(1), and compares the derivative polynomial with the beta mixture at exact
rational probe points.  This tiny quotient lies outside the honest-cell
scope of the topological proof, but the already established finite oracle
and the direct trace checks verify the same identity configurationwise here.
