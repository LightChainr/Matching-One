# Krawtchouk decomposition: the full K-only information content of the matching curve

Date: 2026-09-14. Exact finite product-measure/Fourier theorem plus a conditional scaling consequence. Addendum to draft #773 / #789 / #784 / #578.

## 1. The K sigma-algebra is one symmetric mode from every Fourier level

Let `N` sites be iid Bernoulli(p), `q=1-p`, and

    phi_i=(omega_i-p)/sqrt(pq),
    phi_S=prod_(i in S) phi_i.

For each `n=0,...,N`, define the normalized symmetric level-n mode

    Psi_n
      = [binom(N,n)]^(-1/2) sum_(|S|=n) phi_S.            (1)

`Psi_n` is a degree-n Krawtchouk polynomial in the total occupation count

    K=sum_i omega_i.

The family `{Psi_0,...,Psi_N}` is an orthonormal basis for all square-integrable functions of K.

Thus conditioning on K keeps exactly one fully symmetric/zero-momentum direction from each Boolean Fourier level and discards every nonsymmetric mode.

## 2. Common-p derivatives are coherent Fourier sums

For the rank observable `X=r-1`, p-biased Fourier coefficients satisfy for every subset S

    Xhat_p(S)
      = (pq)^(|S|/2) E_p[Delta_S X].                     (2)

Let

    M(p)=E_p[X].

Differentiating the common product parameter n times gives

    M^(n)(p)
      = n! sum_(|S|=n) E_p[Delta_S X].                   (3)

Combining (2)--(3),

    sum_(|S|=n) Xhat_p(S)
      = (pq)^(n/2) M^(n)(p)/n!.                          (4)

Therefore the coefficient of X on the Krawtchouk mode Psi_n is exactly

    a_n
      := E[(X-E X) Psi_n]
       = (pq)^(n/2) M^(n)(p)
         / [n! sqrt(binomial(N,n))].                     (5)

So thermal derivatives are not a separate hierarchy from the Boolean spectrum: `M^(n)` is the **coherent zero-momentum sum of Fourier level n**.

## 3. Exact variance explained by arbitrary functions of K

Orthogonal projection onto the K sigma-algebra yields

    E[X|K]-E X = sum_(n=1)^N a_n Psi_n(K).               (6)

Hence

    boxed:
    Var(E[X|K])
      = sum_(n=1)^N
          (pq)^n [M^(n)(p)]^2
          / [(n!)^2 binom(N,n)].                         (7)

The normalized correlation ratio

    eta_K^2
      = Var(E[X|K])/Var(X)                               (8)

is the **maximum fraction of X variance explainable by any nonlinear function of total occupancy K**.

This is stronger than a linear correlation statement and requires only the one-dimensional polynomial `M(p)` as a function of p.

## 4. The level-1 term is the thermal/topological Fisher angle

For n=1, (7) gives

    eta_(K,level1)^2
      = pq [M'(p)]^2 / [N Var(X)].                       (9)

This is exactly the squared correlation between X and total K, previously denoted the thermal-topological Fisher overlap. It is also exactly the spectral-sample mass on singleton sets:

    P(|S_X|=1)=pq [M']^2/[N Var(X)].                     (10)

Thus the earlier Fisher-information result has a direct Fourier meaning: total K sees precisely the fully symmetric singleton spectral mode.

Higher nonlinear functions of K access the symmetric Krawtchouk mode at each higher Fourier level, but no directional/spatial modes within the same level.

## 5. Coherence fractions and K-only information

Let the total Fourier power at level n be

    P_n=sum_(|S|=n) Xhat(S)^2.

Define the coherent fraction

    Gamma_n
      = |sum_(|S|=n) Xhat(S)|^2
        / [binom(N,n) P_n].                              (11)

Then

    a_n^2 = Gamma_n P_n,                                 (12)

and

    eta_K^2
      = sum_n Gamma_n P_n / Var(X).                      (13)

So K-only information is exactly the coherent part of the full spectral sample. The incoherent remainder

    sum_n (1-Gamma_n)P_n                                 (14)

is topological variance that cannot be recovered from K, however nonlinear the K-only readout is.

This is the precise spectral meaning of the invisible-spatial-information phenomenon in #737.

## 6. Exact small controls

Using the exact rank-sector polynomials at the finite balance roots:

| model | L | eta_K^2 | level-1 share | higher symmetric share |
|---|---:|---:|---:|---:|
| square NN |3|0.7109291021|0.6539241900|0.0570049121|
| square NN |4|0.6179563853|0.5803884127|0.0375679726|
| triangular self-matching |3|0.7773235800|0.7124435241|0.0648800559|
| triangular self-matching |4|0.6577250033|0.6154710928|0.0422539105|

For self-matching triangular site all even Fourier levels vanish exactly, so its K-only decomposition uses odd Krawtchouk modes only.

The leading normalized contributions at L=4 are:

square:

    n=1: 0.5803884
    n=2: 0.0020245
    n=3: 0.0293998
    n=4: 0.0019773
    n=5: 0.0027878

triangular:

    n=1: 0.6154711
    n=3: 0.0362608
    n=5: 0.0050548
    n=7: 0.0007965
    n=9: 0.0001209.

These tiny sizes are structural controls, not asymptotic measurements.

## 7. Conditional near-critical consequence

Assume the leading continuum rank mismatch has an odd thermal scaling function in a smooth microscopic thermal coordinate, with `y_t=3/4`. Then at the balance point the derivative hierarchy has the candidate orders

    M^(2m+1)=O(L^((2m+1)y_t)),
    M^(2m)  =O(L^((2m-1)y_t)).                           (15)

For fixed n, `binom(N,n)~L^(2n)/n!`. Substituting (15) into (7) gives

odd n:

    a_n^2 = O(L^(-n/2)),                                 (16)

even n:

    a_n^2 = O(L^(-n/2-3/2)).                             (17)

Hence the K-only explained variance should be asymptotically dominated by the linear singleton mode:

    eta_K^2
      ~ pq [M']^2/[N Var(X)]
      = O(L^-1/2),                                       (18)

while genuinely nonlinear K-only improvements vanish faster.

This sharpens the previous information-loss conjecture: **not only does K become weak; asymptotically almost all K-only information should already be contained in the linear total-occupancy score.**

## 8. Two scalar p-curves contain all K-only rank information

`M(p)=E[X]` across all p determines the fixed-cardinality conditional mean

    m_k=E[X|K=k]

by Bernstein inversion.

A second scalar curve

    S(p)=E[X^2]=P0+P2=1-P1                            (19)

determines

    s_k=E[X^2|K=k].

Since X has only three values,

    P(r=2|K=k)=(s_k+m_k)/2,
    P(r=0|K=k)=(s_k-m_k)/2,
    P(r=1|K=k)=1-s_k.                                   (20)

Thus the pair of one-dimensional curves `(M(p),S(p))` is equivalent to the complete K-only three-rank table `C[r,k]` after Bernstein inversion. Anything not determined by these two curves is necessarily spatial/connectivity information beyond K.

## 9. Interfaces

- #789: use (7)--(18) as the main exact information decomposition; report `eta_K^2` and Krawtchouk level contributions for #775 widths.
- #784: report both total Fourier level power `P_n` and coherence `Gamma_n`; this separates spatial spectral mass from K-visible mass.
- #578: a K-only Rao--Blackwell state cannot outperform the exact projection (6); geometry-rich conditioning is required to access the incoherent spectral remainder.
- #737: invisible spatial marks live precisely outside the two scalar K-only curves described in section 8.

## 10. Boundaries

- Equations (1)--(14),(19)--(20) are exact finite product-measure identities.
- Scaling statements (15)--(18) require the near-critical derivative hierarchy and are conjectural for square-site until those regularity inputs are proved.
- `eta_K^2` is variance explained for X, not mutual information; converting its asymptotic order into a KL/MI law needs additional distributional regularity.
