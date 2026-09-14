# Exact L=3 two-birth copula control and marginal-only transport bounds

Date: 2026-09-14. Process-level addendum to draft #773 / issue #778. Exact all-permutation control; no Monte Carlo.

## 1. The static rank polynomials determine both birth marginals exactly

Let `J1<=J2` be the insertion indices at which homology rank first reaches 1 and 2 in a uniformly random site permutation of N sites. For a uniformly random k-subset,

    q0(k)=P(rank=0 | K=k)=C[0,k]/binom(N,k),
    q2(k)=P(rank=2 | K=k)=C[2,k]/binom(N,k).

Monotonicity gives

    P(J1>k)=q0(k),
    P(J2<=k)=q2(k).                                      (1)

Hence the complete one-dimensional marginals of J1 and J2 are already encoded in the static rank-sector occupation polynomials.

The gap

    D=J2-J1

satisfies the previously derived exact identity

    E D = sum_k P(rank=1 | K=k).                         (2)

But (1)--(2) do **not** determine the gap distribution, its variance or the simultaneous-birth atom `P(D=0)`.

## 2. Exact full permutation control at L=3

All `9!=362880` site permutations were enumerated using independent lifted-homology rank lookups for both square NN and the standard six-neighbour triangular torus.

### Square NN L=3

    D histogram counts:
      D=0 :  31104
      D=1 : 181440
      D=2 :  98496
      D=3 :  41472
      D=4 :  10368

Thus

    P(D=0)=3/35,
    E D=3/2,
    E D^2=43/14,
    Var D=23/28.                                         (3)

The birth-index moments are

    E J1=5,
    E J2=13/2,
    Corr(J1,J2)=0.3692744729....                         (4)

### Triangular self-matching L=3

    D histogram counts:
      D=0 :  20736
      D=1 : 186624
      D=2 : 112752
      D=3 :  38880
      D=4 :   3888

Thus

    P(D=0)=2/35,
    E D=3/2,
    E D^2=81/28,
    Var D=9/14.                                          (5)

Here

    E J1=17/4,
    E J2=23/4,
    Corr(J1,J2)=1/5 exactly.                             (6)

Both models have the same mean persistence length at this size but distinct copulas and second moments. One-time rank curves therefore do not determine process-level persistence even when their first integrated rank-one mass agrees.

## 3. What can be said if paired raw data are absent

Given only the exact marginals from (1), any admissible joint law is a nonnegative transport plan

    pi(i,j),  i<=j,

with the fixed J1/J2 marginals. Linear programming gives sharp marginal-only bounds on linear costs of the copula.

For the simultaneous-birth atom:

    square L3:      0 <= P(D=0) <= 3/7,
    triangular L3:  0 <= P(D=0) <= 5/14.                 (7)

The true values 3/35 and 2/35 occupy only a small part of these broad intervals.

For the second moment:

    square L3:      5/2 <= E D^2 <= 59/14,
    triangular L3:  5/2 <= E D^2 <= 27/7,               (8)

while the actual values are 43/14 and 81/28.

Thus a marginal-only archive can support honest optimal-transport bounds, but not a reconstructed copula.

## 4. Chain-count hierarchy: where higher persistence moments live

For a random permutation, rank one persists for exactly D insertion levels. Therefore

    D = sum_k 1{rank at k is 1}.

Consequently

    (D)_m
      = m! sum_{k1<...<km}
          1{rank is 1 at every k1,...,km}.                (9)

Taking expectation shows:

- `E D` is a one-time statistic and is determined by the static rank-one coefficients;
- `E[D(D-1)]` requires a two-time nested-subset count;
- the m-th factorial moment requires an m-time chain count in the Boolean lattice.

This identifies exactly how much extra process information is needed beyond `C[j,k]`. A future transfer/DP need not store full permutations if it can count nested rank-one chains.

## 5. Continuous gap law remains exact once D is known

With iid Uniform[0,1] labels, conditional on `D=d`, the continuous persistence width satisfies

    T2-T1 | D=d ~ Beta(d,N+1-d),                           (10)

with a zero atom for d=0. Hence a paired D histogram determines the full raw-p gap distribution without further simulation.

For example,

    E(T2-T1)=E D/(N+1)=3/20                              (11)

for both L=3 controls, despite their different higher gap laws.

## 6. Consequences for #778 and the arm question

1. If production stores paired `(J1,J2)`, use it directly.
2. If only marginals survive, report transport bounds such as (7)--(8), not a fabricated copula.
3. `P(D=0)` is an explicitly process-level observable. It is the integrated probability that one site insertion jumps rank 0->2, and therefore directly tests the theta/6-arm versus 8-arm fusion geometry from #768/#769.
4. A 6-arm signal in `P(D=0)` can coexist with an L^-4 matching-root shift if the 6-arm contribution cancels in the one-point matching-odd projection. The two claims live at different levels and should be tested separately.

## 7. Boundaries

- L=3 exact permutation results are controls, not an exponent measurement.
- LP bounds use only static marginals and the order constraint `J1<=J2`; actual percolation dynamics imposes additional structure.
- No new production, threshold or field identification is claimed.
