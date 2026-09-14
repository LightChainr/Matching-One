# Rank-sector monotone likelihood ratio and the asymptotic weakness of total occupancy

Date: 2026-09-14. Exact finite order theorem plus scaling conjecture. Addendum to draft #773 / interfaces to #578 and #775.

## 1. Exact fixed-cardinality monotonicity

Let `C[j,k]` count k-site subsets with homology rank j on a finite N-site transitive torus. Under a uniformly random site permutation, `R_k` is monotone nondecreasing in k. Therefore

    q0(k)=C[0,k]/binom(N,k) = P(R_k=0)

is nonincreasing, while

    q2(k)=C[2,k]/binom(N,k) = P(R_k=2)

is nondecreasing. Hence wherever both are positive,

    C[2,k]/C[0,k] = q2(k)/q0(k)                         (1)

is nondecreasing in k.

This is an exact combinatorial statement; no product-measure FKG or asymptotic input is required beyond monotonicity of rank under adding occupied sites.

## 2. Conditional occupation counts have monotone likelihood ratio order

At any Bernoulli parameter p in (0,1), the conditional mass functions are

    mu_j(k)=P(K=k | rank=j)
           = C[j,k] p^k(1-p)^(N-k) / P_j(p).

Thus

    mu_2(k)/mu_0(k)
      = [P_0/P_2] [C[2,k]/C[0,k]],                       (2)

which is nondecreasing by (1). Therefore

    K | rank=2  >=_MLR  K | rank=0.                      (3)

Consequences include first-order stochastic dominance and, for every increasing function phi,

    E[phi(K)|2] >= E[phi(K)|0].                          (4)

In particular the log-sector-odds thermal derivative is strictly positive whenever both sectors have positive mass.

## 3. The thermal odds slope is an optimal-transport distance

For probability laws on the line, stochastic dominance implies that the Wasserstein-1 distance equals the difference of means. Hence

    W1(mu_2,mu_0)
      = E[K|2]-E[K|0]
      = g1.                                               (5)

In logit coordinate z,

    g(z)=log[Z_2(z)/Z_0(z)],
    g_z=g1,

while the canonical odd coordinate is `b=-g/2`. Therefore

    -2 db/dz = W1(mu_2,mu_0).                             (6)

The speed of the physical rank law along its matching-odd direction is exactly half a one-dimensional optimal-transport separation between the two extreme rank-conditioned occupation laws.

## 4. Bayes classification from K is a single threshold

At a balance root `P0=P2`, the two sector priors are equal. The likelihood ratio for observing K=k is simply `C[2,k]/C[0,k]`, monotone in k. Therefore the Bayes-optimal classifier between rank 0 and rank 2 using **only K** is a single threshold (with a possible tie randomisation at one k).

The optimal error is

    err_K = (1-TV(mu_0,mu_2))/2.                          (7)

No more complicated nonlinear function of total K can improve this binary extreme-sector classification.

This gives a clean information-theoretic meaning to the exact rank-sector polynomials: they determine the complete best performance achievable by any K-only classifier.

## 5. Scaling paradox: absolute separation grows while relative information dies

Near criticality the current scaling program predicts

    g1 = W1(mu_2,mu_0) = O(L^(3/4)).                      (8)

But the unconditional Bernoulli occupation count has standard deviation

    sd(K)=O(L).                                           (9)

If the two conditional laws have the same leading Gaussian variance and differ primarily by the mean shift (8), their standardized separation is only

    delta_L = g1/sd(K) = O(L^-1/4).                       (10)

This suggests the hierarchy

    W1(mu_2,mu_0)          ~ L^(3/4),
    TV(mu_2,mu_0)          ~ L^(-1/4),
    Hellinger^2 / KL scale ~ L^(-1/2),
    I(rank;K)              ~ L^(-1/2).                    (11)

The last line agrees with the independent Fisher-information calculation on #773. The TV/Hellinger statements require a conditional local-CLT/regularity input and are conjectural.

Thus the raw conditional means can move farther apart in absolute site count while K becomes **less** useful as a normalized topology classifier.

## 6. Consequence for conditioning / Rao--Blackwell design

Any estimator state which retains only total K can capture at most the information in the K sigma-algebra. If (10)--(11) hold, this information about the topological rank sector vanishes asymptotically.

This does not mean conditioning on K is useless at finite size. It says that a production strategy seeking an asymptotically nonvanishing topological variance reduction must retain geometry beyond K: cut connectivity, separator data, homology-relevant marks, or an equivalent structured state.

This supplies a theory-level ordering for the #578 hierarchy:

    rank/K-only summaries: calibration / cheap finite-size gain,
    cut-network geometry: candidate asymptotically informative state.

## 7. Small exact controls

At the finite balance roots, the K-only extreme-sector Bayes errors are still small because L=3,4 are far from the weak-separation regime:

| lattice | L | TV(mu0,mu2) | Bayes error |
|---|---:|---:|---:|
| square |3|0.944523|0.027739|
| square |4|0.872118|0.063941|
| triangular |3|1 exactly|0|
| triangular |4|0.868297|0.065851|

The movement **toward** larger error is qualitatively consistent with (11), but these widths are not an asymptotic test.

## 8. New exact checks for #775

Given exact C[j,k], verify for every returned width:

1. `q0(k)` nonincreasing;
2. `q2(k)` nondecreasing;
3. `C2/C0` nondecreasing on common support;
4. report the K-only Bayes threshold and error at the balance root;
5. report W1, TV and Hellinger distance.

A violation of 1--3 is an implementation/convention error, not interesting new physics.

## 9. Boundaries

- MLR order is exact finite combinatorics; the asymptotic TV/KL rates are conditional conjectures.
- Total occupancy K is not the same as the original-U source or a spatial field.
- Vanishing relative K information does not imply rank itself becomes trivial; the near-critical rank law remains O(1).
