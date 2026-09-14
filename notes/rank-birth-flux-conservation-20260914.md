# Rank-birth flux conservation: the one-step process information hidden from static rank marginals

Date: 2026-09-14. Exact finite monotone-process addendum to draft #773 / issues #769 and #778.

## 1. One-step rank dynamics has one hidden flux per insertion level

Let `R_k in {0,1,2}` be the homology rank after the first k sites of a uniformly random permutation have been inserted. Put

    q_j(k)=P(R_k=j).

Monotonicity allows only

    0->0, 0->1, 0->2, 1->1, 1->2, 2->2.

Define the unconditional step-k fluxes

    A_k=P(R_k=0,R_{k+1}=1),
    B_k=P(R_k=0,R_{k+1}=2),
    C_k=P(R_k=1,R_{k+1}=2).

Then the static marginals give exactly

    e0_k := q_0(k)-q_0(k+1) = A_k+B_k,
    e2_k := q_2(k+1)-q_2(k) = B_k+C_k,                  (1)

and

    q_1(k+1)-q_1(k)=e0_k-e2_k=A_k-C_k.                 (2)

Therefore **all one-step rank dynamics beyond the static marginals is carried by the single scalar B_k**:

    A_k=e0_k-B_k,
    C_k=e2_k-B_k.                                       (3)

Positivity gives the immediate marginal-only bound

    0 <= B_k <= min(e0_k,e2_k).                         (4)

The expected rank increment does not expose B_k:

    E[R_{k+1}-R_k]
      = A_k+2B_k+C_k
      = e0_k+e2_k.                                      (5)

This explains algebraically why every one-time rank expectation can be correct while the simultaneous-birth process remains unidentified.

## 2. B_k is the exact rank-jump-two pivotal observable at fixed cardinality

For a transitive N-site torus, fix one site v. Let `n_{2,k}` be the number of k-subsets `S` of the other N-1 sites such that

    rank(S)=0,
    rank(S union {v})=2.

Then

    B_k = n_{2,k}/binom(N-1,k).                          (6)

Proof: among all ordered pairs `(S,v)` with `|S|=k`, there are `binom(N,k)(N-k)=N binom(N-1,k)` possibilities; translation makes the jump-two count N times the fixed-v count.

Consequently

    P(D=0)=sum_{k=0}^{N-1} B_k,                         (7)

where `D=J2-J1` is the rank-one persistence length in insertion steps.

The Bernoulli-p formulation is identical:

    P(D=0)
      = N integral_0^1 P_p(Delta_v X=2) dp,             (8)

because `integral_0^1 p^k(1-p)^(N-1-k) dp = 1/[N binom(N-1,k)]`.

Thus #769's cardinality-resolved jump-two histogram is exactly the missing one-step process flux required by #778.

## 3. Exact L=3 flux controls

Using the exact rank-sector coefficients and independent lifted-homology transition enumeration:

### Square NN L=3

| k | e0_k | e2_k | B_k | A_k | C_k |
|---:|---:|---:|---:|---:|---:|
|2|1/14|0|0|1/14|0|
|3|3/14|0|0|3/14|0|
|4|5/14|1/14|1/70|12/35|2/35|
|5|5/14|5/14|1/14|2/7|2/7|
|6|0|4/7|0|0|4/7|

Hence

    sum B_k = 1/70+1/14 = 3/35,

exactly reproducing the all-permutation simultaneous-birth atom.

### Triangular self-matching L=3

| k | e0_k | e2_k | B_k | A_k | C_k |
|---:|---:|---:|---:|---:|---:|
|2|3/28|0|0|3/28|0|
|3|15/28|0|0|15/28|0|
|4|5/14|5/14|2/35|3/10|3/10|
|5|0|15/28|0|0|15/28|
|6|0|3/28|0|0|3/28|

Thus

    sum B_k=2/35.

Self-matching symmetry makes the nonzero direct-jump flux sit exactly at the middle cardinality k=4 for N=9; square-site asymmetry spreads it over k=4,5.

## 4. The full one-step transition kernel once B_k is known

For q_j(k)>0,

    P(R_{k+1}=2 | R_k=0) = B_k/q_0(k),
    P(R_{k+1}=1 | R_k=0) = [e0_k-B_k]/q_0(k),
    P(R_{k+1}=2 | R_k=1) = [e2_k-B_k]/q_1(k),           (9)

with the remaining mass staying in the same rank.

These are exact **one-step conditional probabilities**. They do not make `R_k` a Markov sufficient state for future persistence: the law of the remaining waiting time inside rank one can still depend on hidden geometry. Higher persistence moments require the multi-time chain counts described in the copula note.

## 5. A hierarchy of missing process information

The static rank-sector polynomials determine:

    q_j(k), J1 marginal, J2 marginal, E D.

Adding the one-step jump-two flux sequence `{B_k}` determines:

    the complete one-step rank transition table,
    P(D=0),
    cardinality location of simultaneous births.

It still does not determine:

    P(D=d), d>=1,
    E[D(D-1)],
    the J1/J2 copula.

Those begin at two-time nested rank-one chain counts. Hence there is a clean information ladder:

    one-time C[j,k]
      -> + jump-two flux B_k
      -> + two-time rank-one chain kernel
      -> + higher chain counts / full copula.             (10)

This is a more economical process representation than storing every permutation trajectory.

## 6. Interfaces

**#769.** Report `B_k` (or equivalently the fixed-site jump-two cardinality histogram) in addition to p-weighted jump-two probabilities. It is a process observable with exact combinatorial semantics.

**#778.** If paired raw data are missing but #769 supplies B_k, the simultaneous-birth atom and full one-step rank kernel are still recoverable exactly. Do not confuse this with recovery of the full copula.

**#768.** The scaling of the cardinality-resolved B_k through the critical window directly tests whether the abundant theta/T3 jump-two spines have 6-arm or 8-arm cost before any matching-odd one-point cancellation.

## 7. Boundaries

- The process rank alone is not asserted Markov.
- The L=3 tables are exact controls, not scaling evidence.
- B_k measures direct rank-0 to rank-2 insertion, not by itself a continuum arm count.
