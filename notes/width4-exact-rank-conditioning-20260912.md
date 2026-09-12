# Rejection-free rank-conditioned configurations from exact backward messages

2026-09-12. A concrete sampler and conditional score on the certified width-four
rank automaton. This is the standard finite-horizon conditional change of
measure, specialized and verified here; not a new general Doob-transform theorem.
The implementation uses the source-safe 509 states, not a homogeneous quotient
for an inhomogeneous input.

## 1. Integer weights and the backward recursion

For each site use a rational p_(y,j)=a_(y,j)/d_(y,j). A row-mask weight can be
written w_y(b)=W_y(b)/D_y, where

    W_y(b)=product_j [a_(y,j) if b_j=1 else d_(y,j)-a_(y,j)],
    D_y=product_j d_(y,j).

Let c_r(s)=1{close_rank(s)=r} be the final rank event. Define h_m=c_r and,
for y=m-1,...,1,

    h_y(s)=sum_b W_y(b) h_(y+1)(tau(s,b)).

The exact unnormalized rank probability is

    Z_r=sum_b W_0(b) h_1(initial(b)),
    P_r=Z_r/product_y D_y.

If Z_r=0 the conditioning is impossible and the program raises an error. Zeros
at individual states are not inverted; the sampler can only reach states
with positive conditional mass.

The first row is drawn with integer weights W_0(b)h_1(initial(b)). After that,
at a state s before row y, the next mask is drawn with integer weights
W_y(b)h_(y+1)(tau(s,b)). Exact uniform integer selection avoids rounded
transition probabilities, even when P_r is tiny.

Multiplying these conditional probabilities along a completed word makes all
h factors telescope, giving

    product_y W_y(b_y) / Z_r   if the final rank is r,
    zero otherwise.

That is exactly the original site product measure conditioned on rank r.
The construction neither repeats percolation configurations until one passes
nor uses burn-in, an asymptotic cylinder eigenvector, or an approximate tilt.
It does require the exact backward normalizers. Their construction is the
nontrivial resource, not free conditional information.

The state count is fixed at 509 here. The elementary arithmetic operation
count is O(509*16*m) for a backward table and O(16*m) integer sampling decisions
per word. Integer bit lengths grow with m and the denominators. These are NOT
unit-cost timing bounds or complexity claims as circumference increases.

## 2. Conditional occupation moments and a thermal likelihood score

Let K be the total occupied site count (not a rank-birth time). A second forward
pass propagates the exact weight, K-weight and K^2-weight by

    (z,z1,z2) -> W_b*(z, z1+kz, z2+2kz1+k^2z).

Summing by final rank yields exact conditional means and variances.
For homogeneous interior p, direct differentiation of the product measure gives

    d_p log P_r = [E(K | r)-Np]/[p(1-p)],
    d_p log(P_2/P_0) = [E(K|r=2)-E(K|r=0)]/[p(1-p)].

This is also an exact normalizer/thermal-jet map. It does not estimate a ratio
from an unobserved rare denominator. The rank-2 event is increasing and rank-0
event decreasing, so their conditional mean difference is strictly positive
in the nondegenerate interior setting; the finite polynomial check confirms it.
For general inhomogeneous fields the corresponding per-site conditional score
is obtained from its own occupied indicator, not from a homogeneous K formula.

## 3. Executed controls, not a new probability production

On a nonuniform 4x3 product measure with twelve specified rational site
probabilities, all 4,096 configurations were physically classified. Their
rank normalizers and conditional first/second K moments agree with the forward
DP. For all three target ranks, all 12,288 sequential word-probability identities
agree exactly with the product-weight conditional formula. This validates the
ENTIRE finite conditional distribution at the control size, not just empirical
frequency agreement from a few draws.

A long control uses p=591417/1000000 and m=128, N=512. This is a rational near
q4, not the exact cylinder root and not an estimate of two-dimensional p_c.
Forward and backward integer normalizers agree:

| target rank | exact-probability diagnostic | conditional mean K | conditional variance K |
|---|---:|---:|---:|
| 0 | 4.2444165285438005e-10 | 265.21829715618895 | 78.51788627252415 |
| 2 | 4.2441877885758733e-10 | 341.4414449341080 | 73.32447740787667 |

Their total mass is 8.488604317119674e-10, and the thermal log-odds derivative
is 315.4371084311523. Exact rational quantities, rather than just these
floating diagnostics, are retained. The large derivative rational is stored
with hexadecimal numerator/denominator to avoid Python's decimal-conversion
size limit; this is lossless, not a numerical approximation.

Four fixed PRNG seeds were used per rare rank as algorithm controls. Every
resulting 512-site configuration was reclassified by independent physical
lift traversal and has the required rank. Its full conditional path probability
was independently multiplied and checked against the telescoping formula.
The eight row-mask words and seeds are retained. They are NOT new evidence
about sector probabilities, covariance, or finite-size scaling. The law is
exact given uniform integer bits; seeded pseudorandom bits serve reproducibility.

## 4. Scientific use and limits

The prior observation that direct one-p snapshots can miss both rare sectors
remains correct. It is not an unavoidable cost of the underlying physical
question: a finite exact sufficient-state representation supports a different
conditional algorithm. This example makes that limitation concrete.

This routine is useful as an exact oracle for testing future importance or
conditional algorithms and for generating configurations with a declared
rank. It does not solve a large-width sampling problem, lower-bound the cost
of other estimators, beat an optimized existing transfer implementation, or
justify a new production. No parameter was tuned to obtain a favourable rank
verdict. Source-safe states and normalizers are mandatory.

Background: Corstanje--van der Meulen--Schauer, arXiv:2111.11377v2, §2,
Definitions 2.3 and Example 2.4, describes the established h-transform
conditioning principle for Markov processes. The discrete finite recursion
above is proved directly by its path probabilities; the continuous-time
paper is attribution/background, not a claim that its generator formula is
our discrete update. Primary §2 read at
https://arxiv.org/html/2111.11377v2 . No literature novelty claim is made.
