# Finite-window geometric susceptibility from the original uniform source

2026-09-30. New finite response identity; the consumer of the
[clock-preserving geometric intervention](birth-selection-intervention-20260930.md).
This changes the mechanism question, not the precision target of the old
birth-history contrast.

## Identity

Fix a rank-one prefix A at count b. Let m=N-b, c its completing vacancies,
s=m-c, and d_v the synergy degree of a safe vacancy. Set mean_d=2e/s.
Apply the previously specified source kick once, then return to the uniform
source. Let S_h(theta|A) be survival through h insertions, h>=1.

At theta=0 the first-site likelihood score is d_v-mean_d for safe v and
zero for completing v. Given the unordered set B of the next h sites, all
orders are equiprobable. Survival depends only on A union B. If it survives,
every member of B was safe at A; otherwise its contribution is zero.
Consequently

    S_h'(0|A) = E_uniform [1{rank(A union B)=1}
                          * (sum_{v in B} d_v/h - 2e/s)].       (1)

This is an exact likelihood derivative, with an exact conditional average
over the h possible first sites. It does not approximate a small nonzero
theta by finite differences or hold geometry fixed during continuation.
Degrees are those at the starting prefix, not degrees recomputed later.
For s=0 survival and response are zero. Assume 1<=h<=m.

Proof at the subset level: the perturbed probability of B relative to its
uniform value is (m/h)*sum_{v in B} P_theta(first=v|A). Differentiate and
multiply by the endpoint survival indicator. The endpoint only sees the
set, so internal-order averaging leaves the expectation unchanged. For
h=1 the safe scores sum to zero. For h=2 it reproduces

    S_2'(0|A) = -s Var_safe(d)/[m(m-1)],

but h=735 uses its actual continuation, not that expression extrapolated.
No assumption about the sign for arbitrary h is made.

## Concrete consumer and decision

The input is the existing square-L512 safe-insertion block: 14x10000 full
Fisher-Yates permutations, b=155385. Replaying the same seeds also recovers
their unused continuation. Use **h=735**, the unchanged destination156120
from the original completion-hazard contract, with no lag or subgroup scan.

Primary: average (1) over the original rank-one risk population at b.
Accompanying quantities are its ordinary endpoint survival and the exact
two-step susceptibility. Save their complete aligned 14-batch covariance.
The relevant question is whether a kick invisible to the occupation clock
and immediate exit can influence this specified much longer rank window.

- A resolved response supplies a finite-size transmission beyond local
  curvature, for this explicit source and observer.
- A small/unresolved response despite a nonzero local susceptibility leaves
  rapid loss, cancellation and insufficient precision as distinct possibilities;
  it is not an invitation to search a more significant lag.
- Neither outcome identifies original U, continuum scaling or the fraction
  of passive early/late association mediated by this intervention.

These prefixes have already informed the project; the new readout is not
independent evidence for the old history effect. The other worker's
independent-archive replication is a separate question and remains separate.

[Contract](../analysis/geometric-source-window-20260930/contract.json) and
[engine](../analysis/geometric-source-window-20260930/engine.cpp) are fixed
before this response is computed. The old topology/pair algorithms are
reused; the only change to their source is an entry-point guard for inclusion.
Medium-size replay goes to Huawei, not this Mac.
