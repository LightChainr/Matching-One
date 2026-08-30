# Transport bound for the remaining fixed-line ULC gap

## Radon--Nikodym decomposition

For uniform layer measure `mu` and the marginal of a uniform internal edge,

`d nu_k^up/d mu_k = u/E[u] = (1-h_x)/(1-xi_k)`,

`d nu_(k+1)^down/d mu_(k+1) = d/E[d] = (1-h_beta)/(1-beta_(k+1))`.

Consequently the degree-bias term is not opaque:

`xi_(k+1)-xi_k = Delta_edge - Var_k(h_x)/(1-xi_k) + Cov_(k+1)(h_beta,h_x)/(1-beta_(k+1)).`

This proves a conditional variance--association lemma: uniform exit hazard increases whenever edge slack dominates the lower-layer variance penalty and the upper-layer birth/exit covariance is nonnegative.

## Exact bounded geometry

All 984 audited carrier-layer pairs satisfy both conditional hypotheses. Among the 484 negative-bias pairs, the maximum exact ratio `(-bias)/Delta_edge` is 5/14, achieved by 8 N=12 quotient/line realizations.
There are 78 uniform-hazard equalities; every one is trivial in the transport decomposition (`Delta_edge=bias=0`). There is no nontrivial saturation of the candidate inequality.

## Stronger transports that fail

First-order stochastic dominance fails minimally at N=12, `[[3, 0], [0, 4]]`, carrier `matching`, line `[0, 1]`, layer 4. At threshold `3/4`, the lower uniform tail is `4/19` while the upper edge tail is only `1/5`; the mean inequality nevertheless survives.
Pointwise birth/exit comonotonicity already fails at N=8, `[[2, 0], [0, 4]]`, layer 5: masks 87 and 91 have degree pairs `[0, 2]` and `[2, 1]`. The aggregate covariance remains positive (1/90).
A sign-free Cauchy bound is also too strong: it passes 980/984 pairs but fails on four N=12 realizations.

## Status

- **Proved:** the Radon--Nikodym and variance/covariance decomposition, and the two-hypothesis conditional lemma.
- **Exact finite evidence:** variance domination and nonnegative aggregate birth/exit association hold on all 984 existing pairs, with substantial slack.
- **Disproved as proof routes:** first-order stochastic dominance, pointwise comonotonicity, and worst-case Cauchy control.
- **Still open:** prove the two aggregate hypotheses for arbitrary quotients by a two-step path sum; the present result does not claim general ULC.
