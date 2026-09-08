# A marked birth law gives the complete single-site knockout response

The exact `(birth time, final insertion site)` distribution is more than a
localization map. It determines the entire response to permanently blocking
one site, **without another network solve**. The response is not proportional
to a site's winning probability alone: the birth time and the thermal kernel
give distinct weights to the same local mark.

Applied to the original two N425 witnesses, the saved marked archive already
determines all 346 one-site interventions. This is a new interpretation and
readout of those exact laws, not a new random block or 346 new simulations.

## The inert-clock convention and exact identity

Fix a monotone birth event on d selectable sites with uniform random insertion
order. Let T be its first arrival and V its final insertion site. For one
site v, F0(z) counts safe subsets of the other sites with v absent; F1(z)
counts safe subsets when v is occupied. The original and blocked polynomials
are

`F(z)=F0(z)+zF1(z)`,
`F_inert_v(z)=(1+z)F0(z)`.

The blocked site remains an **inert dummy among the original d insertion
positions**. Deleting it from the permutation instead would define a different
clock. With D_v=F0-F1, the polynomial change is exactly z D_v. The coefficient
D_v(k-1) counts safe (k-1)-subsets which would become unsafe upon insertion of v.
Uniform ordering therefore gives

\[
P(T=k,V=v)=\frac{D_v(k-1)}{d\binom{d-1}{k-1}}
=\frac{D_v(k-1)}{k\binom d k},
\]

and hence the complete fluctuation--response relation

\[
\boxed{S_{\rm inert\ v}(k)-S(k)=kP(T=k,V=v).}
\tag{1}
\]

If eventual birth remains possible after blocking v, summing (1) yields

\[
\boxed{E[T_{\rm inert\ v}]-E[T]=E[T\,1_{\{V=v\}}]
=\pi_v E[T\mid V=v].}
\tag{2}
\]

If blocking makes birth impossible, the true mean is infinite; the same sum
including k=d instead describes a clock censored at d+1. All individual
knockouts in the two current archives retain eventual birth. The saved
last coefficient makes this distinction explicit.

The identity needs monotonicity and uniform insertion order, not the torus
embedding, a low treewidth, a pair-only event, or a new intervention sample.
Summing over v gives another exact consequence:

`sum_v Delta E[T]_v = E[T]`,
`sum_v Delta S_v(k) = k P(T=k)`.

Thus normalized one-site mean impacts form the **time-size-biased final-site
distribution**. This is a budget across separate singleton interventions;
it is not an additive formula for blocking all sites simultaneously.

## A direct local-to-thermal response formula

In a fixed N-site prefix at k0, K2=k0+T and
`F2(p)=E[Pr(Bin(N,p)>=k0+T)]`. Applying (1) to its Bernstein coefficients gives

\[
\boxed{F_{2,\rm inert\ v}(p)-F_2(p)
=-E[T\,1_{\{V=v\}} B_{N,k0+T}(p)].}
\tag{3}
\]

Here B_(N,j) is the binomial **probability mass**, not the tail. Integration
immediately gives

`integral Delta F2(p) dp = -Delta E[T]/(N+1)`.

Equation (3) separates three measures of a site's role: how often it wins,
how much blocking it changes the whole waiting clock, and how much it changes
the near-reference thermal observable. Since the prefix and its earlier
history are unchanged, the isolated second-birth contribution loads equally
into A and E; this does not invent a new independent A/E coordinate or a
cross-orientation covariance.

## Actual knockout landscape from the old exact marks

Inputs are pinned to `1c06230b8f7e13be98f128361ad72b23c0c425ae`, counters
A=43042514269 and B=43042505280, each N425/k0=252/d173. These are **not** the
new 8631/14803 triple-bridge pair; their analyses are complementary.
All probability/clock increments below use exact integer/rational marked
counts. Canonical values evaluate the exact weights against binomial masses
at `p_ref=.59274605079`.

| Individual permanent knockout | Original winning probability | Mean-wait increase | S(40) increase | F2(p_ref) change |
|---|---:|---:|---:|---:|
| A, site6 | .04450993 | .7719247 | .00767269 | -.00637315 |
| B, site121 | .09752434 | 1.8167002 | .02430951 | -.01300597 |
| B, site8 | .09392448 | 1.7825427 | .02430424 | -.01233000 |
| B, interior144 | .02437087 | .6401188 | .01353032 | -.00198895 |
| B, interior413 | .01770113 | .4890574 | .01105231 | -.00127389 |

A's eight equal strongest port sites are 6,27,140,162,251,274,296,409; listing
site6 is not a uniqueness claim. B's strongest single knockout is site121.
The top-five site sets in the saved deterministic tie convention do not
change between winning frequency, mean impact, step40 impact and the reference
thermal loss. No rank reversal is manufactured where the data give none.

The *relative roles* nevertheless differ materially:

| Interior-site share of | A | B |
|---|---:|---:|
| Final-birth probability | 12.4462% | 23.6534% |
| Total individual mean-wait impact | **14.6034%** | **29.7108%** |
| Total individual F2(p_ref) loss | **10.6203%** | **16.6261%** |

In B, interior sites account for almost 30% of the complete-clock impact but
only 16.6% of the near-reference thermal impact. The latter uses a binomial
mass centered near the checkpoint and emphasizes a different part of the
arrival-time distribution. Thus a local structure's importance for integrated
clock delay cannot be identified with its thermal-observable importance.

The inverse-Simpson effective site counts for mean impact are 38.5867 and
23.8353, compared with the previously reported winning-probability counts
36.4300 and21.9519. B has more sites with nonzero participation, yet its
one-site delay sensitivity remains more concentrated. This is an exact
property of these two conditional networks, not an ensemble robustness law.

## Several bridge knockouts are not a sum of singleton effects

For two sites v,w, let F_ab condition their occupancies. The interaction is

`F_inert_vw - F_inert_v - F_inert_w + F
 = z^2 (F00-F10-F01+F11)`.

Its coefficients need not have one sign: serial and alternative routes can
produce opposite effects. The new marked-triple lemma says each genuine
minimal triple has exactly one internal middle. For a set of those middle
sites, the knockout changes are consequently additive through degree3, but
nonadditivity may start at degree4. The parallel 8631/14803 intervention uses
this fact; it does not assume the full clock is additive in the middle marks.

## Scientific card

- Change: a complete finite single-site intervention response, including
  canonical second-birth loading, is recoverable from the original marked
  event law. No counterfactual network replay is required when that law exists.
- Source/dependency: the same two fixed N425 prefix laws in1c06230b; no new
  random samples and no whole-population or continuum-operator assertion.
- New observation: interior structure has substantially different weights
  in time-integrated and near-reference responses, even with the same final
  birth marks and unchanged leading sites.
- Output: `results/p334-clock-knockout-from-marked-law/score.json`, exact
  per-site clock/tail/integral changes, canonical changes, type shares and
  deterministic orderings. Original D_v arrays are referenced, not duplicated.
- Next use: record the final site of a prefix-defined intervention mark in a
  production stream and use (1)--(3). A mark selected from the future suffix
  instead would change the estimand; a simultaneous multi-site intervention
  requires its higher-order interaction data.

```sh
/Users/lc/python-envs/research-py311/bin/python scripts/p334_clock_knockout_from_marked_law.py
```
