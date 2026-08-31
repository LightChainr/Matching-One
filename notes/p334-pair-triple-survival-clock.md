# P334: genuine triples cut the long tail without reversing the pair-clock order

## Result

The complete pair-plus-triple clocks of the two saved real N425 checkpoints
are now exact. **B remains later than A at every comparable step**, but its
mean delay advantage falls from **10.14476 to 5.11084 insertions**. This is not
a small correction: genuine triples remove 49.62% of the pair-only mean gap.

| Fixed real checkpoint | A: counter 43042514269 | B: counter 43042505280 |
| --- | ---: | ---: |
| Minimal trigger pairs | 108 | 108 |
| Genuine minimal triples | 583 | 509 |
| Pair-only mean first-trigger step | 21.60624182 | 31.75099733 |
| Pair+triple mean first-trigger step | 18.00285614 | 23.11369450 |
| Mean clock shortening by triples | 3.60338568 | 8.63730283 |
| Pair-only median / 90% absorption step | 19 / 40 | 26 / 65 |
| Pair+triple median / 90% absorption step | 16 / 32 | 20 / 44 |
| Largest pair-only safe subset | 162 | 168 |
| Largest pair+triple safe subset | 155 | 162 |

The mechanistic finding is that **fewer minimal triples do not imply a weaker
finite-time effect**. B has 74 fewer triples, yet their inclusion shortens its
mean clock 2.397 times as much. Their contribution is concentrated in the
long pair-survival tail, not captured by the third-order count alone.

More precisely, put
`L_X(k) = S_pair,X(k) - S_pair+triple,X(k)`. Exact integer comparisons give

- `L_A(k) > L_B(k)` for `3 <= k <= 18`;
- `L_B(k) > L_A(k)` for `19 <= k <= 168`.

Thus there **is a crossing in the triple-induced survival loss at steps
18/19**, but there is **no crossing in the resulting A/B survival or hazard**.
The tail effect removes 5.03391715 steps from B's relative advantage.

## Same real starting state; no new Monte Carlo

Both checkpoints are in the same N425 second-orientation quotient,
`period_matrix=[[425,268],[0,1]]`, with seed `20260831430425`, `k0=252`,
age 10, ambient line `(12,-19)`, and `d=b1=173`, `H2=0`, `b2=14770`.
They are the previous same-scalar-state witnesses, not new graphs or draws.

The old `safe_triple_census.json` retained only the first eight minimal
triples. This increment replays **only those two occupied prefixes** through
the existing potential-union-find engine to export all 583/509 triples.
The new `full_triples.json` includes the complete occupied prefixes, counters,
full triple lists, and the existing first-triple subset-rank certificate.
The original production/archive files are untouched.

We then consider a fresh uniform, without-replacement order of the 173
vacant sites, conditional on the fixed saved checkpoint. This is a finite
conditional combinatorial distribution, not a population estimate or an
additional sampled continuation of the original seed.

## Exact truncated clock

Let `E2` be all minimal trigger pairs and `E3` all triples whose singleton
and pair subsets remain rank one but whose full insertion reaches rank two.
Define the independence polynomial

\[
 I_{\le3}(z)=\sum_{U\subseteq V:\ \nexists e\in E_2\cup E_3,\ e\subseteq U}
 z^{|U|}=\sum_{k=0}^{173}i_kz^k.
\]

The complete survival and conditional hazard are

\[
 S_{\le3}(k)=\frac{i_k}{\binom{173}{k}},\qquad
 h_{\le3}(k)=1-\frac{k\,i_k}{(174-k)i_{k-1}},
 \qquad E[T_{\le3}]=\sum_{k=0}^{172}S_{\le3}(k).
\]

All integer coefficients, exact rational probabilities, and hazard
cross-products are in `pair_triple_survival.json`. The first coefficients are:

| k | A: pair+triple safe k-sets | B: pair+triple safe k-sets |
| --- | ---: | ---: |
| 0 | 1 | 1 |
| 1 | 173 | 173 |
| 2 | 14770 | 14770 |
| 3 | 829921 | 830535 |
| 4 | 34536807 | 34632298 |
| 5 | 1135712900 | 1143077173 |
| 6 | 30750255347 | 31125789469 |
| 7 | 705306504518 | 719552995192 |

For orientation on the full curve:

| k | A survival | B survival | A hazard | B hazard |
| --- | ---: | ---: | ---: | ---: |
| 10 | 0.74612109 | 0.79062351 | 0.05347783 | 0.03905919 |
| 20 | 0.34412207 | 0.48457894 | 0.08892896 | 0.05322617 |
| 40 | 0.03172400 | 0.12977562 | 0.12986032 | 0.07336329 |
| 65 | 0.00054767 | 0.01339199 | 0.16829452 | 0.09953307 |

The exact mean difference is

\[
 E[T_{\le3,B}]-E[T_{\le3,A}]
 =\frac{6012045525603544194713991}{1176332549747380259634400}
 =5.110838365302944\ldots.
\]

## Small exact computation, not a large engine

The pair+triple hypergraphs involve only 67/79 nonisolated sites. Group sites
with identical hypergraph links. Such sites never occur together in one
minimal edge: if `u,v` did, the link of `u` would contain `v`, but the link of
`v` never contains `v`. Therefore any nonempty subset of a group of size `s`
imposes the same external constraints and carries weight `(1+z)^s-1`.

This reduces the two *actual* hypergraphs to 30/33 weighted variables.
Weighted deletion/contraction, forced exclusion for singleton constraints,
and component factorization give the full polynomials in 104/130 memoized
states. No new Monte Carlo, expanded quotient search, or full-suite rerun
was performed. The only prefix replay obtains triple labels absent from the
old artifact; it does not repeat the old census as a validation campaign.

There is a short coefficientwise dominance certificate:

\[
 I_{\le3,B}(z)-I_{\le3,A}(z)
 =z^3(1+z)^{105}Q_{54}(z),
\]

where **all 55 coefficients of `Q54` are positive**. They are saved under
`comparison.dominance_factor`; the first are
`614,31021,754628,11823082,134540304`, and the leading coefficient is 2.
Consequently `S_B>S_A` exactly at `k=3,...,162`; they agree at `0,1,2`
and are both zero at `163,...,173`.

For every step with both prior survivals positive, the saved integer
cross-products

\[
 i_{A,k}i_{B,k-1}-i_{B,k}i_{A,k-1}
\]

are negative at `k=3,...,156` and zero at `k=1,2`. Hence B's conditional
hazard is also strictly lower throughout the common domain. A's hazard is
undefined from step 157 onward because A is already surely absorbed.

## Physical boundary and next useful observable

This is the **pair+triple truncation**, not the full ambient-rank-two clock.
It is physically exact through three insertions. At four or more insertions,
unexported minimal triggers of size at least four can cause earlier birth:

\[
 T_{\mathrm{true}}\le T_{\le3}\le T_{\mathrm{pair}},\qquad
 S_{\mathrm{true}}\le S_{\le3}\le S_{\mathrm{pair}}.
\]

These inequalities do **not** give a pointwise hazard bound or transfer the
A/B ordering to full physical birth. Nor do two selected witnesses establish
a population law or continuum exponent.

What is established is a concrete collective effect missed by ranking
checkpoints using the bare minimal-triple count:
**the placement and overlap of minimal triples relative to the pair-safe
long tail control their finite-time effect**. The minimal next physical
increment is an exact four-insertion boundary on these same two prefixes,
or extraction of the minimal fourth-order triggers, not another estimate
of the already known third-order counts.

## Reproduce

Parent: `ad6c595a70c66ea4421c816b4c65b1cfe3d9c803`.

```sh
clang++ -std=c++17 -O3 -pthread src/p334_checkpoint_safe_triples.cpp -o /tmp/p334-full-triples
/tmp/p334-full-triples results/p334-pair-triple-clock/full_triples.json
/Users/lc/python-envs/research-py311/bin/python scripts/p334_pair_triple_survival.py
```

The existing local research Python environment is reused; no dependency
installation or remote host is involved.
