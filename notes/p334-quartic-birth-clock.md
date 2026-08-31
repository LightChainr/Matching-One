# P334: the physical fourth-order increment keeps eroding the pair-clock gap

## New result

All minimal quartic triggers of the two real N425 same-scalar-state witnesses
are now exported: **1,178 in A and 2,866 in B**. The complete fourth-order
physical safety coefficients and the full `<=4` truncated clocks are exact.

| Quantity | A: counter 43042514269 | B: counter 43042505280 |
| --- | ---: | ---: |
| Pair+triple-safe four-sets | 34,536,807 | 34,632,298 |
| New minimal quartics | 1,178 | 2,866 |
| **True rank-one-safe four-sets** | **34,535,629** | **34,629,432** |
| True fourth-step survival | 0.9582063182 | 0.9608089239 |
| True fourth-step conditional hazard | 0.0208670039 | 0.0189333806 |
| Mean first hit, pair-only | 21.60624182 | 31.75099733 |
| Mean first hit, pairs+triples | 18.00285614 | 23.11369450 |
| Mean first hit, pairs+triples+quartics | 17.75453001 | 21.17202754 |
| Additional mean shortening from quartics | 0.24832613 | 1.94166696 |

The mean waiting-time gap evolves as

\[
 10.14475552\quad\longrightarrow\quad5.11083837\quad\longrightarrow\quad3.41749753.
\]

Quartics remove a further **33.13% of the preceding gap**; only 33.69% of the
original pair-only gap remains. Their survival loss is greater in B at every
step `k=4,...,162`. Nevertheless **B remains later than A throughout the
complete `<=4` curve**, for both survival and all commonly defined hazards.
This is evidence of a material, asymmetric hierarchy of collective triggers,
not a reason to extrapolate the pair-only 10-step gap to physical birth.

## The new physical coefficient

All 173 vacant sites are singleton-safe. The existing complete pair/triple
lists filter four-sets containing a lower-order trigger. Only the remaining
34.54/34.63 million sets are sent to the existing potential-union-find rank
oracle, reusing the one-, two-, and three-insertion snapshots in the nested
enumeration. A rank-two result is therefore a genuine minimal quartic.

Two explicit physical examples are:

- A: site labels `{1,50,93,361}`;
- B: site labels `{8,22,23,24}`.

Their 15 proper subsets have ambient rank one, while each full four-set has
rank two. The complete lists and first-witness subset ranks are saved in
`results/p334-quartic-clock/full_quartics.json`, together with the original
occupied prefixes and seed/counter identity.

The fourth-step differences are exact:

\[
 S_{B}(4)-S_{A}(4)=\frac{4937}{1896945}=0.0026026057687\ldots,
\]

\[
 h_{B}(4)-h_{A}(4)=
 -\frac{37762802238}{19529555735825}=-0.0019336232093\ldots.
\]

These are **full physical birth** statements through four insertions, not
merely surrogate-graph counts. The higher-k statements below are truncated.

## Complete fourth-order hypergraph clock

For the uniform without-replacement continuation from a fixed checkpoint,
let `i_k` count k-sets containing no member of `E2 union E3 union E4`. Then

\[
 I_{\le4}(z)=\sum_k i_kz^k,\qquad
 S_{\le4}(k)=i_k/\binom{173}{k},\qquad
 h_{\le4}(k)=1-\frac{k i_k}{(174-k)i_{k-1}}.
\]

The earlier weighted-link deletion/contraction calculation now includes
quartics. The real hypergraphs have 78/113 involved sites, reduced to 43/70
equal-link groups. Only 1,059/1,536 memoized states are needed to obtain all
coefficients; the polynomial calculations take approximately 0.03/0.04 s.

The dominance certificate is

\[
 I_{\le4,B}(z)-I_{\le4,A}(z)=z^3(1+z)^{74}Q_{80}(z),
\]

with **all 81 coefficients of Q positive**. The first are
`614,48367,1882881,48200814,911481167`; the leading coefficient is 6. Full
coefficients and exact rational curves are in `quartic_survival.json`.

Consequently B's survival is strictly greater at `k=3,...,157`; both are zero
from step 158. Integer hazard cross-products show B's hazard is strictly
lower at every common step `k=3,...,155`. A's maximum `<=4`-safe subset is
154, versus 157 for B; their median/90%-absorption steps are 16/31 and 19/39.

## Scope and execution

Same N425 second orientation, period `[[425,268],[0,1]]`, seed
`20260831430425`, k0=252, age=10, ambient line `(12,-19)`, H2=0, b2=14770.
Parent `d5d2cc89e77ebb2ec6252df75dc858e9c240e6ce`. These are the same selected
two witnesses, not new random draws or independent confirmation.

The entire rank census was **one sequential process**, 18.59 s for A and
20.95 s for B on the local Mac. No server connection, new MC, broad search,
or full-suite rerun occurred. The new output is the complete previously
unknown fourth-order layer, not a revalidation of the first three orders.

For full physical birth, only the following bound follows at later times:

\[
 T_{\mathrm{true}}\le T_{\le4}\le T_{\le3}\le T_{\mathrm{pair}},\qquad
 S_{\mathrm{true}}\le S_{\le4}\le S_{\le3}\le S_{\mathrm{pair}}.
\]

Minimal triggers of size five or more are not counted. Thus neither full
physical A/B clock ordering beyond k=4 nor convergence of the shrinking mean
gap is proved. There is no population, field-identity, or exponent claim.

```sh
/Users/lc/python-envs/research-py311/bin/python scripts/p334_quartic_constraints.py
clang++ -std=c++17 -O3 -pthread src/p334_checkpoint_safe_quartics.cpp -o /tmp/p334-full-quartics
/tmp/p334-full-quartics results/p334-quartic-clock/frozen_constraints.txt results/p334-quartic-clock/full_quartics.json
/Users/lc/python-envs/research-py311/bin/python scripts/p334_quartic_survival.py
```
