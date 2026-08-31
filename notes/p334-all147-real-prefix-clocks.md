# P334: complete physical birth clocks on all 147 fixed-source prefixes

## Result

The full conditional rank-one-to-rank-two birth clock is now solved for **all
147 eligible real N425 prefixes**, not only the previous twelve. No new Monte
Carlo was sampled. The 135 new evaluations took 14.46 seconds in total on one
local research-Python process; the old twelve clocks were reused unchanged.

The strongest new mechanism witness is more specific than a broad spread of
waiting times: **two real prefixes agree in singleton survival, pair survival,
and the exact safe-successor degree second moment, but contain respectively
5 and 19 genuine minimal triple triggers.** Their physical clocks consequently
already differ at the third insertion. Second-moment branching information
does not close the collective birth state.

## Fixed source and exact observable

Selection was committed in `8d7ac0e9` before evaluating the new rows:
`analysis/p334_all147_clock_manifest.json`. It retains every source row with
N425, second orientation, k0=252, age=10, ell=(12,-19), and rank one, excluding
the original two case-study counters 43042514269 and 43042505280. Rows are
ordered by original counter, not by their results. The original source is
`results/local-20260831/P334-cooperative-closure/raw/N425.geometry_pilot.csv`;
the selected raw rows and source SHA256 are saved in `selection.json`.

There are d=173 initially vacant sites. For each fixed occupied prefix, let
f_k count the k-subsets whose insertion leaves ambient black-NN rank one. Then

\[
S(k)=P(T>k\mid\text{prefix})=\frac{f_k}{\binom{173}{k}},\qquad
E[T\mid\text{prefix}]=\sum_{k=0}^{172}S(k).
\]

These are **full physical** safe-subset counts, not pair-only or truncated
hypergraph counts. The saved occupied-prefix replay, contracted gains, gauge,
port addresses, and balanced remainder certify each finite mapping. Its
site-disjoint two-port factors supply safety polynomials F_i; off-core sites
supply `(1+z)^r`. Multiplying them gives the full F. We reused the existing
evaluator without changing it to improve pass rate.

All 147 completed under the frozen 20-second/200,000-state per-prefix caps;
none were replaced. Their 535 factors range from 1 to 14 per prefix. Maximum
factor treewidth was 9 and maximum DP boundary-state count 25,882; the slowest
new prefix, counter 43042508795, took 5.153 seconds. This establishes tractable
conditional solves for this fixed stratum, not an all-geometries cost bound.

## Broad physical-clock variation survives the instantaneous gate

| Fixed-source readout | Exact-polynomial result, decimals for display |
| --- | ---: |
| Completed / selected | 147 / 147 |
| Conditional mean waiting time, min / median / max | 5.77269 / 13.73653 / 30.86927 |
| S(40), min / max | 0.00006809 / 0.24899546 |
| All unordered prefix pairs | 10,731 |
| Pairs with a strict survival-order crossing | 3,510 |
| Crossing pairs with first sign switch by step 40 | 1,867 |
| Repeated (H2,b2) keys | 8, all with different full clocks |
| Such differing pairs entirely among the new 135 | 7 |

Crossings use integer signs of `f_B(k)-f_A(k)`, so no numerical threshold is
involved. `scientific_summary.json` saves sign runs and switch certificates,
not 10,731 duplicated difference polynomials. Each complete F is retained and
reconstructs every comparison. The counts describe this selected source
stratum; they are not estimates of population-wide crossing frequencies.

Within H2=0 alone (16 prefixes), mean waiting ranges from 17.55046 to 30.86927
and S(40) from 0.0154376 to 0.2489955. H2=5 (11 prefixes) still spans means
10.66158 to 20.01983. Thus equal immediate birth hazard does not approximately
fix the remaining clock across these real prefixes.

The lexicographically first new/new repeated-(H2,b2) pair is
43042505485 / 43042509612: H2=0, b2=14,806, but means 19.80945 / 22.15515 and
S(40)=0.0379825 / 0.0765829. Its 2.34570-step mean gap is entirely collective,
since neither prefix has an original single-site trigger.

## A stronger exact state collision: genuine triples 5 versus 19

The fixed147 search has exactly one different-clock pair sharing the larger
key `(H2,b2,checkpoint_sum_child_b1_sq)`:

| Same geometry, age and winding line | 43042508631 | 43042514803 |
| --- | ---: | ---: |
| H2 / safe singleton count a | 15 / 158 | 15 / 158 |
| Safe pairs b2 | 12,397 | 12,397 |
| Sum of squared safe-successor counts | 3,890,796 | 3,890,796 |
| True safe triples f3 | 644,020 | 644,006 |
| Genuine minimal triple triggers g3 | **5** | **19** |
| Mean true birth step | 10.48090361 | 10.39836128 |
| S(40) | 0.01078726 | 0.00935215 |
| Original-direct-site winning fraction | 0.95958062 | 0.95094000 |

The source implementation explicitly sums the square of the number of safe
second insertions over each safe first insertion
(`src/threshold_rank_integer_period_mc.cpp`, cooperative-closure block).
Consequently this is the safe-pair graph degree-square sum, not an unrelated
score moment. If the complementary minimal-pair-trigger graph has m edges,

\[
m=\binom{158}{2}-12397=6,\qquad
\sum_v d_{\rm trigger}(v)^2
=3890796-158(157)^2+4(157)(6)=22,
\]

so its wedge count is W2=(22-2m)/2=5. The existing black-NN rank-one torus
bipartite-trigger lemma (`notes/p334-trigger-bipartite-theorem.md`) gives zero
trigger triangles. Therefore

\[
f^{\rm pair}_3=\binom{158}{3}-6(156)+5=644025,
\qquad g_3=f^{\rm pair}_3-f_3=5\ \text{or}\ 19.
\]

This is exact integer algebra from the already stored source moment and the
new full-clock coefficient, not another triple enumeration. The unconditional
third-step survival gap is exactly `7/424023`. The two safe polynomials start

```
A: 1, 158, 12397, 644020, 24919518, 766033288, ...
B: 1, 158, 12397, 644006, 24917380, 765871160, ...
```

It is now possible to say precisely what extra mechanism first breaks this
scalar closure: irreducible three-site triggering, not different singleton
gates, pair abundance, or pair overlap as measured by the degree second
moment. This does not attribute the entire later clock difference to g3;
higher-order structure also differs. No CFT field identification or path-memory
interpretation follows merely from this finite-state witness.

## Direct versus collective competition from the same polynomial

Assign independent uniform insertion labels u to the remaining sites and use
s=-log(1-u). With h=H2 original direct-trigger sites, their first clock is
Exp(h), independent of the collective event on the other sites. Thus

\[
E[s_{\rm birth}]=\sum_{k=0}^{172}\frac{f_k}{173\binom{172}{k}},\qquad
P(\text{original direct site wins})=hE[s_{\rm birth}].
\]

Every original direct site has the same winning share E[s_birth]. The summary
saves exact direct and collective shares for all 147 without any sitewise DP.
Even at H2=2, the direct winning fraction ranges from 0.18625 to 0.32901: the
same number of immediate gates competes differently with the collective clock.

## Deliverables and boundary

Full clocks and factor coefficients were pushed first in `9cca7bc6` so the
parallel thermal-mixture analysis could consume them immediately.
`scripts/p334_all147_prefix_clocks.py` reproduces selection/replay and the 135
new evaluations while reusing old twelve artifacts; the separate summary
script uses saved polynomials only. No new MC, remote job, sitewise pivotal
sweep, or full-repository test campaign was run.

These 147 rows share one old production dependency block. Their complete
conditional clocks can now be reused for exact conditional thermal readouts
and within/between-prefix variance decomposition. Such reuse is a new readout
of this source, not 147 newly independent simulations and not an unconditional
random-prefix ensemble claim.
