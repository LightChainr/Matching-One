# P334: twelve frozen real-prefix clocks reveal genuine time-order crossings

## Result

**All 12 new real prefixes are solved exactly**, with no replacement, new
Monte Carlo, server, or sitewise-pi sweep. The original two worked examples
are excluded. Every selected prefix has the same N425 second geometry,
k0=252, age10 and ambient line (12,-19).

The complete mean waits range from **5.77269 to 24.89885 steps**. The
probability of waiting beyond step40 ranges from **0.00006809 to 0.11284544**.
Most importantly, **11 of the 66 fixed-prefix pairs have crossing complete
survival curves**. The previous two examples' all-time ordering is therefore
not a general feature even within this tightly matched geometry/age/line.

| Counter minus 43042500000 | H2 | Parallel two-port factors | E[T] | P(T>40) | Treewidth upper bound |
| --- | ---: | ---: | ---: | ---: | ---: |
| 48 | 9 | 5 | 13.322716 | 0.01870887 | 5 |
| 83 | 6 | 3 | 17.977368 | 0.06182175 | 6 |
| 106 | 7 | 3 | 11.736249 | 0.00703240 | 5 |
| 622 | 15 | 5 | 9.928906 | 0.00538598 | 6 |
| 835 | 4 | 2 | 19.153133 | 0.06310962 | 6 |
| 881 | 2 | 1 | 21.085662 | 0.06232955 | 4 |
| 904 | 28 | 14 | 5.772688 | 0.00006809 | 4 |
| 1006 | 5 | 2 | 16.804840 | 0.03524933 | 4 |
| 1010 | 3 | 1 | 16.779479 | 0.02582011 | 5 |
| 1013 | 21 | 9 | 7.435020 | 0.00079860 | 4 |
| 1043 | 0 | 1 | 24.898851 | 0.11284544 | 7 |
| 1410 | 2 | 1 | 19.700589 | 0.05323451 | 6 |

The reported completed-result wall times sum to approximately **0.56 s**,
with a maximum of **0.28 s** for a prefix; the initial structure scan adds
only a small overhead. No prefix reaches its frozen 20 s / 200,000-state
limit. Core variable counts range from 87 to151. These are structural core
counts, **not** the exact pivotal support/effective-site counts: pi_v was
deliberately not computed for this batch.

## Frozen selection and provenance

Selection was committed as `b9cbe13e` before any new-prefix mapping or clock
evaluation. The source is the existing
`results/local-20260831/P334-cooperative-closure/raw/N425.geometry_pilot.csv`,
with its original metadata. Required fields were N425, second orientation,
k0=252, age_steps10, ell_u12, ell_v-19, and rank one at the checkpoint
(`k1 <= k0 < k2`). Exclude counters43042514269/43042505280, then sort numeric
counters and take the first12. **H2, b2, runtime, width and outcomes played
no role in selection.** There were147 eligible rows after exclusions.

`selected_prefixes.json` preserves the original rows and source/manifest
SHA256 values. Each occupied prefix is decoded from the original seed
20260831430425 and counter, rather than newly sampled. This is a new
deterministic analysis of an old production stream, not a new independent
Monte Carlo replication block or a population estimate from twelve cases.

## New structural increment: several exact parallel channels

Four prefixes have one double-address component, as in the original examples.
Eight instead have 2 to14 **site-disjoint double-address components** around
the same contracted essential occupied cluster. In all twelve saved maps:

- the entire graph with that cluster removed is q-balanced;
- each component has at most two port addresses;
- every two-address difference is exactly425;
- no three-address or gain-consistency obstruction occurs.

The first single-factor scorer reported these eight as not-single-factor;
that initial pass is retained. They are not failures of the physical map.
The exact extension is simply

\[
 F(z)=(1+z)^{n_{\rm free}}\prod_{j=1}^{m} F_j(z),
 \qquad S(k)=\frac{[z^k]F(z)}{\binom{173}{k}},
\]

where F_j is the safe-subset polynomial of one two-port component's terminal
core. A transverse birth occurs if **any** component connects its ports;
safety requires all components to remain disconnected. Their random site
sets are disjoint, so their generating polynomials multiply exactly. This is
not the incorrect product of fixed-k component survival probabilities.

The original four completed clocks and all exported maps were reused; only
the other eight received the product calculation. Selection and resource
limits were unchanged. Every full polynomial retains the original exact
singleton/pair coefficients without rerunning those censuses.

## Concrete time-order reversal: stronger immediate gate, slower later clock

For a readable post-hoc example, define

- A = counter43042501006: H2=5, b2=13994;
- B = counter43042500083: H2=6, b2=13845.

B has **more immediate triggers** and fewer safe pairs. Nevertheless,

\[
S_B(k)<S_A(k)\quad(k=1,\ldots,10),\qquad
S_B(k)>S_A(k)\quad(k=11,\ldots,153).
\]

Its mean wait is longer by exactly

\[
\frac{18201264770939272531633063727}{15523108758334473814211393700}
=1.1725270404465133\ldots,
\]

and its step40 tail is 1.75384 times A's. This is a material crossing around
survival0.65, not merely an asymptotic zero-probability tail crossing. The
instantaneous hazard order has its own changes: h_B>h_A at steps1..5,
h_B<h_A at6..83, then h_B>h_A again at84..153. Cumulative survival dominance
must not be confused with permanent hazard dominance.

This example was chosen descriptively **after** computing all frozen clocks.
All eleven crossing pairs, including additional very-late crossings, are
reported with exact integer sign certificates in `scientific_summary.json`;
11/66 is not an estimated ensemble crossing probability.

## Exact direct-versus-collective decomposition explains the crossing

Let a_i=173-H2_i and f_i(k) be the true safe-set coefficient. Then

\[
 S_i(k)=
 \underbrace{\frac{\binom{a_i}{k}}{\binom{173}{k}}}_{G_i(k):\ \text{avoid original single-site triggers}}
 \underbrace{\frac{f_i(k)}{\binom{a_i}{k}}}_{C_i(k):\ \text{remaining collective survival}}.
\]

For this pair a_A=168, a_B=167, so the direct-avoidance ratio is especially
simple: `G_B/G_A=(168-k)/168`. It always works against B. The actual
collective survival advantage overcomes it between steps10 and11:

| k | Direct ratio G_B/G_A | Collective ratio C_B/C_A | Total S_B/S_A |
| --- | ---: | ---: | ---: |
| 1 | 0.99404762 | 1.00000000 | 0.99404762 |
| 5 | 0.97023810 | 1.01297088 | 0.98282294 |
| 10 | 0.94047619 | 1.06021678 | 0.99710863 |
| 11 | 0.93452381 | 1.07406181 | 1.00373633 |
| 20 | 0.88095238 | 1.27172951 | 1.12033314 |
| 40 | 0.76190476 | 2.30191715 | 1.75384164 |

One early structural clue is already in the old exact counts: after excluding
the single-site triggers, A has `C(168,2)-13994=34` minimal trigger pairs,
whereas B has only `C(167,2)-13845=16`. The additional immediate gate in B is
opposed by weaker collective completion among the remaining safe sites.
The full coefficient calculation establishes how that opposition develops
over time; neither H2 alone nor the count of parallel channels orders the
whole clock.

## Scope and next use

These are complete physical conditional clocks on twelve fixed real
rank-one prefixes, with full gain/gauge/port certificates saved under `maps/`.
No population law, CFT field identity, or universal small-width claim follows.
The core sizes and time bound support trying exact conditional continuations
on a small production prefix cost map, but do not guarantee all N900 states
are tractable. No such production launch or new pi_v analysis was made here.

```sh
/Users/lc/python-envs/research-py311/bin/python scripts/p334_twelve_prefix_clocks.py
/Users/lc/python-envs/research-py311/bin/python scripts/p334_twelve_prefix_summary.py
```
