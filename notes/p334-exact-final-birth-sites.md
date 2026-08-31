# P215 / P334: exact final-birth site allocation on the full physical clock

## Outcome: concentrated ordinary births, but a more distributed B tail

The final insertion site is now resolved exactly for both real N425 prefixes.
This uses the full physical two-terminal event of `6358ba49`, not a trigger
truncation or a sampled continuation.

| Conditional finite-prefix result | A | B |
| --- | ---: | ---: |
| Sites with strictly positive final-birth probability | 83 | 127 |
| Exactly irrelevant vacant sites, out of 173 | 90 | 46 |
| Inverse-Simpson effective number of birth sites | 36.4300 | 21.9549 |
| Probability assigned to the five highest-probability sites | 22.2550% | 43.0173% |
| Share of T>40 births assigned to the five largest late contributors | 22.9465% | 32.3775% |
| Interior-site share of all births | 12.4462% | 23.6534% |
| Interior-site share among T>40 births | 18.4034% | 40.8146% |

**B has a larger exact influential support but a more concentrated ordinary
birth allocation. Its long tail does not concentrate further on those main
ports: it shifts toward interior completion.**

In B the leading sites are `{121,8,166,277,390}`, all on port 0. Together
they account for 43.0173% of all births but only 32.3775% of births after step
40. The other sites finish 67.6225% of that tail. Site 121 is individually
largest, but its mean birth time conditional on being the final site is
18.6282, below B's overall mean 20.7788. It is not a uniquely late bottleneck.

| B final site | Whole-clock probability | E[T given final site] |
| --- | ---: | ---: |
| 121 | 0.09752434 | 18.62817 |
| 8 | 0.09392448 | 18.97847 |
| 166 | 0.09392448 | 18.97847 |
| 277 | 0.09392448 | 18.97847 |
| 390 | 0.05087476 | 22.87999 |
| 211 | 0.03382598 | 23.93871 |
| 257 | 0.03382598 | 23.93871 |
| 144, interior | 0.02437087 | 26.26574 |
| 413, interior | 0.01770113 | 27.62860 |

For A, the eight sites `{6,27,140,162,251,274,296,409}` are tied at
0.04450992824 each. Its reported top-five number uses the deterministic label
ordering among that tie; it does not identify a special five-site geometry.

## Exact joint law, not just integrated influence

Fix one of the original rank-one prefixes and a vacant site v. For the other
172 sites, let

\[
 a_{v,0}(k)=\#\{U: |U|=k,\ \mathrm{rank}(B+U)=1\},
\]

\[
 a_{v,1}(k)=\#\{U: |U|=k,\ \mathrm{rank}(B+U+v)=1\}.
\]

The nonnegative integer difference `D_v(k)=a_v,0(k)-a_v,1(k)` counts sets
for which adding v creates the first rank-two birth. By monotonicity, every
ordering of a safe U has been safe before v is added. Under a fresh uniform
permutation of all d=173 vacant sites,

\[
\Pr(T=k+1,\ V_{\rm final}=v)=
\frac{D_v(k)}{d\binom{d-1}{k}}.
\]

Thus the exact probability that v is the final birth insertion is

\[
\pi_v=\sum_{k=0}^{172}\frac{D_v(k)}{173\binom{172}{k}}
     =\int_0^1 \partial_{p_v}\Pr_p(\text{rank-two connection})\,dp.
\]

The derivative here keeps the other site probabilities equal to p. Integrating
`p^k(1-p)^(172-k)` gives exactly the displayed factorial denominator.

All forced-off/on polynomials use the already constructed physical two-port
network and the same connectivity-partition tree DP. A forced site has no z
weight; the remaining 172 Bernoulli/site choices are counted exactly. Outside
the terminal-path core, D is identically zero without a DP call.

The resulting **full joint law** is saved as the exact integer arrays
`pivotal_count_by_prior_size`, plus rational pi, conditional mean times, and
the T>40/T>65 joint masses for every original site label. For both real states,

\[
 \sum_v\pi_v=1,\qquad
 \sum_v\frac{D_v(k)}{173\binom{172}{k}}=S(k)-S(k+1)
\]

hold exactly at every step. The local marked-birth observer therefore sums
back to the complete global birth clock, not merely to a low-order proxy.

Since every D coefficient is nonnegative and every integration weight is
positive, `pi_v=0` is equivalent to `D_v(k)=0` for every k. The 90/46 zero-pi
sites are consequently **exactly irrelevant to the entire Boolean birth
event**, not just rarely observed under one sampled continuation. In the
previous terminal cores, another 39/19 sites beyond the off-core 51/27 turn
out to be irrelevant because of the fixed occupied background.

## Port-resolved allocation

Ports are the two exact gain-address classes of the previous map, ordered by
their transverse address; they are not a new angular/CFT label. The share
by physical network role is:

| Role | A all births | A T>40 | B all births | B T>40 |
| --- | ---: | ---: | ---: | ---: |
| Port 0 | 49.3538% | 62.0347% | 57.2554% | 53.9591% |
| Port 1 | 38.2000% | 19.5619% | 19.0912% | 5.2263% |
| Interior | 12.4462% | 18.4034% | 23.6534% | 40.8146% |

This pins down a real time-dependent mechanism: the early-to-late completion
mix changes substantially, especially in B. It connects the P215 insertion
idea to the complete P334 clock while retaining the exact site labels and
step distribution needed for a subsequently defined marked observer. No H4
or field-identity interpretation is assigned here.

## Scope, cost and next production use

Same two N425 second-orientation real checkpoints, k0=252, seed
20260831430425, counters43042514269/43042505280, age10, ell=(12,-19).
No new Monte Carlo, no independent production block, no remote host. One
sequential local process took 2.64 s for A and 57.63 s for B.

The proposed next production use is recorded, **not executed here**: after
sampling a first-birth prefix, replace random future completion by its exact
conditional full F2(k|prefix), jointly returning marked-birth probabilities.
This would Rao–Blackwellize the K2 and derived A/E thermal-curve observations.
The next gate is a small empirical prefix-cost map; treewidth 4/6 on these
two examples does not establish tractability for all prefixes or N900.

The reported T>40 comparisons are a finite-instance descriptive decomposition,
not a population hypothesis test. The integral allocation is exact for these
two conditional events, not a universal spatial law.

```sh
/Users/lc/python-envs/research-py311/bin/python scripts/p334_exact_marked_birth.py
```
