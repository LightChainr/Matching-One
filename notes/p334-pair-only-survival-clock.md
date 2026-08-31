# Two real trigger graphs: more branching, uniformly later pair absorption

**Result:** the high-branching N425 checkpoint B does not switch between early
and late pair triggering. Its complete pair-only survival curve dominates A,
and its discrete pair-only hazard is smaller at every common nontrivial step.
This is exact for the two saved graphs, not a population result or a statement
about their full ambient-rank birth times.

The source is the same pair of real second-orientation N425 checkpoints as
`6147e22` and `c1fbcc4`: A counter `43042514269`, B counter `43042505280`, seed
`20260831430425`, period matrix `[[425,268],[0,1]]`, k0=252, age=10,
line `(12,-19)`, H2=0, a=d=173 and 108 minimal trigger edges each.

## 1. The independent-set polynomial is the exact pair-only clock

Keep this checkpoint's minimal-trigger-pair graph T fixed. With `i_k` its
number of independent k-vertex sets, define

```text
I_T(z)=sum_k i_k z^k,
S_pair(k)=P(T_pair>k)=i_k/C(d,k),
h_pair(k)=1-S_pair(k)/S_pair(k-1)
         =1-k*i_k/[(d-k+1)*i_(k-1)].
```

The hazard is undefined after survival has become zero. A uniform future
permutation has a uniform k-subset as its first k sites, so these equations
need no independence approximation between trigger edges.

For bipartition L/R and f isolated vertices,

```text
I_T(z)=(1+z)^f * sum_(U subset L) z^|U| (1+z)^(|R|-|N(U)|).    (1)
```

Enumerating the smaller side takes only **4096 subsets for A and 32 for B**.
The output retains all 174 exact integer counts, rational survival/hazard
values, the (|U|,free-opposite-sites) histogram and the actual site labels.
No new MC, graph reconstruction, H2/b2 or triple census is run.

## 2. Both real graphs have a short bottleneck decomposition

The committed adjacency gives an even simpler form:

- A has 147 isolates. Site **159** is adjacent to all 12 opposite-side sites;
  deleting it leaves `K_(11,8) disjoint-union K_(2,4)`.
- B has 139 isolates. Site **121** is adjacent to all 29 opposite-side sites;
  deleting it leaves `K_(3,25) disjoint-union K_(1,4)`.

These are exact statements about the two graphs, not proposed universal
classes. In particular they are not Ferrers graphs. Write `x=1+z`. Splitting
on whether the universal site is included gives the complete polynomials:

```text
I_A(z)=x^147 * [z*x^13 + (x^11+x^8-1)*(x^2+x^4-1)],
I_B(z)=x^139 * [z*x^4  + (x^3+x^25-1)*(x+x^4-1)].             (2)
```

The first coefficients are:

| k | i_k(A) | i_k(B) |
|---:|---:|---:|
| 0 | 1 | 1 |
| 1 | 173 | 173 |
| 2 | 14770 | 14770 |
| 3 | 830504 | 831044 |
| 4 | 34625825 | 34709231 |
| 5 | 1142459004 | 1148862594 |
| 6 | 31088585941 | 31414423490 |

## 3. No crossing, even in the stronger hazard order

The difference has the exact positive factorization

```text
I_B(z)-I_A(z)=z^3*(1+z)^139*Q(z),
```

where Q has degree26 and, in increasing-power order, coefficients

```text
[540,8346,64356,331519,1275568,3882024,9669438,20154212,
 35688718,54262286,71364333,81592028,81342874,70808877,
 53810642,35637729,20500296,10190102,4345000,1573406,
 477319,119055,23776,3655,406,29,1].
```

Every coefficient is strictly positive. Hence B survival is strictly larger
for **all k=3,...,168**, equal at k=0,1,2 and zero for both from k=169 onward.
The largest pair-safe subset has size162 for A and168 for B.

For hazard, the exact sign is the sign of the integer

```text
C_k=i_k(A)*i_(k-1)(B)-i_k(B)*i_(k-1)(A).
```

All values are saved. They are zero at k=1,2 and strictly negative at every
**k=3,...,163**, the entire remaining range where both hazards are defined.
After k=163, A has no pair survivors and no conditional hazard to compare.
B finally has hazard1 at k=169. No floating-point crossing search is used.

| pair-only statistic | A | B |
|---|---:|---:|
| mean first-trigger step | 3755998829/173838600 = 21.60624 | 28111063/885360 = 31.75100 |
| 10% / 50% / 90% quantile | 7 / 19 / 40 | 7 / 26 / 65 |
| S(10) | 0.7839011 | 0.8229963 |
| S(20) | 0.4507272 | 0.5945725 |
| S(40) | 0.0961174 | 0.2923020 |
| S(65) | 0.0079358 | 0.0997152 |
| h(20) | 0.0618899 | 0.0329185 |
| h(40) | 0.0837309 | 0.0373596 |

The exact mean delay is `765380741759/75445952400 = 10.1447555` steps.
Thus the larger degree-overlap/branching response of B corresponds to a
uniformly delayed **pair** absorption clock, not merely a third-step effect
or an early/late crossing. Equation (2) makes the structural reason visible:
many candidate pairs share a small set of gate-side sites, leaving a broad
collection of subsets that avoid completing any pair.

## 4. Genuine triples remain a separate clock-shortening mechanism

The true rank-two time can precede the first original trigger edge, because
three or more individually pair-safe sites can already complete rank two.
Monotonicity gives the pathwise inequality

```text
T_true <= T_pair,   S_true(k) <= S_pair(k).
```

This makes the pair-only mean an upper bound on each true mean. It does **not**
give a pointwise hazard bound or transfer the A/B survival ordering to the
complete birth process.

The existing exact triple census gives

```text
i_3(A)-safe_true_3(A)=830504-829921=583,
i_3(B)-safe_true_3(B)=831044-830535=509.
```

At k=3 the true B-A survival gap is therefore614/C(173,3), compared with the
pair-only540/C(173,3). This is the already known triple correction, retained
as a distinct component. At k>=4 the complete minimal-trigger hypergraph is
not enumerated here. Whether the **full** birth curves cross remains open;
the present exact result settles the entire pair-only question.

## Scientific boundary and reproduction

This is zero-new-sample exact analysis of **two specified real configurations**,
reusing the original cooperative N425 dependency block. It does not sample the
population, identify M loading or a continuum field, or prove that arbitrary
high-W2 graphs have this ordering. No old tests or full regression are rerun.
The one calculation retains exact sign certificates and the existing k<=3
integer identities as inline bookkeeping.

```bash
/Users/lc/python-envs/research-py311/bin/python scripts/p334_pair_only_survival.py
```

Output: `results/p334-pair-only-clock/pair_only_survival.json`.
