# Bipartite capacity forces the branching order: sharp overlap envelopes

**New result:** the two saved real N425 checkpoints have disjoint admissible
two-star intervals even before their detailed trigger edges are specified.
This replaces a random-graph expectation comparison by deterministic bounds.
It reuses `119cb5f` bipartiteness and the existing 22 saved graphs; there are
no new configurations, random paths, H2/b2 measurements or triple enumerations.

## 1. The matched real pair cannot reverse its branching order

The N425 second-orientation witnesses `43042514269` (A) and `43042505280` (B)
share `d=a=173`, `m=108`, age 10, line `(12,-19)`, H2=0 and safe-pair count
14770. Their nonisolated trigger-component capacities differ.

| fixed capacity | sharp minimum W2 | observed W2 | sharp maximum W2 | position in range |
|---|---:|---:|---:|---:|
| A: 14 x 12, 108 edges | 796 | 926 | 1046 | 52.000% |
| B: 5 x 29, 108 edges | 1263 | 1466 | 1578 | 64.444% |

**B's minimum exceeds A's maximum by 217.** Every edge arrangement in these
declared capacity classes therefore has

```text
W2(B)-W2(A) >= 217,
branch_survival(B)-branch_survival(A)
  = Delta_coop(B)-Delta_coop(A)
  >= 217/2559016 = 0.0000847982193.
```

The actual W2 gap is 540 and the actual branching probability gap is
135/639754. Thus **217/540=40.1852%** of the observed contrast cannot be removed
by any such edge rewiring. This is not the earlier 84% decomposition of a
uniform-bipartite expectation; it is a worst-case ordering statement.

Scope matters: the bounds are **sharp over all simple bipartite graphs with
the specified vertex slots, side capacities and edge counts**, permitting
isolates or disconnection inside each slot block. We do not assert that their
extremizers are realizable lattice checkpoints. They are therefore valid,
possibly conservative bounds for the narrower physical class. For multiple
observed components we forbid new cross-component edges and sum each block's
extrema, but do not impose continued connectedness inside each block.

## 2. Exact lower bound: balance both sets of degrees

For `s` slots and `m` incidences define

```text
q=floor(m/s), r=m mod s,
B(s,m)=s*C(q,2)+r*q.
```

For a simple bipartite graph with capacities L,R and m edges,

```text
W2_min(L,R,m) = B(L,m)+B(R,m).                  (1)
```

Convexity of `C(degree,2)` gives the lower bound separately on each side.
Both bounds can be attained simultaneously: give the L rows lengths
`floor(m/L)` or `ceil(m/L)`; place their consecutive edge slots in columns
`0,1,...,R-1,0,1,...` without restarting between rows. No row has more than R
slots, so the graph is simple; column counts also differ by at most one.
Equation (1) is therefore exact, not a Jensen-only relaxation.

## 3. Exact upper bound: Ferrers is an extremizer, not the physical model

Order row degrees from largest to smallest. If two row neighbourhoods are
not nested, move an edge from the lower-degree row to a vacant position in
the higher-degree row in the same column. Column degrees stay fixed and W2
increases by `d_high-d_low+1>0`. Thus some maximum has nested rows, represented
by an integer partition

```text
R >= lambda_1 >= ... >= lambda_L >= 0,
sum lambda_i=m.
```

Its exact objective is

```text
W2_max = max_lambda sum_i [C(lambda_i,2)+(i-1)*lambda_i].       (2)
```

The second term counts column two-stars. The short dynamic program in the
scorer solves (2) exactly with state `(row,previous_degree,remaining_edges)`;
it enumerates partitions, not graphs or lattice configurations. For A an
extremizing right-degree partition is `(14,14,14,14,14,14,14,10,0,0,0,0)`;
for B an extremizing left partition is `(29,29,29,21,0)`.

The fact that an **extremizer** is Ferrers does not revive Ferrers as a model
for the data: A and B are the previously certified non-Ferrers graphs.

## 4. What collectively amplifies the finite-time fluctuations?

Let `t_v` be a safe vertex's trigger degree, `c_v=a-1-t_v` its safe-successor
degree, `d` the number of all vacant sites, and `m` the trigger-edge count.
The existing same-checkpoint branching identity becomes

```text
D=d*(d-1)^2,
D*Delta_coop = sum_v t_v^2 - 4*m^2/a
             = 2m+2W2-4m^2/a.                              (3)
```

For fixed capacity blocks `(Lj,Rj,mj)`, equation (3) has the exact decomposition

```text
D*Delta_coop =
  [sum_j mj^2*(1/Lj+1/Rj) - 4*m^2/a]
  + sum_j [sum_(v in Lj)(t_v-mj/Lj)^2
           +sum_(v in Rj)(t_v-mj/Rj)^2].                    (4)
```

The first term is forced by the support/side capacities; the second is
within-side degree organization. Applying the integer-sharp lower bound adds
the unavoidable rounding variance to the first term.

| checkpoint | minimum Delta_coop | observed | maximum | minimum / observed |
|---|---:|---:|---:|---:|
| A | 0.0003005671 | 0.0003513679 | 0.0003982609 | **85.5420%** |
| B | 0.0004830591 | 0.0005623865 | 0.0006061533 | **85.8945%** |

Thus these high-variance states do not need nearly maximal within-side
clustering. Most of their cooperative variance is already unavoidable when
many safe sites lie outside a relatively narrow, asymmetric two-sided trigger
support. The actual two-star counts are only 52% and 64.44% through their
available ranges. This is an algebraic fixed-state attribution, not causality
or a claim about population-wide percentages.

For the **same previously selected 22 graphs**, the sharp lower floor accounts
for 82.1889%--100% of each observed Delta_coop. Eight have zero-width W2
envelopes: their capacity/edge data already determine W2 exactly. This is a
bounded descriptive result, not an estimate from all production checkpoints.

## 5. More branch dispersion need not mean faster absorption

The `119cb5f` theorem removes trigger triangles. With `c3` the genuine minimal
three-site trigger count, the existing exact survival identity is

```text
s2=[C(a,2)-m]/C(d,2),
s3=[C(a,3)-m*(a-2)+W2-c3]/C(d,3),
h3=1-s3/s2,
d h3 / d W2 = -3/[(d-2)*(C(a,2)-m)]    at fixed c3.          (5)
```

W2 raises clone dispersion but **lowers** the third-step conditional exit
hazard at fixed c3: overlapping trigger edges cover fewer distinct triples.
The sharp W2 bounds give sharp pair-clique triple-count envelopes; subtracting
a specified c3 gives an algebraic slice, not a statement that c3 is freely
rewirable or that every endpoint is physically attainable.

Using only the already measured `c3(A)=583`, `c3(B)=509`, the real B-A
third-step hazard difference is

```text
h3(B)-h3(A) = -307/420945 = -0.00072931143.
```

If these two c3 values are additionally held fixed, even the worst allowed W2
rewiring leaves `h3(B)-h3(A)<=-97/280630`. **Capacity alone forces branching
order; it does not force three-step hazard order when c3 is allowed to change.**
The new mechanism is therefore concentrated two-sided trigger support plus
degree-overlap organization, with an explicitly separate higher-order layer,
not a single scalar "more cooperation means earlier birth" narrative.

## Reproduction and limits

`results/p334-sharp-overlap-envelope/sharp_envelopes.json` retains every selected
graph, exact fractions, component capacity, extremizing partition and response
envelope. Source is the committed `trigger_graph_structure_bounded.json` and
`safe_triple_census.json` under the existing cooperative production directory.
All observations retain the original N325/N425 dependency groups.

The three focused new checks compare (1)--(2) against all tiny 1x4/2x3/3x3/3x4
bipartite graphs, inspect extremizing partitions, and check the separated real
pair envelopes. No old graph-property tests or full repository suite are rerun.
No continuum field, M loading, scaling exponent or fresh independent sample is
claimed.

```bash
/Users/lc/python-envs/research-py311/bin/python scripts/p334_sharp_overlap_envelope.py
/Users/lc/python-envs/research-py311/bin/python -m pytest -q tests/test_p334_sharp_overlap_envelope.py
```
