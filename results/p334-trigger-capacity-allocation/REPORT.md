# P334: side imbalance, not bipartiteness alone, explains the selected overlap contrast

**The saved N425 pair separates two effects that the previous 84% capacity figure combined.** Expanding the active support from 26 to 34 vertices lowers the fixed-edge simple-graph expected two-star count; the much more unequal 5×29 bipartition raises it strongly enough to reverse that change. Both real graphs retain overlap beyond their capacity benchmarks.

This is exact zero-new-sample arithmetic on the 22 already saved trigger graphs. It does not repeat graph replay, bipartite/Ferrers screening, the full-production overlap scorer, or the minimal-triple census.

## Three additive structural comparisons

Edges are minimal rank-two-triggering pairs on individually safe sites; W2 counts shared-endpoint edge pairs. Keep a safe sites and m trigger edges. Compare: (i) G(a,m); (ii) uniform simple edges independently within each observed nonisolated component's vertex block with its observed edge count; (iii) uniform cross edges within that same block's observed L/R sides; (iv) the observed graph.

The last two benchmarks may disconnect a block or introduce new isolates. They do not preserve connectedness or the complete component decomposition. The first increment is support/block/edge-allocation localization, not purely support size or a causal physical intervention.

For one block with s=L+R, exact expectations are `E_simple W2=2m(m-1)/(s+1)` and `E_bip W2=m(m-1)(s-2)/[2(LR-1)]`. Hence

```text
E_bip W2 - E_simple W2
  = m(m-1) [(L-R)^2-s+2] / [2(LR-1)(s+1)].
```

For m≥2, bipartition constraints increase expected overlap only when `(L-R)^2 > s-2`. Near-balanced bipartite capacities can reduce overlap. The exceptional one-slot case has both expectations zero.

## Matched real N425 witnesses

A/B remain the previously selected counters 43042514269 / 43042505280, with identical a=d=173, m=108, k0=252, age=10 and ell=(12,-19). Their single nonisolated blocks are 14×12 and 5×29.

| W2 quantity | A | B | B−A |
|---|---:|---:|---:|
| global_null_W2 | 132.827586 | 132.827586 | +0.000000 |
| block_simple_null_W2 | 856.000000 | 660.342857 | -195.657143 |
| block_bipartite_null_W2 | 830.371257 | 1284.000000 | +453.628743 |
| observed_W2 | 926.000000 | 1466.000000 | +540.000000 |

| Contribution to observed B−A = 540 | exact | decimal |
|---|---:|---:|
| localization | -6848/35 | -195.657143 |
| bipartition_constraint | 3795076/5845 | +649.285885 |
| residual_organization | 14424/167 | +86.371257 |

Thus the already reported +453.628743 capacity contribution is **−195.657143 localization +649.285885 side constraint**. Residual organization contributes +86.371257. This refines, rather than repeats, the earlier 84.0053% arithmetic: two-sided imbalance overcompensates an opposing support effect. A's side constraint lowers its expected W2 by 25.628743; B's raises it by 623.657143.

Multiplying each contribution by `2/[173·172²]` gives an exact decomposition of the saved double-clone probability difference `135/639754`; the JSON stores these rational contributions. There is no new independent branching evidence.

## Existing first-five extension: fixed-set descriptive check

The two specially chosen witnesses are excluded from this table. Each row uses the already fixed five lowest eligible counters in that environment. Entries are sums of W2 increments across those five saved graphs, not population means, prevalence estimates or a new sampling test.

| fixed environment | support/block localization | side constraint | residual organization | positive side increments / 5 | positive residuals / 5 |
|---|---:|---:|---:|---:|---:|
| N325_first | +193.0237 | +100.0561 | +7.6407 | 2/5 | 1/5 |
| N325_second | +117.4221 | -7.0691 | -2.8222 | 1/5 | 1/5 |
| N425_first | +1288.7169 | +80.8429 | +101.2903 | 4/5 | 5/5 |
| N425_second | +470.5864 | -2.2601 | +13.9348 | 2/5 | 2/5 |

Only 9/20 of these fixed extension graphs have strictly positive residual W2 beyond the bipartite-capacity expectation. The N325-second residual sum is negative, whereas the other three are positive. This does not estimate population prevalence; it shows why the full-production excess above G(a,m) cannot be interpreted as excess beyond every geometry-aware benchmark. The selected A/B positive residuals are not representative evidence for that stronger statement.

## Scientific consequence and next output

The microscopic target is now sharper than 'a bipartite graph' or 'more hidden memory': explain the physical two-sided capacity and the residual nonexchangeable degree organization. Ferrers nesting is already false; none of this proves universal bipartiteness. The separately measured genuine triple-trigger counts 583/509 remain needed for three-step survival and are not supplied by this pair-graph benchmark.

**Next actual output:** use the saved graph/site labels and period matrix to identify the two sides with explicit topological cut or boundary-landing classes on the existing matched configurations, and predict their capacities before consulting degrees. That would turn the presently observed L/R partition into a physical explanation. Another graph census, generic third-clone moment, or repeated support-baseline score is not required to reach this question.

Input: `1b5a9dea07e1c62f69798fddbf4899ff986c0b72` / `results/local-20260831/P334-cooperative-closure/trigger_graph_structure_bounded.json`, SHA256 `2765a8181051c519be89bf04219239b4d49dd2bc4df788216ae62047d4664300`. Background full-production overlap source: `85af4981628541882fd038ef07234b9a2f2f0266`. All rows remain in the original P334 N325/N425 dependency groups.

All additions and benchmark identities use exact rational arithmetic. No graph simulations, topology runner, parent scorer, or test suite were run. One script produces both outputs:

```bash
/Users/lc/python-envs/research-py311/bin/python scripts/p334_trigger_capacity_allocation.py
```
