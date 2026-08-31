# N25 canonical Kreg joint-source producer and completed run

`scripts/p337_regular_pair_joint_exact.cpp` provides the missing joint-source moments for the fixed contract `4ce4dfe894c9fe96f268c61cf21eb6585dba5418:analysis/p337_regular_pair_joint_contract.json`. The pre-data scorer is `5da4749245450048625a2da43e8f73da1ee9275c`. After acceptance of both theory gates and the root's explicit GO, the producer was compiled once and both complete traversals ran once. Raw outputs and the complete execution receipt are in `results/p337-regular-pair-joint/`. No score was calculated here.

## Fixed source and integer units

For the full occupation population, `b16(A)=sum_{y != 0} g16(0,y;A)` if the origin is vacant; occupied origin and occupied y give zero source. For translation-invariant `K,q,E`, the homogeneous second-source moment uses `s2=b16/(16*N)=b16/400`. The ordered-pair/translation reduction already accounts for the second derivative's factor of two: **do not multiply by two or by N again**. This is the signed joint tensor kernel, not a product or covariance of the one-site activation mark.

The only decomposition is fixed before production: y is an original NN neighbor of the origin (`adj`) or another distinct site (`far`). All configurations remain in the count, including occupied-origin configurations. The CSV schema is

```text
k,count,sum_q,sum_e,sum_b16,sum_qb16,sum_eb16,sum_b16_adj,sum_qb16_adj,sum_eb16_adj,sum_b16_far,sum_qb16_far,sum_eb16_far
```

Every source/cross-source sum uses `int64_t`. Total is formed as adjacent plus nonadjacent at each configuration, before aggregation. Divide all nine source columns by 400 exactly once. There is no scan over radii, angles, support classes or other sources.

## Physical edge identification and preserved backend

The N25 `(5,0)` and `(4,3)` quotient construction, black/white rollback DSUs, face counting, binary recursion and original `q=CB-CW-(K-B+F)`, `E=q^2` are copied without changing their semantics from `scripts/p337_regular_pair_activation_exact.cpp` at base `a237968f1d7a82d26b46e83c58179dbba7f1a908` (blob `df4ca4db1c894ea649a660afe6b14ad0923b2ae3`). No inserted virtual wiring enters q/E.

Ports are strictly `x:N,E,S,W,y:N,E,S,W`. At a vacant mark, occupied neighbors receive their original occupied-NN component root. A vacant neighbor leaves an isolated **physical edge-node**, labeled `N+edgeID`. Edge IDs are global on the quotient: `2*v` for the north edge out of v and `2*v+1` for its east edge. South and west ports refer respectively to the north/east edge out of their neighbor. Thus two adjacent vacant marks assign their shared edge exactly the same singleton ID; different edges incident to the same vacant vertex remain distinct. Black-component root IDs `0..N-1` cannot collide with edge-node IDs `N..3*N-1`.

The eight IDs are canonicalized in first-occurrence restricted-growth order and packed as `sum(label[i]<<(3*i))`. The signed exact Bell8 kernel is reused unchanged; sparse omissions are zero. A dense int16 lookup supports the fixed 24-bit key domain while keeping every accumulated source sum int64. Occupied-component roots are cached once per configuration; only origin-to-y pairs are inspected, never all O(N^2) pairs.

## Provenance and executed invocation

- New producer SHA256: `1db71a287d106f0c3eb60a3d02987737059c37842bcdf4fcba1003e95b227418`.
- Kernel: `analysis/regular_pair_spatial_kernel.tsv`, SHA256 `36ae069d370b1d7a4398861c928afb41aa76885c8895c696b1bc0c97e9c314fd`.
- Preparation performed only a compiler syntax check. After GO, one compilation and one full traversal per geometry produced the fixed joint moments. No benchmark, small-size run, old result scoring, Monte Carlo, root solve or cloud work was performed.

```sh
clang++ -std=c++17 -O3 scripts/p337_regular_pair_joint_exact.cpp -o /private/tmp/p337_regular_pair_joint_exact
/private/tmp/p337_regular_pair_joint_exact 5 0 analysis/regular_pair_spatial_kernel.tsv results/p337-regular-pair-joint/axis.csv
/private/tmp/p337_regular_pair_joint_exact 4 3 analysis/regular_pair_spatial_kernel.tsv results/p337-regular-pair-joint/tilted.csv
```

These commands each enumerated exactly `2^25` configurations; the producer refuses existing output files. Scoring and importing the locked old root/D/U are the parent task's responsibility. The new source moments use the same exact finite occupation populations as earlier N25 work, not independent confirmation.
